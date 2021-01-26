#!/usr/bin/env python
# coding: utf-8

import grpc
from grpc._cython.cygrpc import CompressionAlgorithm, CompressionLevel
import tensorflow as tf
from tensorflow_serving.apis import predict_pb2
from tensorflow_serving.apis import prediction_service_pb2_grpc
from tensorflow_serving.config.model_server_config_pb2 import (
    ModelConfig,
    ModelConfigList,
    ModelServerConfig,
)
from random import randrange
from tempfile import TemporaryDirectory
from pathlib import Path
from subprocess import Popen, PIPE, DEVNULL, TimeoutExpired
from tools.evil import pin
from weakref import finalize
from time import sleep
from socket import create_connection, getfqdn
from functools import cached_property
from os import environ
from shlex import quote
from urllib.parse import urlparse


class Client:
    def __init__(self, address, **kwargs):
        if isinstance(address, BaseServer):
            self.server = address
            address = address.address
        pin(locals())

    @cached_property
    def channel(self):
        return grpc.insecure_channel(self.address, **self.kwargs)

    def __getstate__(self):
        state = dict(self.__dict__)
        state.pop("channel", None)
        server = state.pop("server", None)
        if server:
            state["address"] = server.host_address
        return state

    @property
    def stub(self):
        return prediction_service_pb2_grpc.PredictionServiceStub(self.channel)

    def __call__(
        self,
        *args,
        inputs=None,
        outputs=None,
        model="default",
        signature="serving_default",
    ):
        assert inputs or args
        request = predict_pb2.PredictRequest()
        request.model_spec.name = model
        request.model_spec.signature_name = signature
        if args:
            inputs = {"input_%d" % i: v for i, v in enumerate(args, start=1)}
        for key, value in inputs.items():
            request.inputs[key].CopyFrom(tf.make_tensor_proto(value))
        if outputs:
            request.output_filter.extend(outputs())
        try:
            result = self.stub.Predict(request)
        except Exception as e:
            raise RuntimeError("Predict failed") from e
        return {key: tf.make_ndarray(tproto) for key, tproto in result.outputs.items()}


class BaseServer:
    exe = "tensorflow_model_server"

    def __init__(self, command, /, auto_stop=True, **kwargs):
        self.command = command
        self.process = Popen(command, **kwargs)
        if auto_stop:
            finalize(self, self._stop, self.process)

    @property
    def running(self):
        return self.process.returncode is None

    def stop(self, wait=False, kill=False):
        if self.running:
            self.process.terminate()
            if wait is not False:
                try:
                    self.process.wait(None if wait is True else wait)
                except TimeoutError:
                    if self.running and kill:
                        self.process.kill()
                        self.process.wait()
        return self.process.returncode

    @staticmethod
    def _stop(proc):
        if proc.returncode is None:
            proc.terminate()

    def client(self, **kwargs):
        return Client(self, **kwargs)

    @property
    def unix_path(self):
        return self.kwargs.get("grpc_socket_path", None)

    @property
    def unix_address(self):
        u = self.unix_path
        return "unix://%s" % u if u else None

    @property
    def local_address(self):
        return "localhost:%d" % self.port

    @property
    def address(self):
        return self.unix_address or self.host_address


class Server(BaseServer):

    def __init__(self, /, port=0, auto_stop=True, env=None, cuda_exe=None, **kwargs):
        self.port = kwargs["port"] = port or self.random_port()
        self.kwargs = kwargs
        exe = kwargs.pop("exe", self.exe)
        if cuda_exe and environ.get("CUDA_VISIBLE_DEVICES", None) != "-1":
            exe = cuda_exe
        super().__init__(
            [exe] + ["--%s=%s" % (k, v) for k, v in kwargs.items()],
            auto_stop=auto_stop,
            env=env,
        )

    @classmethod
    def relieable(cls, /, delay=0.5, retry=0, **kwargs):
        attempts = retry + 1
        for i in range(attempts):
            s = cls(**kwargs)
            while s.running:
                sleep(delay)
                try:
                    create_connection((None, s.port)).close()
                except ConnectionRefusedError:
                    if delay < 5:
                        delay += 0.1
                    continue
                else:
                    return s
        raise RuntimeError("failed %d attemps to start %r" % (attempts, s))

    @classmethod
    def auto(cls, model_path, **kwargs):
        if not isinstance(model_path, dict):
            model_path = dict(default=model_path)
        tmp = TemporaryDirectory(prefix="tf_server_ahoc")
        tmp_path = Path(tmp.name)
        dirs = {}
        for model, path in model_path.items():
            path = Path(path)
            if not path.is_dir():
                raise RuntimeError("%s: %s is not a dir" % (model, path))
            b = tmp_path.joinpath(model)
            dirs[model] = b.absolute().as_posix()
            if path.joinpath("saved_model.pb").is_file():
                b.mkdir()
                b = b.joinpath("1")
            b.symlink_to(path.absolute(), True)
        del model, path, b
        if len(dirs) > 1:
            mc = tmp_path.joinpath("model_config.txt.pb")
            kwargs["model_config_file"] = mc.absolute()
            with mc.open("w") as f:
                f.write(
                    str(
                        ModelServerConfig(
                            model_config_list=ModelConfigList(
                                config=[
                                    ModelConfig(
                                        name=key, base_path=value, model_platform="tensorflow"
                                    )
                                    for key, value in dirs.items()
                                ]
                            )
                        )
                    )
                )
        else:
            ((kwargs["model_name"], kwargs["model_base_path"]),) = dirs.items()
        kwargs.setdefault("grpc_socket_path", tmp_path.joinpath("sock").absolute())
        self = cls.relieable(**kwargs)
        pin(locals(), cls)
        return self

    @classmethod
    def random_port(self, randint=randrange):
        """ randint is expected to behave like numpy.random.randint or random.randrage """
        return randint(
            *map(int, open("/proc/sys/net/ipv4/ip_local_port_range", "r").read().split())
        )

    @property
    def host_address(self):
        return "%s:%d" % (getfqdn(), self.port)


class RemoteServer(BaseServer):
    def __init__(self, command, /, auto_stop=True, env=None, timeout=None, silent=False):
        self.kwargs = {}
        super().__init__(
            command,
            auto_stop=auto_stop,
            env=env,
            stdout=PIPE,
            stderr=DEVNULL if silent else None
        )

        # nasty blocking read
        self.host_address = self.process.stdout.readline().decode("utf8").strip()

        assert self.process.poll() is None
        assert self.host_address
        parsed = urlparse(f"//{self.host_address}")
        self.host = parsed.hostname
        self.port = parsed.port
        if timeout is not False:
            create_connection((self.host, self.port), timeout=timeout).close()

    @classmethod
    def cmdline(self, model_path, **kwargs):
        cmd = [__file__]
        for k, v in kwargs.items():
            cmd += ["--" + k.replace("_", "-"), str(v)]
        if isinstance(model_path, dict):
            assert not any("=" in k for k in model_path.keys())
            cmd += ["%s=%s" % kv for kv in model_path.items()]
        else:
            cmd.append(model_path)
        return cmd


# environment specific helper
def autoClient(
    model_path,
    /,
    remote=False,
    sg="openports",
    ram=10000,
    vram=10000,
    sargs=dict(
        cuda_exe="/net/scratch/BFischer/tfsrv/gpu/tms_gpu"
    ),
    cargs=dict(
        options=[
            ("grpc.max_send_message_length", 1 << 29),
            ("grpc.max_receive_message_length", 1 << 29),
            ("grpc.default_compression_algorithm", CompressionAlgorithm.gzip),
            ("grpc.default_compression_level", CompressionLevel.low),
        ]
    )
):
    if remote:
        cmd = RemoteServer.cmdline(model_path, **sargs)
        cmd = ["sg", sg, "-c", " ".join(map(quote, cmd))]
        cmd = ["submit", "-f", "-m", str(ram), "-M", str(vram)] + cmd
        srv = RemoteServer(cmd)
    else:
        srv = Server.auto(model_path, **sargs)
    return srv.client(**cargs)


if __name__ == "__main__":
    from argparse import ArgumentParser
    from sys import stderr
    import signal
    import re

    ap = ArgumentParser()
    ap.add_argument("-p", "--port", default=0, type=int)
    ap.add_argument("-e", "--exe", default=None, type=Path)
    ap.add_argument("-c", "--cuda-exe", default=None, type=Path)
    ap.add_argument("-r", "--retry", default=3, type=int)
    ap.add_argument("-d", "--delay", default=0.5, type=float)
    ap.add_argument("model", nargs="+")

    args = ap.parse_args()

    assert 0 <= args.port < 0xFFFF
    assert not args.exe or args.exe.exists()
    assert not args.cuda_exe or args.cuda_exe.exists()
    assert 0 <= args.retry
    assert 0 <= args.delay

    model = {}
    for m in map(re.compile(r"(?:(\w+)=)?(.*)").match, args.model):
        k, v = m.groups()
        if not v:
            if k:
                print(f"key w/o value interpreted as value: {k}", file=stderr)
                k, v = None, k
            else:
                continue
        if not k:
            k = "default"
        if k in model:
            print(f"redefinition of key: {k}", file=stderr)
        model[k] = v
    stderr.flush()

    def _exit(sig_no, sig_frame):
        exit(2)

    signal.signal(signal.SIGTERM, _exit)

    srv = Server.auto(
        model,
        port=args.port,
        exe=args.exe or Server.exe,
        cuda_exe=args.cuda_exe or None,
        retry=args.retry,
        delay=args.delay,
    )

    print(srv.host_address, flush=True)

    while True:
        sleep(1234)
