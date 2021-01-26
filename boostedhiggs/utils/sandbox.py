# coding: utf-8

import law
import os
from subprocess import Popen
from law.util import quote_cmd


class SimpleSandbox(law.Sandbox):
    def run(self, cmd, stdout=None, stderr=None):
        p = Popen(cmd, shell=True, env=self.env)

        try:
            p.wait()
        except KeyboardInterrupt:
            pass
        p.wait()

        return (p.returncode,) + p.communicate()

    @property
    def env(self):
        env = dict(os.environ)
        env.update(self._get_env())
        return env


class SGSandbox(SimpleSandbox):
    sandbox_type = "sg"

    def cmd(self, proxy_cmd):
        return law.util.quote_cmd(["sg", self.name, "-c", proxy_cmd.build()])


class OpenportsTask(law.SandboxTask):
    sandbox = "sg::openports"
