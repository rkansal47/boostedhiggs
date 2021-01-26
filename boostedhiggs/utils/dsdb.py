# coding: utf-8

from __future__ import absolute_import

import re
import json
import logging
from os import environ
from argparse import ArgumentParser, Namespace
from os.path import exists, commonprefix
from subprocess import check_output, CalledProcessError
from order import Dataset, DatasetInfo, UniqueObjectIndex
from operator import methodcaller

logger = logging.getLogger(__name__)


def dasquery(query, raw=False):
    args = [
        environ.get("DASGOCLIENT", "/cvmfs/cms.cern.ch/common/dasgoclient"),
        "-query=%s" % query,
        "-limit=0",
    ]
    logger.debug("query: %s", query)
    if not raw:
        args.append("-json")
    out = check_output(args)
    if raw:
        out = [l.decode("utf8") for l in map(methodcaller("strip"), out.splitlines()) if l]
    else:
        out = json.loads(out)
    return out


def accdict(iterable, onto=None):
    onto = onto or {}
    for other in filter(None, iterable):
        for key, value in other.items():
            if key in onto:
                onto[key] = onto[key] + value
            else:
                onto[key] = value
    return onto


class DSDB(dict):
    query = False

    def load(self, file, **kwargs):
        self.update(json.load(file, **kwargs))
        return self

    def dump(self, file, **kwargs):
        kwargs.setdefault("indent", 2)
        kwargs.setdefault("sort_keys", True)
        json.dump(self, file, **kwargs)

    def clean_missing(self):
        for key, value in self.items():
            if value is None:
                del self[key]

    def apply(self, obj):
        ob = getattr(obj, "datasets", obj)
        if isinstance(ob, UniqueObjectIndex):
            for o in ob.values():
                self.apply(o)
        elif isinstance(obj, Dataset):
            for o in obj.info.values():
                self.apply(o)
        elif isinstance(obj, DatasetInfo):
            info = accdict(map(self.parse, obj.keys))
            obj.n_files = info.pop("n_files", -1)
            obj.n_events = info.pop("n_events", -1)
            for key, value in info.items():
                obj.set_aux(key, value)
        else:
            logger.warning("cant process %r of type %s", obj, type(obj))

    def parse(self, k):
        v = self[k]
        if isinstance(v, list):
            return accdict(map(self.parse, v))
        elif isinstance(v, dict):
            v["keys"] = [k]
            if "lfns" in v and "lfn_prefix" in v:
                v = dict(v)
                lfn_prefix = v.pop("lfn_prefix", "")
                v["lfns"] = [lfn_prefix + lfn for lfn in v["lfns"]]
            return v

    def __missing__(self, key):
        ret = None
        if self.query:
            try:
                if "*" in key or "?" in key:
                    sub = dasquery("dataset=%s" % key, raw=True)
                    if (
                        len(sub) > 1
                        and len(set(re.sub(r"_ext[1-9]\d*", "", x, 1) for x in sub)) > 1
                    ):
                        logger.warning("skip ambiguous expansion for query: %s", key)
                    else:
                        ret = sub
                else:
                    ret = accdict(
                        dict(
                            lfns=[f["name"]],
                            n_events=f["nevents"],
                            n_bytes=f["size"],
                            n_files=1,
                        )
                        for item in dasquery("file dataset=%s" % key)
                        for f in item["file"]
                    )
                    if not ret:
                        logger.warning("dasquery returned no files for key: %s" % key)
                        lfn_prefix = ""
                    else:
                        lfn_prefix = commonprefix(ret["lfns"])
                    if "/" in lfn_prefix and len(ret["lfns"]) > 1:
                        lfn_prefix = lfn_prefix.rsplit("/", 1)[0] + "/"
                        len_prefix = len(lfn_prefix)
                        if len_prefix > 1:
                            ret["lfn_prefix"] = lfn_prefix
                            ret["lfns"] = [lfn[len_prefix:] for lfn in ret["lfns"]]
            except CalledProcessError as err:
                logger.warning(
                    "dasquery failed with code %d and output %s", err.returncode, err.output
                )
            else:
                logger.debug("add missing key: %s", key)
        else:
            logger.debug("missing key: %s", key)
        self[key] = ret
        return ret


def auto_apply(obj, name, file):
    """

    Usage in any campaign defining file:

    .. code-block:: python

        from utils.dsdb import auto_apply
        auto_apply(campaign, __name__, __file__)

    """
    if file.endswith(".pyc"):
        file = file[:-1]

    if name == "__main__":
        ap = ArgumentParser(description="dataset database tool for: %r" % obj)
        ap.add_argument("-l", "--load", action="store_true", help="load database")
        ap.add_argument("-u", "--update", action="store_true", help="update missing")
        ap.add_argument("-s", "--save", action="store_true", help="save database")
        ap.add_argument("-f", "--file", default=file, help="db file basename")
        ap.add_argument("-v", "--verbose", action="store_true", help="verbose output")
        args = ap.parse_args()

        from rich.logging import RichHandler

        logging.basicConfig(
            level="NOTSET",
            format="%(message)s",
            datefmt="[%X]",
            handlers=[RichHandler()],
        )
        logger.setLevel(level=logging.DEBUG if args.verbose else logging.INFO)
    else:
        args = Namespace(file=file, load=True, update=False, save=False)

    fn = args.file + ".json"
    db = DSDB()

    # load
    if args.load:
        if exists(fn):
            with open(fn, "r") as f:
                db.load(f)
            db.clean_missing()
            logger.info("loaded %d entries", len(db))
        else:
            logger.warn("no db file found")

    # apply & update
    if args.update:
        db.query = True
    db.apply(obj)

    if not args.save:
        if None in db.values():
            n_missing = list(db.values()).count(None)
            logger.error("number of missing keys: %i", n_missing)
    # save
    if args.save:
        db.clean_missing()
        with open(fn, "w") as f:
            db.dump(f)
        logger.info("saved %d entries", len(db))
