# -*- coding: utf-8 -*-
from itertools import chain
from operator import itemgetter
import re


def getPGroot(config, process, group="default"):
    pgroups = config.aux["process_groups"][group]
    for p in chain((process,), map(itemgetter(0), process.walk_parent_processes())):
        if any(p == pg for pg in pgroups):
            return p


def NanoAOD_version(dataset):
    (ret,) = set(re.search(r"NanoAODv(\d+)", key).group(1) for key in dataset.keys)
    return int(ret)
