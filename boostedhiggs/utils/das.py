# -*- coding: utf-8 -*-

import re
import logging

import numpy as np
from utils.dsdb import dasquery


def get_name(str):
    str = str.split("_dipoleRecoil")[0].split("_Tune")[0].split("_amcatnlo")[0].split("_13TeV")[0]
    str = str.replace("/", "")
    return str


def has_expr(dataset, expression="ext"):
    name = get_das_name(dataset)
    return expression in name


def get_dataset(dataset):
    ds = dataset["dataset"]
    if len(ds) == 1:
        return ds[0]
    else:
        raise NotImplementedError


def get_das_name(dataset):
    ds = get_dataset(dataset)
    return ds["name"]


def get_creation_date(dataset):
    if dataset is None:
        return -1
    else:
        ds = get_dataset(dataset)
        return ds["creation_date"]


def get_gensim(dataset, identifier="GEN-SIM"):
    name = get_das_name(dataset)
    while identifier not in name:
        query = dasquery("parent dataset=%s" % name)
        if len(query) > 0:
            dataset = query[0]["parent"][0]  # % FIXME! - dataset can have more than one parent
            name = dataset["name"]
        else:
            return None
    return dataset


def dataset_has_name(dataset, name="new_pmx", match_case=True):
    das_name = get_das_name(dataset)
    if not match_case:
        name, das_name = name.lower(), das_name.lower()
    return name in das_name


def get_datasets_by_name(datasets, name="new_pmx", parent="", match_case=True):
    if len(parent) > 0:
        out = []
        for dataset in datasets:
            parent_dataset = get_gensim(dataset, parent)
            if parent_dataset is None:
                logging.info(f"dropping {get_das_name(dataset)} because it has no {parent}")
                continue
            a, b = name, parent_dataset["name"]
            if not match_case:
                a, b = a.lower(), b.lower()
            if a in b:
                out.append(dataset)
    else:
        out = [dataset for dataset in datasets if dataset_has_name(dataset, name, match_case)]
    return out


def get_version(dataset, regex="v(\d+)/"):
    if not isinstance(dataset, str):
        dataset = get_das_name(dataset)

    version = re.search(regex, dataset)
    if version:
        return int(version.group(1))
    else:
        return 0


def sort_datasets_by_expressions(datasets, expressions):
    idx = np.argsort(-np.array(expressions).astype(np.int))
    return list(np.array(datasets)[idx])


def get_highest_version(datasets):
    sorted_datasets = sort_by_version(datasets, regex=".+ver(\d+).+")
    return sorted_datasets[0]


def sort_by_version(datasets, regex="v(\d+)/"):
    versions = []
    for dataset in datasets:
        versions.append(get_version(dataset, regex=regex))
    return sort_datasets_by_expressions(datasets, versions)


def sort_by_ps_weights(datasets):
    ps_weights = []
    for dataset in datasets:
        ps_weights.append(has_expr(dataset, expression="PSWeights"))
    return sort_datasets_by_expressions(datasets, ps_weights)


def sort_by_ext(datasets):
    ext = []
    for dataset in datasets:
        ext.append(has_expr(dataset, expression="ext"))
    return sort_datasets_by_expressions(datasets, ext)


def sort_by_new_pmx(datasets):
    new_pmx = []
    for dataset in datasets:
        new_pmx.append(has_expr(dataset, expression="new_pmx"))
    return sort_datasets_by_expressions(datasets, new_pmx)


def get_latest(datasets):
    out = []

    datasets = sort_by_ps_weights(datasets)
    datasets = sort_by_ext(datasets)
    datasets = sort_by_version(datasets)
    datasets = sort_by_new_pmx(datasets)

    gen_sim = [get_gensim(dataset) for dataset in datasets]

    if any(_ is None for _ in gen_sim):
        for dataset in datasets:
            dataset["dataset"][0]["name"] += "_nogen"
            out.append(dataset)
    else:
        creation_dates = []
        for _ in gen_sim:
            _dataset = dasquery(_["name"])[0]
            creation_dates.append(get_creation_date(_dataset))

        sorted_idx = np.argsort(-np.array(creation_dates).astype(np.int))
        burned = set([])
        for idx, creation_date in zip(sorted_idx, creation_dates):
            if creation_date == -1 or creation_date in burned:
                continue
            out.append(datasets[idx])
            burned.add(creation_date)
    return out


def clear_datasets(datasets):
    out = []
    keywords = [
        "percentMaterial",
        "FlatPU",
        "BSandPUSummer16",
        "pilot",
        "ForMUOVal",
        "MUOTrackFix",
        "PU2017RECOPF",
        "PU2017RECOSIM",
    ]
    for dataset in datasets:
        name = get_das_name(dataset)
        if any(keyword in name for keyword in keywords):
            continue
        out.append(dataset)
    return out


def same_dataset_different_version(datasets):
    unified_das_names = []
    for dataset in datasets:
        das_name = get_das_name(dataset)
        unified_das_name = re.sub("v(\d+)/", "", das_name)
        unified_das_names.append(unified_das_name)
    return len(set(unified_das_names)) == 1


def clear_pu(datasets):
    out = []
    for dataset in datasets:
        if (
            len(get_datasets_by_name([dataset], "new_pmx")) > 0
            or len(get_datasets_by_name([dataset], "new_pmx", parent="/MINIAODSIM")) > 0
            or len(get_datasets_by_name([dataset], "PU2017", parent="/MINIAODSIM")) > 0
            or len(get_datasets_by_name([dataset], "PU2017", parent="/AODSIM")) > 0
        ):
            out.append(dataset)
        else:
            logging.info(f"pu cleaning dropped:, {get_das_name(dataset)}")
    return out


def get_datasets(das_name, year):
    datasets_out = []
    datasets = dasquery(das_name)
    datasets = clear_datasets(datasets)
    if len(datasets) == 0:
        logging.info(f"No datasets avalable for: {das_name}")
    elif len(datasets) == 1:
        datasets_out = [datasets[0]]
    else:
        logging.info("== Make decision to choose out of:")
        logging.info("\n".join([get_das_name(ds) for ds in datasets]))
        # if "TuneCP5_PSweights" in das_name:
        #     from IPython import embed
        #
        #     embed()
        # has_psweights = get_datasets_by_name(datasets, "PSWeights", match_case=False)
        # if len(has_psweights) > 0:
        #     datasets = has_psweights

        datasets_out = get_latest(datasets)
        if same_dataset_different_version(datasets_out):
            datasets_out = datasets_out[:1]
        else:
            logging.info("NONTRIVIAL")
        logging.info("== chosen:")
        logging.info("\n".join([get_das_name(ds) for ds in datasets_out]))
    dataset = {"keys": [get_das_name(ds) for ds in datasets_out]}
    return dataset


def split_das_name(das_name):
    values = das_name.split("/")
    if len(values) == 4:
        _, name, tag, sim = values
    else:
        print(das_name)
        print(values)
        raise NotImplementedError
    return name, tag, sim


def merge_das_name(*args):
    return "/".join(["", *args])


def get_20xx(name, year, ppd, sim="nano"):
    tag = ppd[year][sim]["name"] + "*" + ppd[year][sim]["gt"] + "*"
    das_name = merge_das_name(name, tag, sim.upper() + "AODSIM")
    datasets = get_datasets(das_name, year)
    return datasets


def change_campaign(datasets, year, ppd, sim="nano"):
    names = []
    xs = {}
    changed_datasets = []
    for dataset in datasets:
        name, tag, _ = split_das_name(dataset["keys"][0])
        if name not in names:
            names.append(name)
            xs[name] = dataset["xs"]
    for name in names:
        print(name)
        dataset = get_20xx(name, year, ppd, sim)
        dataset["xs"] = xs[name]
        changed_datasets.append(dataset)

    return changed_datasets


def change_campaign_data(datasets, year, ppd, sim="nano"):
    changed_datasets = []
    for dataset in datasets:
        name, tag, _ = split_das_name(dataset["keys"][0])
        run = tag.split("-")[0]
        das_name = merge_das_name(
            name, "-".join([run + "*", ppd[year][sim]["name"] + "*"]), sim.upper() + "AOD"
        )
        if all(_ in das_name for _ in ["MINIAOD", "Run2018D"]):
            if any(_ in das_name for _ in ["EGamma", "SingleMuon"]):
                das_name = merge_das_name(
                    name, "-".join([run + "*", "22Jan2019" + "*"]), sim.upper() + "AOD"
                )
            else:
                das_name = merge_das_name(
                    name, "-".join([run + "*", "Prompt" + "*"]), sim.upper() + "AOD"
                )
        _ = dasquery(das_name)
        _dataset = get_highest_version(_)
        if _dataset is None:
            raise NotImplementedError
        dataset["keys"] = [get_das_name(_dataset)]
        changed_datasets.append(dataset)
    return changed_datasets


simple_channels = {
    "SingleElectron": ["e"],
    "EGamma": ["e", "ee"],
    "SingleMuon": ["mu"],
    "DoubleEG": ["ee"],
    "DoubleMuon": ["mumu"],
    "MuonEG": ["emu"],
}


def to_dict(dataset, is_data=False):
    aux = {}
    xs = None
    keys = dataset["keys"]
    if "xs" in dataset:
        xs = dataset["xs"]
    if "aux" in dataset:
        aux = dataset["aux"]

    if len(keys) == 0:
        return {"name": "NON_EXISTENT"}
    das_name = keys[0]
    if dataset["is_data"]:
        physics_name, run, format = das_name.split("/")[1:]
        channels = simple_channels[physics_name]
        run_name = re.search("Run\d+(\w)[-_]", run).group(1)
        aux.update({"run": run_name, "channels": channels})
        name = "_".join(["data", run_name, "+".join(channels)])
    else:
        name = get_name(das_name)
    return {
        "name": name,
        "is_data": dataset["is_data"],
        "aux": aux,
        "keys": keys,
        "misc": {"xs": xs},
    }


def dictify_datasets(datasets):
    dictified = []
    for dataset in datasets:
        dictified.append(to_dict(dataset))
    return dictified
