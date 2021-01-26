# -*- coding: utf-8 -*-

import re

from markdown import markdown
from bs4 import BeautifulSoup


class MDPyTree:
    def __init__(self, mother, level=0, el=""):
        self.mother = mother
        self.el = el
        self.level = level
        self.children = []
        self.content = []

    def add_child(self, child):
        self.children.append(child)

    def add_content(self, content):
        self.content.append(MDPyTree(self.mother, self.level, content))

    def get_children(self):
        return self.children

    def get_mother(self):
        return self.mother

    def tell_mother(self):
        self.mother.add_child(self)

    @staticmethod
    def fromMarkdown(md, *args, **kwargs):
        """
        Creates abstraction using path to file

        :param str path: path to markdown file
        :return: MDPyTree object
        """
        return MDPyTree.fromHTML(markdown(md, *args, **kwargs))

    @staticmethod
    def fromHTML(html, *args, **kwargs):
        """
        Creates abstraction using HTML

        :param str html: HTML
        :return: MDPyTree object
        """
        source = BeautifulSoup(html, "html.parser", *args, **kwargs)

        tree = MDPyTree(mother=None, level=0)
        for el in source:
            if el.name is None:
                continue
            level = get_heading_level(el.name)
            if level:
                if level <= tree.level:
                    moved_tree = tree
                    for i in range(1 + tree.level - level):
                        moved_tree = moved_tree.get_mother()
                    tree = moved_tree

                new_tree = MDPyTree(mother=tree, level=level, el=el)
                new_tree.tell_mother()
                tree = new_tree
            else:
                tree.add_content(el)

        for i in range(tree.level):
            tree = tree.get_mother()

        return tree


def join_contents(l):
    a = ""
    for s in l:
        a += str(s)
    return a


def get_heading_level(s):
    for i in range(10):
        if s == "h" + str(i):
            return i


def strip_remnants(value):
    kilo = 1000
    mio = 1000 * kilo

    if not isinstance(value, str):
        return value
    elif str == "":
        return None
    else:
        values = re.findall("(?:<sup>|<sub>)(.+)(?:<\/sup>|<\/sub>)", value)
        if len(values) > 0:
            value = values[0]
        values = re.findall("^(\d+\.*\d*[MK]*).*", value)
        if len(values) > 0:
            value = values[0]
            if value.endswith("K"):
                value = float(value[:-1]) * kilo
            elif value.endswith("M"):
                value = float(value[:-1]) * mio
            else:
                value = float(value)
        return value


def parse_table(table_string, is_data=False):
    erases = [
        "\n<summary>Click me to collapse/expand.</summary>\n\n",
        "\n\n\\* appears to be fully inclusive\n\n",
    ]

    datasets = []
    ids = []
    for erase in erases:
        table_string = table_string.replace(erase, "")
    table = table_string.split("\n")
    for id in table[0].split("|"):
        id = id.strip()
        if id == "":
            continue
        if id == "XS [pb]":
            id = "XS"
        ids.append(id)
    num_ids = len(ids)
    for dataset in table[2:]:
        values = [a.strip() for a in dataset.split("|")[1:][:num_ids]]
        if all(value == "" for value in values):
            continue
        if any(value == "" for value in values[: num_ids - 1]):
            print("Problem in parsing:", values)
            continue
        if "~~" in values[0]:
            print("Discarded:", values)
            continue
        values = [strip_remnants(value) for value in values]
        _dataset = {id: group for (id, group) in zip(ids, values)}
        _dataset["keys"] = [_dataset["DAS name" if is_data else "Sample"]]
        if not is_data:
            if "BR" in _dataset:
                _dataset["xs"] = _dataset["BR"]
            else:
                _dataset["xs"] = _dataset["XS"]
        datasets.append(_dataset)
    return datasets
