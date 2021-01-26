# coding: utf-8


class Shape(object):

    # histogram naming template:
    # uproot does not support subdirectories yet...
    template = "{category}#{process}#{variable}#{systematic}"

    # template for extracting shapes for combine:
    ch_template_nom = "$BIN#$PROCESS#{variable}#nominal"
    ch_template_sys = "$BIN#$PROCESS#{variable}#$SYSTEMATIC"

    def __init__(self, name):
        self.extract(name)

    @property
    def name(self):
        return self._name

    @property
    def category(self):
        return self._properties["category"]

    @property
    def process(self):
        return self._properties["process"]

    @property
    def variable(self):
        return self._properties["variable"]

    @property
    def systematic(self):
        return self._properties["systematic"]

    def __str__(self):
        return f"Shape(category={self.category}, process={self.process}, variable={self.variable}, systematic={self.systematic})"

    def extract(self, name):
        assert len(name.split("#")) == 4
        self._name = name
        self._properties = {
            self.template.split("#")[i].strip("{}"): name.split("#")[i]
            for i in range(len(name.split("#")))
        }
