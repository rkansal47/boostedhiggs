from uproot.write.TKey import TKey
from uproot.write.sink.cursor import Cursor
from uproot.write.objects.TTree import TTree
from uproot.write.TFile import TFileUpdate
from uproot.write.objects.util import Util

import uproot_methods.convert


def bulk_write(self, content):
    assert isinstance(self, TFileUpdate)

    self.util = Util()

    cursor = Cursor(self._fSeekFree)
    for where, what in content:
        where, cycle = self._normalizewhere(where)

        isTTree = what.__class__.__name__ in ("newtree", "TTree")
        assert not isTTree  # prevent TTree writing, otherwise migth invoke nasty magic
        if not isTTree:
            what = uproot_methods.convert.towriteable(what)
        elif what.__class__.__name__ == "newtree":
            what = TTree(where, what, self)

        newkey = TKey(
            fClassName=what._fClassName,
            fName=where,
            fTitle=what._fTitle,
            fObjlen=0,
            fSeekKey=cursor.index,
            fSeekPdir=self._fBEGIN,
            fCycle=cycle if cycle is not None else self._rootdir.newcycle(where),
        )
        if isTTree:
            # Need to (re)attach the cycle number to allow getitem to access writable TTree
            tree_where = where + b";" + str(newkey.fCycle).encode("utf-8")
            self._treedict[tree_where] = what

        newkeycursor = Cursor(newkey.fSeekKey)
        newkey.write(cursor, self._sink)
        what._write(self, cursor, where, self.compression, newkey, newkeycursor, self.util)

        # prevent overwrite (migth be excessive and actually work fine)
        assert (newkey.fName, newkey.fCycle) not in self._rootdir.keys

        self._rootdir.headkey.fObjlen += newkey.fKeylen
        self._rootdir.keys[(newkey.fName, newkey.fCycle)] = newkey

    # write (root) TDirectory
    self._rootdir.fNbytesKeys = self._rootdir._nbyteskeys()
    while self._rootdir.fNbytesKeys > self._rootdir.allocationbytes:
        self._rootdir.allocationbytes *= self._rootdir.growfactor

    self._rootdir.writekeys(cursor)

    self._expandfile(cursor)
