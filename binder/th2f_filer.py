import ROOT

import coffea.hist as hist
import coffea
import sys

name = sys.argv[1]
template = coffea.util.load(name)

print(template.values()[()])

name = sys.argv[1].split(".coffea")[0]
outfile = ROOT.TFile(name + ".root", "recreate")
outfile.cd()

h1 = ROOT.TH2F("h1", "h1", 100, 200, 1000, 100, -8, -0.5)

print(h1.GetNbinsX())
print(h1.GetNbinsY())

for i in range(h1.GetNbinsX() + 2):
    for j in range(h1.GetNbinsY() + 2):
         h1.SetBinContent(i, j, template.values(overflow="allnan")[()][i][j])

h1.Write()
outfile.Close()
