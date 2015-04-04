import sys
import ROOT

the_file = ROOT.TFile(sys.argv[1],'read')
the_tree = the_file.Get('RunTree')

for event in the_tree:
    print event.L1_Edep + event.L2_Edep + event.L3_Edep + event.L4_Edep
