import ROOT
import sys
import looks

looks.looks()

the_file = ROOT.TFile(sys.argv[1],'read')
the_tree = the_file.Get('RunTree')

histo = ROOT.TH2D('histo','Electron in Water;E_{dep} (MeV);Layer;Events',30,1,0,8,0,8)

for event in the_tree:
    histo.Fill(event.L1_Edep,0)
    histo.Fill(event.L2_Edep,1)
    histo.Fill(event.L3_Edep,2)
    histo.Fill(event.L4_Edep,3)
    histo.Fill(event.L5_Edep,4)
    histo.Fill(event.L6_Edep,5)
    histo.Fill(event.L7_Edep,6)
    histo.Fill(event.L8_Edep,7)

c = ROOT.TCanvas()
c.SetLogz()
    
histo.Draw('lego')

raw_input('')
