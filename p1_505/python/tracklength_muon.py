import ROOT
import looks
looks.looks()

file_F_Pb = ROOT.TFile('../data/FTFP_Pb_muon.root')
file_F_PP = ROOT.TFile('../data/FTFP_PP_muon.root')
file_F_H2 = ROOT.TFile('../data/FTFP_H2O_muon.root')

tree_F_Pb = file_F_Pb.Get('RunTree')
tree_F_PP = file_F_PP.Get('RunTree')
tree_F_H2 = file_F_H2.Get('RunTree')

h_F_Pb = ROOT.TH1D('h_F_Pb',';Avg. Track Length/Layer;N events',100,40,200)
h_F_PP = ROOT.TH1D('h_F_PP',';Avg. Track Length/Layer;N events',100,40,200)
h_F_H2 = ROOT.TH1D('h_F_H2',';Avg. Track Length/Layer;N events',100,40,200)


def filler(the_tree,the_hist):
    for event in the_tree:
        tlsum = event.L1_Tlen + event.L2_Tlen + event.L3_Tlen + event.L4_Tlen + event.L5_Tlen + event.L6_Tlen + event.L7_Tlen + event.L8_Tlen
        tlavg = tlsum/8.0;
        the_hist.Fill(tlavg)
    the_hist.Scale(1./the_hist.Integral())
        
filler(tree_F_Pb,h_F_Pb)
filler(tree_F_PP,h_F_PP)
filler(tree_F_H2,h_F_H2)

h_F_Pb.SetLineColor(ROOT.kBlack)
h_F_PP.SetLineColor(ROOT.kRed)
h_F_H2.SetLineColor(ROOT.kAzure-2)

legend = ROOT.TLegend(.7,.7,.9,.9)
legend.SetHeader('2 GeV #mu^{-}')
legend.SetTextSize(0.042)
legend.SetTextFont(42)
legend.AddEntry(h_F_Pb,'Lead','l')
legend.AddEntry(h_F_PP,'Poly','l')
legend.AddEntry(h_F_H2,'Water','l')

c = ROOT.TCanvas()

h_F_PP.Draw()
h_F_Pb.Draw('same')
h_F_H2.Draw('same')

legend.Draw('same')

c.SaveAs('TL_muon.pdf')

raw_input('')
