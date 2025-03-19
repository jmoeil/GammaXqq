import ROOT
from array import array
import numpy as np

# Load ROOT files and retrieve histograms

S = 0 # If 1, source simulation. If 0, gamma+jets simulation.

if S == 1:
	inputfile = ROOT.TFile('../output_pt.root','open') #Source
else:
	inputfile = ROOT.TFile('../GammaAndJets.root',"open") #Background

# Define histograms for leading and subleading jet pT
h_jet1_pt = ROOT.TH1F("h_jet1_pt", "Leading Jet p_{T};p_{T} (GeV);Entries", 50, 0, 500)
h_jet2_pt = ROOT.TH1F("h_jet2_pt", "Subleading Jet p_{T};p_{T} (GeV);Entries", 50, 0, 500)


for flavour in range(1,6):
	hist_name = f'jet_pt_partonflavour{flavour}_zg;1'
	hh = inputfile.Get(hist_name)
	
	# Check if histogram was plotted properly
	if not hh:
		print('Hstogram not found')
		inputfile.ls() #List contents of inputfile for debugging
		exit()
	
	# Convert histogram to RooDataHist
	pt = ROOT.RooRealVar("pt", "pt", 0, 500)
	dh = ROOT.RooDataHist(f"dh_{flavour}", f"dataset from histogram {hist_name}", ROOT.RooArgList(pt), hh)

	# Fill histograms for leading and subleading jets
	for bin in range(1, hh.GetNbinsX() + 1):
		bin_center = hh.GetBinCenter(bin)
		bin_content = hh.GetBinContent(bin)
		if bin_content > 0:
			if bin_center > 50:  # Example threshold for leading jet
				h_jet1_pt.Fill(bin_center, bin_content)
			else:
				h_jet2_pt.Fill(bin_center, bin_content)

	frame = pt.frame(Title=f'pT for Flavour {flavour}')
	dh.plotOn(frame)
	c = ROOT.TCanvas(f"c_{flavour}", f"pT for Flavour {flavour}", 800, 400)
	frame.Draw()
	
	if S == 1:
		c.SaveAs(f"Source_pt_plot{flavour}.png")
	else:
		c.SaveAs(f"Background_pt_plot{flavour}.png")
	# Plot the combined distributions for leading and subleading jets
	c_combined = ROOT.TCanvas("c_combined", "Jet pT Distribution", 800, 600)
	h_jet1_pt.SetLineColor(ROOT.kRed)
	h_jet2_pt.SetLineColor(ROOT.kBlue)
	h_jet1_pt.Draw("HIST")
	h_jet2_pt.Draw("HIST SAME")

	legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
	legend.AddEntry(h_jet1_pt, "Leading Jet pT", "l")
	legend.AddEntry(h_jet2_pt, "Subleading Jet pT", "l")
	legend.Draw()

	if S == 1:
	    c_combined.SaveAs(f"Source_jet_pt_distribution{flavour}.png")
	else:
	    c_combined.SaveAs(f"Background_jet_pt_distribution{flavour}.png")
