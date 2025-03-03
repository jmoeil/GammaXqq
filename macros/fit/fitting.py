import ROOT
from array import array
import numpy as np

# Load ROOT files and retrieve the histograms
inputfile = ROOT.TFile("../output_100to200and200.root", "open")

# Open a text file to store results
output_filename = "fit_results.txt"
with open(output_filename, 'w') as output_file:
	# Loop over flavour 1 to 5
	for flavour in range(1,6):
		hist_name = f'mjj_partonflavor{flavour}_zg;1'
		hh = inputfile.Get(hist_name)
		
		# Check if the histogram was loaded properly
		if not hh:
			print(f"Error: Histogram {hist_name} not found in file!")
			inputfile.ls()  # List contents of the ROOT file for debugging
			exit()

		# Define the variable
		x = ROOT.RooRealVar("mjj", "mjj", 0, 1000)
		x.setRange("fitRange", 60, 350)
		#dh = ROOT.RooDataHist("dh", "dh", [x], Import=hh)

		# Define the signal (gaussian)
		mean = ROOT.RooRealVar("mean", "mean", 91.2, 85, 110)
		sigma = ROOT.RooRealVar("sigma", "sigma", 10, 5, 25) 
		gx = ROOT.RooGaussian("gauss", "gauss", x, mean, sigma)


		# Define the background (Chebyshev polynomial)
		a = 4
		a0 = ROOT.RooRealVar("a0", "a0", 0, -a, a)
		a1 = ROOT.RooRealVar("a1", "a1", 0, -a, a)
		a2 = ROOT.RooRealVar("a2", "a2", 0, -a, a)
		a3 = ROOT.RooRealVar("a3", "a3", 0, -a, a)
		px = ROOT.RooChebychev("px", "px", x, ROOT.RooArgList(a0, a1, a2, a3))

		# Just trying another background (exponential)

		#p0 = ROOT.RooRealVar('p0', 'p0', -0.01, -1, 0)
		#p2 = ROOT.RooExponential('p2', 'p2', x, p0)


		# Construct composite model : signal + background 
		f = ROOT.RooRealVar("f", "f", 0.5, 0.0, 1.0)
		model = ROOT.RooAddPdf("model", "model", ROOT.RooArgList(gx, px), ROOT.RooArgList(f))
		#model = ROOT.RooAddPdf('model', 'model', ROOT.RooArgList(gx, p2), ROOT.RooArgList(f))

		# Convert histogram to RooDataHist
		dh = ROOT.RooDataHist(f"dh_{flavour}", f"dataset from histogram {hist_name}", ROOT.RooArgList(x), hh)

		# Perform a fit

		#print:(f'Total entries in histogram: {hh.Integral()}')
		#hh.Scale(1/hh.Integral()) # To normalize the data
		fitResult = model.fitTo(dh, ROOT.RooFit.Range('fitRange'), PrintLevel=-1, Save=True)
		fitResult.Print()

		# Plot model on data

		frame = x.frame(Title=f"Fit for Flavour {flavour}")
		dh.plotOn(frame)
		model.plotOn(frame, Components='gx', LineStyle=2, LineColor=ROOT.kRed)
		model.plotOn(frame, Components="px", LineStyle=2, LineColor=ROOT.kBlue)

		# Add legend
		leg = ROOT.TLegend(0.6, 0.7, 0.88, 0.88)
		leg.SetBorderSize(0)  # Removes the border to get a neater result
		leg.SetTextSize(0.03)  # Edits text length
		leg.SetFillStyle(0)  # Transparent background

		leg.AddEntry(frame.findObject("mij"), "Total Fit", "l")
		leg.AddEntry(frame.findObject("gx_Norm[mjj]"), "Signal (Z)", "l")
		leg.AddEntry(frame.findObject("px_Norm[mjj]"), "Bruit de fond (Chebyshev)", "l")

		c = ROOT.TCanvas(f"c_{flavour}", f"Fit for Flavour {flavour}", 800, 400)
		ROOT.gPad.SetLeftMargin(0.15)
		frame.GetYaxis().SetTitleOffset(1.4)
		frame.GetYaxis().SetTitleOffset(1.4)
		frame.Draw()
		c.SaveAs(f"fit_{flavour}.png")

		resid = frame.pullHist()  # Get pull histogram
		pull_frame = x.frame(Title=f"Pull distribution for Flavour {flavour}")
		pull_frame.addPlotable(resid, "P")

		c2 = ROOT.TCanvas(f"c2_{flavour}", f"Pulls for Flavour {flavour}", 800, 400)
		pull_frame.Draw()
		c2.SaveAs(f"pull_distribution_flavour{flavour}.png")

		chi2 = frame.chiSquare()
		
		# Save fit results in the text file
		output_file.write(f"\nFlavour {flavour} Results:\n")
		output_file.write(f"Chi-Squared/Ndof: {chi2:.3f}\n")
		output_file.write(f"a0 = {a0.getVal():.3f} ± {a0.getError():.3f}\n")
		output_file.write(f"a1 = {a1.getVal():.3f} ± {a1.getError():.3f}\n")
		output_file.write(f"a2 = {a2.getVal():.3f} ± {a2.getError():.3f}\n")
		output_file.write(f"a3 = {a3.getVal():.3f} ± {a3.getError():.3f}\n")
		output_file.write(f"f  = {f.getVal():.3f} ± {f.getError():.3f}\n")
		output_file.write(f"mean = {mean.getVal():.3f} ± {mean.getError():.3f}\n")
		output_file.write(f"sigma = {sigma.getVal():.3f} ± {sigma.getError():.3f}\n")
		output_file.write(f"Total entries in histogram: {hh.Integral():.0f}\n")
print('\n All fits completed')
