import ROOT
from array import array
import numpy as np
from scipy.stats import chi2 as chi2_stats

# Load ROOT files and retrieve the histograms
S = 1

if S == 1:
	inputfile = ROOT.TFile("../Source.root", "open") # Signal
else:
	inputfile = ROOT.TFile('../Background.root',"open") # Background

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
		#if S == 1:
		x = ROOT.RooRealVar("mjj", "mjj", 0, 200)
		x.setRange("fitRange", 40, 160)
		#else:
		#	x = ROOT.RooRealVar("mjj", "mjj", 0, 200)
		#	x.setRange("fitRange", 20, 60)
		#dh = ROOT.RooDataHist("dh", "dh", [x], Import=hh)

		# Define the signal (gaussian)
		if S == 1:
			mean = ROOT.RooRealVar("mean", "mean", 90, 85, 110)
			sigma = ROOT.RooRealVar("sigma", "sigma", 10, 5, 25) 
		else:
			mean = ROOT.RooRealVar('mean','mean',35,40,60)
			sigma = ROOT.RooRealVar('sigma','sigma',5,0,20)
		
		gx = ROOT.RooGaussian('gauss','gauss',x,mean,sigma)

		# Define the background (Chebyshev polynomial)
		
		a = 10000
		a0 = ROOT.RooRealVar("a0", "a0", 0, -a, a)
		a1 = ROOT.RooRealVar("a1", "a1", 0, -a, a)
		#a2 = ROOT.RooRealVar("a2", "a2", 0, -a, a)
		#a3 = ROOT.RooRealVar("a3", "a3", 0, -a, a)
		#px = ROOT.RooChebychev("px", "px", x, ROOT.RooArgList(a0, a1, a2, a3))
		px = ROOT.RooChebychev("px", "px", x, ROOT.RooArgList(a0, a1))
		
		# Just trying another background (exponential)

		p0 = ROOT.RooRealVar('p0', 'p0', -0.01, -1, 0)
		p2 = ROOT.RooExponential('p2', 'p2', x, p0)
		
		# Construct composite model : signal + background 
		f = ROOT.RooRealVar("f", "f", 0.5, 0.0, 1.0)
		model = ROOT.RooAddPdf("model", "model", ROOT.RooArgList(gx, px), ROOT.RooArgList(f))
		#model = ROOT.RooAddPdf('model', 'model', ROOT.RooArgList(gx, p2), ROOT.RooArgList(f))
		#model = ROOT.RooAddPdf("model", "model", ROOT.RooArgList(gx, bkg), ROOT.RooArgList(f))
		
		# Convert histogram to RooDataHist
		dh = ROOT.RooDataHist(f"dh_{flavour}", f"dataset from histogram {hist_name}", ROOT.RooArgList(x), hh)

		# Perform a fit

		#print:(f'Total entries in histogram: {hh.Integral()}')
		#hh.Scale(1/hh.Integral()) # To normalize the data
		fitResult = model.fitTo(dh, ROOT.RooFit.Range('fitRange'), PrintLevel=-1, Save=True)
		fitResult.Print()
		Status = np.array([])
		if fitResult.status() != 0:
			print(f'Warning: fit did not converge for flavour {flavour}')
			Status = np.append(Status,flavour)
		# Plot model on data

		frame = x.frame(Title=f"Fit for Flavour {flavour}")
		dh.plotOn(frame)
		model.plotOn(frame, Components="gauss", LineStyle=2, LineColor=ROOT.kRed, Name='Signal')
		model.plotOn(frame, Components="px", LineStyle=2, LineColor=ROOT.kBlue, Name='Background')
		model.plotOn(frame,LineStyle=2,LineColor=ROOT.kMagenta, Name='TotalFit')

		chi2 = frame.chiSquare() # Computes chi2/ndof (reduced chi2)
		ndof = frame.GetNbinsX() - fitResult.floatParsFinal().getSize()
		p_value = 1 - chi2_stats.cdf(chi2*ndof,ndof)

		# Add legend
		leg = ROOT.TLegend(0.6, 0.7, 0.88, 0.88)
		leg.SetBorderSize(0)  # Removes the border to get a neater result
		leg.SetTextSize(0.03)  # Edits text length
		leg.SetFillStyle(0)  # Transparent background

		leg.AddEntry(frame.findObject("TotalFit"), "Total Fit", "l")
		leg.AddEntry(frame.findObject("Signal"), "Signal (Z)", "l")
		leg.AddEntry(frame.findObject("Background"), "Background (Chebyshev)", "l")
		#leg.AddEntry('', f'Chi^{{2}}/ndof = {chi2:.2f}','')

		if S == 1:
	
			c = ROOT.TCanvas(f"c_{flavour}", f"Fit for Flavour {flavour}", 800, 400)
			ROOT.gPad.SetLeftMargin(0.15)
			frame.GetYaxis().SetTitleOffset(1.4)
			frame.GetYaxis().SetTitleOffset(1.4)
			frame.Draw()
			leg.Draw()
			c.SaveAs(f"Source_fit_{flavour}.png")

			resid = frame.pullHist()  # Get pull histogram
			pull_frame = x.frame(Title=f"Pull distribution for Flavour {flavour}")
			pull_frame.addPlotable(resid, "P")

			c2 = ROOT.TCanvas(f"c2_{flavour}", f"Pulls for Flavour {flavour}", 800, 400)
			pull_frame.Draw()
			c2.SaveAs(f"Source_pull_distribution_flavour{flavour}.png")
		else:

			c = ROOT.TCanvas(f"c_{flavour}", f"Fit for Flavour {flavour}", 800, 400)
			ROOT.gPad.SetLeftMargin(0.15)
			frame.GetYaxis().SetTitleOffset(1.4)
			frame.GetYaxis().SetTitleOffset(1.4)
			frame.Draw()
			leg.Draw()
			c.SaveAs(f"Background_fit_{flavour}.png")

			resid = frame.pullHist()  # Get pull histogram
			pull_frame = x.frame(Title=f"Pull distribution for Flavour {flavour}")
			pull_frame.addPlotable(resid, "P")

			c2 = ROOT.TCanvas(f"c2_{flavour}", f"Pulls for Flavour {flavour}", 800, 400)
			pull_frame.Draw()
			c2.SaveAs(f"Background_pull_distribution_flavour{flavour}.png")


	
		# Save fit results in the text file
		output_file.write(f"\nFlavour {flavour} Results:\n")
		output_file.write(f"Chi-Squared/Ndof: {chi2:.3f}\n")
		output_file.write(f'ndof: {ndof:3f}\n')
		output_file.write(f'p_value : {p_value:3f}\n')
		output_file.write(f"a0 = {a0.getVal():.3f} ± {a0.getError():.3f}\n")
		output_file.write(f"a1 = {a1.getVal():.3f} ± {a1.getError():.3f}\n")
		#output_file.write(f"a2 = {a2.getVal():.3f} ± {a2.getError():.3f}\n")
		#output_file.write(f"a3 = {a3.getVal():.3f} ± {a3.getError():.3f}\n")
		output_file.write(f"f  = {f.getVal():.3f} ± {f.getError():.3f}\n")
		output_file.write(f"mean = {mean.getVal():.3f} ± {mean.getError():.3f}\n")
		output_file.write(f"sigma = {sigma.getVal():.3f} ± {sigma.getError():.3f}\n")
		output_file.write(f"Total entries in histogram: {hh.Integral():.0f}\n")
print('\n All fits completed')
