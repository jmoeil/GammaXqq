import ROOT
#from array import array
import numpy as np
import matplotlib
matplotlib.use("Agg")  # Use non-interactive backend (no display needed)
import matplotlib.pyplot as plt
from scipy.stats import chi2 as chi2_stats
import csv

def plot_hist(hist_name,x_label):
	hist = inputfile.Get(hist_name)
	if not hist:
		print(f'Histogram {hist_name} not found !')
		return
	
	# Convert ROOT histogram to numpy
	bins = hist.GetNbinsX()
	x_values = np.array([hist.GetBinCenter(i) for i in range(1, bins+1)])
	y_values = np.array([hist.GetBinContent(i) for i in range(1, bins+1)])
	
	# Plot with matplotlib
	plt.figure(figsize=(8,6))
	plt.scatter(x_values, y_values, color="b", label=hist_name, s=10)
	plt.xlabel(x_label)
	plt.ylabel("Entries")
	plt.title(hist_name)
	plt.grid()

	# Save as png
	png_filename = hist_name.replace(';1','') + '.png'
	plt.savefig(png_filename,dpi=300)
	print(f'Saved {png_filename}')
	plt.close()

# Process both Source and background
for S in [0,1]:
	if S == 1: 
		inputfile = ROOT.TFile('../Source.root','open')
		# Dictionnary linked histogram names to x-axis labels.
		hist_names = {
			"Jet_delta_eta_zg;1": r"$\Delta\eta$",
			"Jet_delta_phi_zg;1": r"$\Delta\phi$",
			"Jet_delta_pT_zg;1": r"$\Delta p_T$",
			'Jet_pT2pT1_zg;1':r"$\frac{p_{T2}}{p_{T1}}$",
			"Jet_delta_R_zg;1":r"$\Delta R = \sqrt{\Delta\eta^2 + \Delta\phi^2}$",
			"Jet_delta_eta_cut_pt2pt1_zg;1":r'$\Delta\eta$',
			"Jet_delta_phi_cut_pt2pt1_zg;1":r'$\Delta\phi$',
			"Jet_delta_pT_cut_pt2pt1_zg;1":r'$\Delta p_T$',
			"Jet_delta_R_cut_pt2pt1_zg;1":r"$\Delta R = \sqrt{\Delta\eta^2 + \Delta\phi^2}$"
		}
		hist_eta_name = 'Jet_delta_eta_zg;1'
		hist_phi_name = 'Jet_delta_phi_zg;1'
		label_suffix = 'zg'
	else:
		inputfile = ROOT.TFile('../Background.root','open')
		hist_names = {
			"Jet_delta_eta_gjets;1": r"$\Delta\eta$",
			"Jet_delta_phi_gjets;1": r"$\Delta\phi$",
			"Jet_delta_pT_gjets;1": r"$\Delta p_T$",
			'Jet_pT2pT1_gjets;1':r'$\frac{p_{T2}}{p_{T1}}$',
			"Jet_delta_R_gjets;1":r'$\Delta R = \sqrt{\Delta\eta^2 + \Delta\phi^2}$',
                        "Jet_delta_eta_cut_pt2pt1_gjets;1":r'$\Delta\eta$',
                        "Jet_delta_phi_cut_pt2pt1_gjets;1":r'$\Delta\phi$',
                        "Jet_delta_pT_cut_pt2pt1_gjets;1":r'$\Delta p_T$',
                        "Jet_delta_R_cut_pt2pt1_gjets;1":r"$\Delta R = \sqrt{\Delta\eta^2 + \Delta\phi^2}$"
		}
		hist_eta_name = 'Jet_delta_eta_gjets;1'
		hist_phi_name = 'Jet_delta_phi_gjets;1'
		label_suffix = 'gjets'
			
	for hist_name,x_label in hist_names.items():
		plot_hist(hist_name, x_label)

	inputfile.Close()
