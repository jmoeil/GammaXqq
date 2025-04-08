import ROOT
#from array import array
import numpy as np
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

def plot_approx_deltaR(hist_eta_name,hist_phi_name,label_suffix):
	hist_eta = inputfile.Get(hist_eta_name)
	hist_phi = inputfile.Get(hist_phi_name)
	if not hist_eta or not hist_phi:
		print(f'Error : {hist_eta} or {hist_phi} not found.')
		return 
	nbins_eta = hist_eta.GetNbinsX()
	nbins_phi = hist_phi.GetNbinsX()
	
	if nbins_phi != nbins_eta:
		print(f'Error: {hist_eta} and {hist_phi} do not have matching binning.')
		return 
	x_vals = []
	y_vals = []
	
	for i in range(1, nbins_eta + 1):
		eta = hist_eta.GetBinCenter(i)
		phi = hist_phi.GetBinCenter(i)
		deltaR = np.sqrt(eta**2 + phi**2)
		
		eta_entries = hist_eta.GetBinContent(i)
		phi_entries = hist_phi.GetBinContent(i)
		entry = np.sqrt(eta_entries**2 + phi_entries**2)

		x_vals.append(deltaR)
		y_vals.append(entry)

	# Save in a CSV file
	csv_filename = f'Jet_DeltaR_{label_suffix}.csv'
	with open(csv_filename,mode='w',newline='') as f:
		writer = csv.writer(f)
		writer.writerow(["Delta R","Entries"])
		writer.writerows(zip(x_vals,y_vals))
		print(f'Saved data to {csv_filename}')

	plt.figure(figsize=(8,6))
	plt.scatter(x_vals,y_vals,color='purple',label=r'Approximate $\Delta R$',s=10)
	plt.xlabel(r'$\Delta R = \sqrt{\Delta\eta^2 + \Delta\phi^2}$')
	plt.ylabel('Approximate entries')
	plt.title(r'Approximate $\Delta R$ distribution')
	plt.grid()
	
	filename=f'Jet_deltaR_{label_suffix}.png'
	plt.savefig(filename,dpi=300)
	print(f'Saved {filename}')
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

	#plot_approx_deltaR(hist_eta_name,hist_phi_name,label_suffix)
	inputfile.Close()

'''
The following code was indented to sum the data and show the overlay. It doesn't work.
'''

## Define histogram pairs for overlaying
#hist_pairs = [
#    ("Jet_delta_eta_gjets;1", "Jet_delta_eta_zg;1"),
#    ("Jet_delta_phi_gjets;1", "Jet_delta_phi_zg;1"),
#    ("Jet_delta_pT_gjets;1", "Jet_delta_pT_zg;1"),
#]
#
## Open ROOT files
#background_file = ROOT.TFile("../Background.root", "open")
#source_file = ROOT.TFile("../Source.root", "open")
#
## Function to extract histogram data
#def extract_hist_data(hist):
#	bins = hist.GetNbinsX()
#	x_values = np.array([hist.GetBinCenter(i) for i in range(1, bins + 1)])
#	y_values = np.array([hist.GetBinContent(i) for i in range(1, bins + 1)])
#	return x_values, y_values
#
## Function to overlay histograms
#def sum_histograms(hist_bg_name, hist_src_name): hist_bg = backgroun.Get(hist_bg_name)
#	hist_src = source_file.Get(hist_src_name)
#
#	if not hist_bg or not hist_src:
#		print(f"Histogram {hist_bg_name} or {hist_src_name} not found!")
#		return
#
#	# Extract data
#	x_bg, y_bg = extract_hist_data(hist_bg)
#	x_src, y_src = extract_hist_data(hist_src)
#
#	# Sum histograms (bin-by-bin)
#	y_sum = y_bg + y_src
#
#	# Plot
#	plt.figure(figsize=(8, 6))
#	plt.scatter(x_bg, y_bg, color="red", label="Background (Gamma+jets)", s=10)
#	plt.scatter(x_src, y_src, color="blue", label="Signal (ZG)", s=10)
#	plt.scatter(x_bg, y_sum, color="purple", label="Summed Signal + Background", s=10)
#	
#	plt.xlabel(hist_bg.GetXaxis().GetTitle())
#	plt.ylabel("Entries")
#	plt.title(hist_bg_name.replace(";1", ""))
#	plt.grid()
#	plt.legend()
#
#	# Save overlayed histogram	
#	png_filename = hist_bg_name.replace(";1", "_comparison.png")
#	plt.savefig(png_filename, dpi=300)
#	print(f"Saved {png_filename}")
#	plt.close()
#
## Loop through histogram pairs
#for hist_bg_name, hist_src_name in hist_pairs:
#	sum_histograms(hist_bg_name, hist_src_name)
#
## Close files
#background_file.Close()
#source_file.Close()	 
