import ROOT
from array import array 
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2 as chi2_stats

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
	if S == 1 : 
		inputfile = ROOT.TFile('../Source.root','open')
		# Dictionnary linked histogram names to x-axis labels.
		hist_names = {
			"Jet_delta_eta_zg;1": r"$\Delta\eta$",
			"Jet_delta_phi_zg;1": r"$\Delta\phi$",
			"Jet_delta_pT_zg;1": r"$\Delta p_T$",
		}
	else:
		inputfile = ROOT.TFile('../Background.root','open')
		hist_names = {
			"Jet_delta_eta_gjets;1": r"$\Delta\eta$",
			"Jet_delta_phi_gjets;1": r"$\Delta\phi$",
			"Jet_delta_pT_gjets;1": r"$\Delta p_T$",
		}
			
	for hist_name,x_label in hist_names.items():
		plot_hist(hist_name, x_label)

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
