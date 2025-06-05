import ROOT
import matplotlib
matplotlib.use("Agg")  # Force un backend non-interactif (pas de GUI requise)
import matplotlib.pyplot as plt
from array import array
import numpy as np
import os

# --- Import ROOT files ---

file_bg = ROOT.TFile("../Background.root")
file_sg = ROOT.TFile('../Source.root')

cuts = {
	"nocut":{
		"base":"",
		"label":"No cut"	
	},
	"pt2pt1":{
		"base":"cut_pt2pt1",
		"label":r"$p_{T2}/p_{T1} > 0.02$"
	},
	"Delta_R":{
		"base":"cut_Delta_R",
		"label": r'$\Delta R < 3.95$'
	}
 }


def build_mjj_histos_dict(ROOT_file,sample_type):
	assert sample_type in ["zg","gjets"]

	histos = {}
	for cut_key,cut_info in cuts.items():
		cut_name = cut_info["base"]
		
		# Inclusive mjj
		hist_name_inclusive = f'mjj_{sample_type}' if cut_name == '' else f"mjj_{cut_name}_{sample_type}"
		key = 'mjj' if cut_key == "nocut" else f"mjj_{cut_key}"
		histos[key] = ROOT_file.Get(hist_name_inclusive)

		# Flavour sensitive mjj
		for i in range(1,6):
			if cut_key == "nocut":
				flavour_name = f'mjj_partonflavour{i}_{sample_type}'
				key = f"mjj_partonflavour{i}"
			else:
				flavour_name = f'mjj_PartonFlavour_{i}_{cut_name}_{sample_type}'
				key = f"mjj_PartonFlavour{i}_{cut_key}"
			histos[key] = ROOT_file.Get(flavour_name)
	return histos

def plot_mjj_histograms(histos,cut_labels,output_dir="plots"):
	"""
	Plots mjj histograms and mjj by parton flavour for each cut level.

	Parameters:
	- histos	Dictionnary with each histogram name as key and ROOT histogram object as value.
	- cut_labels	Dictionnary mapping cut-keys to human-readable labels.
	- output_dir	Base directory to save plots

	"""

	os.makedirs(output_dir,exist_ok=True)
	
	for cut_key,label in cut_labels.items():
		cut_dir = os.path.join(output_dir,cut_key)
		os.makedirs(cut_dir,exist_ok=True)
		
		# Inclusive mjj histogram
		histo_name = 'mjj' if cut_key == "nocut" else f"mjj_{cut_key}"
		if histo_name in histos:
			hist = histos[histo_name]
			fig,ax = plt.subplots()
			x_vals = [hist.GetBinCenter(i) for i in range(1,hist.GetNbinsX()+1)]
			y_vals = [hist.GetBinContent(i) for i in range(1,hist.GetNbinsX()+1)]
			ax.plot(x_vals,y_vals,label=f'Inclusive',color='black')
			ax.set_title(f'Mjj - {label}')
			ax.set_xlabel(f'Mjj [GeV]')
			ax.set_ylabel(f'Entries')
			ax.legend()
			fig.tight_layout()
			fig_name = f'mjj_{cut_key}.png'
			fig.savefig(os.path.join(cut_dir,fig_name))
			plt.close(fig)
		else:
			print(f'{histo_name} not found !')
		# Flavour-tagged mjj histogram
		for k in range(1,6):
			hist_name = f'mjj_partonflavour{k}' if cut_key == 'nocut' else f'mjj_PartonFlavour{k}_{cut_key}'
			if hist_name not in histos:
				print(f'{hist_name} not found !')
			else:
				hist = histos[hist_name]
				fig,ax = plt.subplots()
				x_vals = [hist.GetBinCenter(i) for i in range(1,hist.GetNbinsX()+1)]
				y_vals = [hist.GetBinContent(i) for i in range(1,hist.GetNbinsX()+1)]
				hist_label = f'Flavour {k}'
				ax.plot(x_vals,y_vals,label=hist_label)
				hist_title = f'Mjj PartonFlavour {k} - {label}'
				ax.set_title(hist_title)
				ax.set_xlabel(f'Mjj [GeV]')
				ax.set_ylabel(f'Entries')
				ax.legend()
				fig.tight_layout()
				fig_name = f'mjj_PartonFlavour{k}' if cut_key == "nocut" else f"mjj_PartonFlavour{k}_{cut_key}" + ".png"
				fig.savefig(os.path.join(cut_dir,fig_name))
				plt.close(fig)

cut_labels = {key:vals["label"] for key,vals in cuts.items()}

histos_bg = build_mjj_histos_dict(file_bg, "gjets")
histos_sg = build_mjj_histos_dict(file_sg, "zg")

plot_mjj_histograms(histos_bg, cut_labels, output_dir="plots/bg")
plot_mjj_histograms(histos_sg, cut_labels, output_dir="plots/sg")
