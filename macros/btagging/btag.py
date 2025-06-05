import os
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import uproot

# Disable interactive mode for clean plotting
plt.ioff()

# Input file
source_file_path = "../Source.root"
background_file_path = '../Background.root'

# Output directory
base_output_dir = "histogram_outputs"

# Working points (efficiencies)
efficiencies = [0.7, 0.8, 0.9]


# Open ROOT file and return histograms
def ROOT_to_histo(file, isSignal=True):
	all_keys = file.keys()
	inclusive_histograms = []
	flavour_histograms = {}

	for key in all_keys:
		name = key.split(';')[0]
		if not name.startswith('Jet_btagPNetB'):
			continue
		# No PartonFlavour -> Inclusive
		if "PartonFlavour" not in name:
			inclusive_histograms.append(name)
			continue
		# PartonFlavour -> Put in dictionnary by flavour
		if isSignal:
			match = re.match(r'Jet_btagPNetB_PartonFlavour(\d+)(_\d+)?(_.*)?', name)
			if match:
				flavour = int(match.group(1))
				flavour_histograms.setdefault(flavour,[]).append(name)
	return (inclusive_histograms, flavour_histograms) if isSignal else inclusive_histograms	

# Function to calculate the cut value
def find_cut_value(hist_values, bin_edges, efficiency):
    """
    Given histogram counts and bin edges, return the threshold (cut value)
    such that 'efficiency' fraction of the total area is kept,
    starting from the highest values (right side of the histogram).
    """
    total = np.sum(hist_values)
    if total == 0:
        return bin_edges[-1]  # Avoid division by zero

    cumulative = np.cumsum(hist_values[::-1])  # Cumulative sum from right
    bin_edges_reversed = bin_edges[::-1]

    for i, value in enumerate(cumulative):
        if value >= efficiency * total:
            return bin_edges_reversed[i+1] if i + 1 < len(bin_edges_reversed) else bin_edges[0]

    return bin_edges[0]



# Function to plot and save histograms
def draw_and_save(hist_values, bin_edges, title, filepath, cut_lines=None):
    plt.figure(figsize=(8, 6))

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    width = np.diff(bin_edges)
    plt.bar(bin_centers, hist_values, width=width, align='center', color='black', edgecolor='black', linewidth=0.5)

    plt.title(title)
    plt.xlabel("Jet_btagPNetB")
    plt.ylabel("Entries")

    if cut_lines:
        for i, (eff, xcut) in enumerate(sorted(cut_lines.items(), reverse=True)):
            plt.axvline(x=xcut, color=f'C{i+1}', linestyle='--', linewidth=2, label=f'{int(eff*100)}% cut @ {xcut:.3f}')

    if cut_lines:
        plt.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig(filepath)
    plt.close()

def process_histograms(file,inclusive_histograms,output_dir,efficiencies,isSignal=True, flavour_histograms=None):
	"""
	Process inclusive and optionally flavour-separated histograms.
	"""
	# Inclusive
	inclusive_dir = os.path.join(output_dir,"Inclusive")
	os.makedirs(inclusive_dir,exist_ok=True)
	
	for hist_name in inclusive_histograms:
		if hist_name in file:
			hist = file[hist_name]
			hist_values, bin_edges = hist.to_numpy()
			cut_lines = {eff: find_cut_value(hist_values,bin_edges, eff) for eff in efficiencies}
			output_filepath = os.path.join(inclusive_dir,f"{hist_name}.png")
			draw_and_save(hist_values,bin_edges,hist_name, output_filepath, cut_lines)
		else:
			print(f'Warning: Inclusive histogram {hist_name} not found in file.')
		
	
	if isSignal and flavour_histograms:
		for flavour, hist_list in flavour_histograms.items():
			flavour_dir = os.path.join(output_dir,f"PartonFlavour{flavour}")
			os.makedirs(flavour_dir,exist_ok=True)
			
			for hist_name in hist_list:
				if hist_name in file:	
					hist = file[hist_name]
					hist_values, bin_edges = hist.to_numpy()
					cut_lines = {eff: find_cut_value(hist_values,bin_edges, eff) for eff in efficiencies}
					output_filepath = os.path.join(flavour_dir,f"{hist_name}.png")
					draw_and_save(hist_values,bin_edges,hist_name,output_filepath,cut_lines)
				else:
					print(f'Warning: Flavour histogram "{hist_name}" not found in file.')

for root_file_path in [source_file_path,background_file_path]:
	isSignal = (root_file_path == source_file_path)

	file = uproot.open(root_file_path)

	if isSignal:
		inclusive_histograms, flavour_histograms = ROOT_to_histo(file)
		output_dir = os.path.join(base_output_dir, "Signal")
		os.makedirs(output_dir, exist_ok=True)
	else:
		inclusive_histograms = ROOT_to_histo(file,isSignal=False)
		output_dir = os.path.join(base_output_dir, "Background")
		os.makedirs(output_dir, exist_ok=True)

	os.makedirs(output_dir,exist_ok=True)
	process_histograms(
		file,
		inclusive_histograms,
		output_dir,
		efficiencies,
		isSignal,
		flavour_histograms if isSignal else None
		)
