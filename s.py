import uproot
import matplotlib.pyplot as plt
import numpy
import math
import numpy as np



bcks = ["Others", "Singletop", "Ttbar", "TtZ", "Wjets", "Zjets"]
files = [None] * 6
for i in range(6):
    files[i] = bcks[i] + "_output.root"

def s_calc(histos, sig, dire, b_err=0.3):
    size = histos[0].size
    bck = [0.] * size
    for i in range(6):
        for j in range(size):
            bck[j] = bck[j] + histos[i][j]
    s = [0.] * size
    if dire:
        for i in range(size):
            signal = sum(sig[i:])
            background = sum(bck[i:])
            background_error = background * b_err
            s[i] = getZnGlenCowen(s=signal, b=background, b_err_abs=background_error)
    else:
        for i in range(size):
            signal = sum(sig[:(i+1)])
            background = sum(bck[:(i+1)])
            background_error = math.sqrt(background) * b_err
            s[i] = getZnGlenCowen(s=signal, b=background, b_err_abs=background_error)
    return s
        
def getZnGlenCowen(s,b,b_err_abs):
    tot = s+b
    b2 = b*b
    b_err2 = b_err_abs*b_err_abs
    b_plus_err2 = b+b_err2
    retval=0
    try:
        retval= math.sqrt(2*((tot) * math.log(tot*b_plus_err2/(b2+tot*b_err2))-b2/b_err2 * math.log(1+b_err2*s/(b*b_plus_err2))))
    except ValueError:
        pass
    return retval

file_uproot = [None] * 6
for i in range(6):
    file_uproot[i] = uproot.open(files[i])

h_mt2 = [None] * 6
h_mt = [None] * 6
h_met = [None] * 6
h_pt_l = [None] * 6
h_n_bjet = [None] * 6
h_pt_bjet = [None] * 6
h_pt_jet = [None] * 6
h_m_bl = [None] * 6
h_dphi_min = [None] * 6
h_n_l = [None] * 6
h_n_jets = [None] * 6
cutflows = [None] * 6

for i in range(6):
    h_mt2[i] = file_uproot[i]["h_mt2"]
    h_met[i] = file_uproot[i]["h_met"]
    h_mt[i] = file_uproot[i]["h_mt"]
    h_n_bjet[i] = file_uproot[i]["h_bjet_n"]
    h_pt_l[i] =  file_uproot[i]["h_pt_l"]
    h_pt_bjet[i] =  file_uproot[i]["h_pt_bjet"]
    h_pt_jet[i] =  file_uproot[i]["h_pt_jet"]
    h_m_bl[i] =  file_uproot[i]["h_m_bl"]
    h_dphi_min[i] =  file_uproot[i]["h_dphi_min"]
    h_n_l[i] =  file_uproot[i]["h_n_l"]
    h_n_jets[i] =  file_uproot[i]["h_n_jets"]
    cutflows[i] = file_uproot[i]["cutflow"]

file_sig = uproot.open("../bck/signal_output.root")
h_mt2_sig = file_sig["h_mt2"]
h_met_sig = file_sig["h_met"]
h_mt_sig = file_sig["h_mt"]
h_n_bjet_sig = file_sig["h_bjet_n"]
h_pt_l_sig =  file_sig["h_pt_l"]
h_pt_bjet_sig =  file_sig["h_pt_bjet"]
h_pt_jet_sig =  file_sig["h_pt_jet"]
h_m_bl_sig =  file_sig["h_m_bl"]
h_dphi_min_sig =  file_sig["h_dphi_min"]
h_n_l_sig =  file_sig["h_n_l"]
h_n_jets_sig =  file_sig["h_n_jets"]
cutflow_sig = file_sig["cutflow"]

file_data = uproot.open("../bck/Data_output.root")
h_mt2_data = file_data["h_mt2"]
h_met_data = file_data["h_met"]
h_mt_data = file_data["h_mt"]
h_n_bjet_data = file_data["h_bjet_n"]
h_pt_l_data =  file_data["h_pt_l"]
h_pt_bjet_data =  file_data["h_pt_bjet"]
h_pt_jet_data =  file_data["h_pt_jet"]
h_m_bl_data =  file_data["h_m_bl"]
h_dphi_min_data =  file_data["h_dphi_min"]
h_n_l_data =  file_data["h_n_l"]
h_n_jets_data =  file_data["h_n_jets"]
cutflow_data = file_data["cutflow"]

def plotting(histos, histo_sig, histo_data, xlab="", tit="", bck=bcks, dir=1):
    # Initialize values and edges
    values = [None] * 6
    edges = [None] * 6
    
    # Process each histogram
    for i in range(6):
        values[i], edges[i] = histos[i].to_numpy()
    
    # Calculate the sums and sort the indices based on these sums
    sums = np.array([np.sum(val) for val in values])
    sorted_indices = np.argsort(sums)
    
    # Sort the values and tags manually
    sorted_values = [values[i] for i in sorted_indices]
    sorted_bck = [bck[i] for i in sorted_indices]

    # Signal values and errors
    values_sig, edges_sig = histo_sig.to_numpy()
    errors = np.sqrt(values_sig)
    
    values_data, edges_data = histo_data.to_numpy()
    
    # Initialize the figure and axes
    fig, axs = plt.subplots(2, 1, figsize=(16, 20), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    
    # Prepare cumulative histograms
    cumulative_values = [np.copy(sorted_values[i]) for i in range(6)]
    for i in range(1, 6):
        cumulative_values[i] += cumulative_values[i-1]
    
    # Plot histograms in reverse order for proper stacking
    for i in range(5, -1, -1):
        axs[0].hist(edges[0][:-1], bins=edges[0], weights=cumulative_values[i], histtype='stepfilled', 
                    label=sorted_bck[i], alpha=1.0)
    x = (edges[0][:-1] + edges[0][1:]) / 2.0
    axs[0].hist(edges[0][:-1], bins=edges[0], weights=values_sig, histtype='step', label='signal', color='black')
    axs[0].scatter(x, values_data, label="data", marker='o', alpha=1., s=300, c='black')
    # Plot the signal histogram
    
    #axs[0].hist(edges[0][:-1], bins=edges[0], weights=values_sig, color='k', histtype='step', label='signal')
    
    # Log scale, labels, and legend
    axs[0].set_yscale('log')
    axs[0].legend(ncol=2)
    axs[0].set_xlabel(xlab)
    axs[0].set_ylabel("NoE")
    
    # Add grid and customize y-axis ticks
    # axs[0].grid(True, which='both', linestyle='--', linewidth=0.5)
    # axs[0].xaxis.set_major_locator(MultipleLocator((edges[0][-1] - edges[0][0])/10.))  

    # Plot significance
    s_values = s_calc(histos=sorted_values, sig=values_sig, dire=dir)
    axs[1].plot(x, s_values, color='red', label='long-ass formula')

    axs[1].legend()
    axs[1].set_xlabel(xlab)
    axs[1].set_ylabel('significance')
    
    # Add grid and customize y-axis ticks to the lower plot as well
    # axs[1].grid(True, which='both', linestyle='--', linewidth=0.5) 
    # axs[1].xaxis.set_major_locator(MultipleLocator((edges[0][-1] - edges[0][0])/10.))  
    
    # Adjust layout and save the figure
    plt.tight_layout()
    plt.savefig(tit)

    
plotting(histos=h_mt2, histo_sig=h_mt2_sig, histo_data=h_mt2_data, xlab="Stranverse mass [GeV]", tit="mt2.png")
plotting(histos=h_met, histo_sig=h_met_sig, histo_data=h_met_data, xlab="Missing Transverse Momentum [MeV]", tit="met.png")
plotting(histos=h_mt, histo_sig=h_mt_sig, histo_data=h_mt_data, xlab="Tranverse mass [GeV]", tit="mt.png")
plotting(histos=h_n_bjet, histo_sig=h_n_bjet_sig, histo_data=h_n_bjet_data, xlab="Number of B-Jets", tit="nbjet.png")
plotting(histos=h_pt_l, histo_sig=h_pt_l_sig, histo_data=h_pt_l_data, xlab="Transverse momentum of lepton [GeV]", tit="ptlep.png")
plotting(histos=h_pt_bjet, histo_sig=h_pt_bjet_sig, histo_data=h_pt_bjet_data, xlab="Transverse momentum of B-Jet [GeV]", tit="ptbjet.png")
plotting(histos=h_pt_jet, histo_sig=h_pt_jet_sig, histo_data=h_pt_jet_data, xlab="Tranverse momentum of jet [GeV]", tit="ptjet.png")
plotting(histos=h_m_bl, histo_sig=h_m_bl_sig, histo_data=h_m_bl_data, xlab="Invariant mass of B-Jet and lepton[GeV]", tit="mbl.png", dir=0)
plotting(histos=h_dphi_min, histo_sig=h_dphi_min_sig, histo_data=h_dphi_min_data, xlab="Delta phi min of jet and MET", tit="dphimin.png")
plotting(histos=h_n_l, histo_sig=h_n_l_sig,histo_data=h_n_l_data ,xlab="Number of lepton", tit="nlep.png")
plotting(histos=h_n_jets, histo_sig=h_n_jets_sig, histo_data=h_n_jets_data, xlab="Number of jets", tit="njet.png")


cut_val_sig, e = cutflow_sig.to_numpy()
cut_val_dat, e = cutflow_data.to_numpy()

cut_val_bck = [None] * 6 

for i in range(6):
    cut_val_bck[i], _ = cutflows[i].to_numpy()
signals = cut_val_sig[-1]
backgrounds = numpy.sum(cut_val_bck[:][-1])
significance = getZnGlenCowen(signals, backgrounds, 0.3 * backgrounds)
print(signals)
print(backgrounds)
print("Significance = " + str(significance))
print("s/sqrt(b) = "+str(signals / math.sqrt(backgrounds)))
sig_alt = signals / math.sqrt(backgrounds + (0.3) * (0.3) * backgrounds * backgrounds)
print("sig_alt = "+str(sig_alt))

cut_names = numpy.array(range(12))

sums = np.array([np.sum(val) for val in cut_val_bck])
sorted_indices = np.argsort(sums)
sorted_values = [cut_val_bck[i] for i in sorted_indices]
sorted_bck = [bcks[i] for i in sorted_indices]

fig, axs = plt.subplots(2, 1, figsize=(16, 20), sharex=True, gridspec_kw={'height_ratios': [3, 1]})

cumulative_values = [np.copy(cut_val_bck[i]) for i in range(6)]
for i in range(1, 6):
    cumulative_values[i] += cumulative_values[i-1]
for i in range(5, -1, -1):
    axs[0].bar(cut_names, cumulative_values[i], label=sorted_bck[i])
axs[0].bar(cut_names, cut_val_sig, label="Signal", color='black', alpha=0.5)
axs[0].scatter(cut_names, cut_val_dat, label="Data", color='black')
axs[0].legend(ncol=2)
axs[0].set_yscale('log')
axs[0].set_ylabel("NoE")
sig = [None] * 12
sig_alt = [None] * 12
sig_alt2 = [None] * 12
for i in range(12):
    sig[i] = getZnGlenCowen(cut_val_sig[i], cumulative_values[5][i], (cumulative_values[5][i])*(0.3))
    sig_alt[i] = cut_val_sig[i] / math.sqrt(cumulative_values[5][i] + (0.3) * (0.3) * cumulative_values[5][i] * cumulative_values[5][i])
    sig_alt2[i] = cut_val_sig[i] / math.sqrt(cumulative_values[5][i])
axs[1].plot(cut_names, sig, color='red', label='Long formula')
#axs[1].plot(cut_names, sig_alt, color='blue', label='s/sqrt(b + sigma^2)')
#axs[1].plot(cut_names, sig_alt2, color='green', label='s/sqrt(b)')
axs[1].legend()
axs[1].set_ylabel("Significance")
plt.tight_layout()
plt.savefig("cutflow.png")