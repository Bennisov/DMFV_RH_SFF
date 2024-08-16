import uproot
import matplotlib.pyplot as plt
import numpy
import math
import numpy as np

cs_sig = 0.3665 #pb
lumi = 29.1 #1/fb
N = 28505
cs_sig = cs_sig * 1000 #fb
w_sig = cs_sig * lumi / N


bcks = ["others", "singletop", "ttbar", "ttZ", "wjets", "zjets"]
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
            background_error = background * b_err
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

file_sig = uproot.open("../DMRH/output_2.0_850.root")
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


def plotting(histos, histo_sig, xlab="", tit="", bck=bcks, dir=1):
    # Initialize values and edges
    values = [None] * 6
    edges = [None] * 6
    
    # Process each histogram
    for i in range(6):
        values[i], edges[i] = histos[i].to_numpy()
    
    # Calculate the sums and sort the indices based on these sums
    sums = np.array([np.sum(val) for val in values])  # Corrected sum calculation
    sorted_indices = np.argsort(sums)
    
    # Sort the values and tags manually
    sorted_values = [values[i] for i in sorted_indices]
    sorted_bck = [bck[i] for i in sorted_indices]
    
    # Signal values and errors
    values_sig, edges_sig = histo_sig.to_numpy()
    errors = np.sqrt(values_sig)
    errors *= w_sig
    values_sig *= w_sig
    
    # Initialize the figure and axes
    fig, axs = plt.subplots(2, 1, figsize=(8, 10), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    
    # Prepare cumulative histograms
    cumulative_values = [np.copy(sorted_values[i]) for i in range(6)]
    for i in range(1, 6):
        cumulative_values[i] += cumulative_values[i-1]
    
    # Plot histograms in reverse order for proper stacking
    for i in range(5, -1, -1):
        axs[0].hist(edges[0][:-1], bins=edges[0], weights=cumulative_values[i], histtype='stepfilled', 
                    label=sorted_bck[i], alpha=1.0)
    
    # Plot the signal histogram
    x = (edges[0][:-1] + edges[0][1:]) / 2.0
    axs[0].hist(edges[0][:-1], bins=edges[0], weights=values_sig, color='k', histtype='step', label='signal')
    
    # Log scale, labels, and legend
    axs[0].set_yscale('log')
    axs[0].legend(ncol=2)
    axs[0].set_xlabel(xlab)
    axs[0].set_ylabel("NoE")
    
    # Plot significance
    s_values = s_calc(histos=sorted_values, sig=values_sig, dire=dir)
    axs[1].plot(x, s_values, color='red')
    axs[1].set_xlabel(xlab)
    axs[1].set_ylabel('significance')
    
    plt.tight_layout()
    plt.savefig(tit)
    
plotting(histos=h_mt2, histo_sig=h_mt2_sig, xlab="Stranverse mass [GeV]", tit="mt2.png")
plotting(histos=h_met, histo_sig=h_met_sig, xlab="Missing Transverse Momentum [MeV]", tit="met.png")
plotting(histos=h_mt, histo_sig=h_mt_sig, xlab="Tranverse mass [GeV]", tit="mt.png")
plotting(histos=h_n_bjet, histo_sig=h_n_bjet_sig, xlab="Number of B-Jets", tit="nbjet.png")
plotting(histos=h_pt_l, histo_sig=h_pt_l_sig, xlab="Transverse momentum of lepton [GeV]", tit="ptlep.png")
plotting(histos=h_pt_bjet, histo_sig=h_pt_bjet_sig, xlab="Transverse momentum of B-Jet [GeV]", tit="ptbjet.png")
plotting(histos=h_pt_jet, histo_sig=h_pt_jet_sig, xlab="Tranverse momentum of jet [GeV]", tit="ptjet.png")
plotting(histos=h_m_bl, histo_sig=h_m_bl_sig, xlab="Invariant mass of B-Jet and lepton[GeV]", tit="mbl.png", dir=0)
plotting(histos=h_dphi_min, histo_sig=h_dphi_min_sig, xlab="Delta phi min of jet and MET", tit="dphimin.png")
#plotting(histos=h_n_l, histo_sig=h_n_l_sig, xlab="Number of lepton", tit="nlep.png")
plotting(histos=h_n_jets, histo_sig=h_n_jets_sig, xlab="Number of jets", tit="njet.png")
