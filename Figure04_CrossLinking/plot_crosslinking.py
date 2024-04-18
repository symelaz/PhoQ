import pandas as pd
import re
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import os
from mpl_axes_aligner import align

def plot_crosslinking(filename):
    with open(filename,"r") as f: content = f.readlines()
    splitted = [re.split(" |\n",line) for line in content]
    filtered = [list(filter(None, line)) for line in splitted]
    df = pd.DataFrame(filtered)
    df.drop_duplicates(inplace=True)

    return df.sort_values(by=0)

def read_data(filename):
    means = []
    stds = []
    df = pd.read_csv(filename,header=None)
    resids = np.unique(df[0])
    for resid in resids:
        # Fit a normal distribution to the data:
        mu, std = norm.fit(list(df[2][df[0]==resid]))
        means.append(mu)
        stds.append(std)

    return pd.DataFrame({"Resid": resids, "Mean": means, "STD": stds})

############################### Inputs ################################

folder = "/home/symela/Documents/Symela/PhD/Project/PhoQ/PhoQ/Analysis/CrossLinking"

exp_file = os.path.join(folder, "ExperimentalData-Cross-Linking/BayesianModeling/modeling/data/expcrosslink_full.dat")
cross_full = plot_crosslinking(exp_file)

af_full = read_data(os.path.join(folder, "CrossLinking_phoqaf_replicas.dat"))
c_full = read_data(os.path.join(folder, "CrossLinking_phoqc_replicas.dat"))

# Read AlphaFold predicted structure
os.system("vmd -dispdev text -e cross_linking_AF_predicted.tcl -args /home/symela/Documents/Symela/PhD/Project/PhoQ/PhoQ/phoq_af_model.pdb AF_predicted.dat > logs")
af_predicted = pd.read_csv("AF_predicted.dat", header=None, names=["Resid", "Mean"])
af_predicted["STD"] = 0


########################### TM2 + HAMP plot ###########################

df_af = af_full[ (af_full.Resid.astype(int) >= 185) & (af_full.Resid.astype(int) < 227) ]
df_c = c_full[ (c_full.Resid.astype(int) >= 185) & (c_full.Resid.astype(int) < 227) ]
df_af_predicted = af_predicted[ (af_predicted.Resid.astype(int) >= 185) & (af_predicted.Resid.astype(int) < 227) ]
exp = cross_full[cross_full[0].astype(int) > 170]

inv=False

fig, axs = plt.subplots(2,figsize=(1.9, 3))
ax1 = axs[1]
ax2 = axs[0]

#ax2 = ax1.twinx()

ax1.errorbar(df_af.Resid.astype(int), df_af.Mean.astype(float), yerr=df_af.STD.astype(float), fmt='-+', color="red",
             ecolor="red", elinewidth=0.5, capsize=0, linewidth=1, markersize=2, label="AlphaFold");
ax1.errorbar(df_c.Resid.astype(int), df_c.Mean.astype(float), yerr=df_c.STD.astype(float), fmt='-+', color="darkorange",
            ecolor="darkorange", elinewidth=0.5, capsize=0, linewidth=1, markersize=2, label="Hybrid");
#ax1.plot(df_af_predicted.Resid.astype(int), df_af_predicted.Mean.astype(float), '-*', color="gray", linewidth=1, markersize=2, label="Predicted Structure")
ax2.errorbar(exp[0].astype(int), 1 - exp[4].astype(float), yerr=exp[5].astype(float), marker='.', color='black',
             ecolor='lightgray', elinewidth=0.5, capsize=0, linewidth=1, markersize=2, label="Experimental Data");

ax2.set_xlabel('Residue ID', fontsize=7)
ax1.set_ylabel(r'Distance ($\AA$)', color='black', fontsize=7, rotation=90)
ax2.set_ylabel('1 - Cross-Linking Fraction', color='black', fontsize=7, rotation=90)
ax1.tick_params(axis='both', which='major', labelsize=6, width=0.5, rotation=90, direction='inout')
ax1.tick_params(axis='both', which='minor', labelsize=4, width=0.5, direction='inout')
ax2.tick_params(axis='both', which='major', labelsize=6, width=0.5, rotation=90, direction='inout')
ax2.tick_params(axis='both', which='minor', labelsize=4, width=0.5, direction='inout')    
ax1.set_ylim(3,32)
ax2.set_ylim(0, 1.2)
if inv: 
    ax1.set_xlim(df_af.Resid.astype(int).max() + 1, df_af.Resid.astype(int).min() - 1)
    ax2.set_xlim(exp[0].astype(int).max() + 1, exp[0].astype(int).min() - 1)
ax1.legend()
ax2.legend()
plt.savefig("TM2+HAMP.svg", transparent=True)


####################### SD interface + TM1 plot #######################

df_af = af_full[af_full.Resid.astype(int) < 65]
df_c = c_full[c_full.Resid.astype(int) < 65]
df_af_predicted = af_predicted[af_predicted.Resid.astype(int) < 65]
exp = cross_full[cross_full[0].astype(int) < 65]

inv=False

fig, axs = plt.subplots(2,figsize=(1.9, 3))
ax1 = axs[1]
ax2 = axs[0]

ax1.errorbar(df_af.Resid.astype(int), df_af.Mean.astype(float), yerr=df_af.STD.astype(float), fmt='-+', color="red",
             ecolor="red", elinewidth=0.5, capsize=0, linewidth=1, markersize=2, label="AlphaFold");
ax1.errorbar(df_c.Resid.astype(int), df_c.Mean.astype(float), yerr=df_c.STD.astype(float), fmt='-+', color="darkorange",
            ecolor="darkorange", elinewidth=0.5, capsize=0, linewidth=1, markersize=2, label="Hybrid");
#ax1.plot(df_af_predicted.Resid.astype(int), df_af_predicted.Mean.astype(float), '-*', color="gray", linewidth=1, markersize=2, label="Predicted Structure")
ax2.errorbar(exp[0].astype(int), 1 - exp[4].astype(float), yerr=exp[5].astype(float), marker='.', color='black',
             ecolor='lightgray', elinewidth=0.5, capsize=0, linewidth=1, markersize=2, label="Experimental Data");

ax2.set_xlabel('Residue ID', fontsize=7)
ax1.set_ylabel(r'Distance ($\AA$)', color='black', fontsize=7, rotation=-90)
ax2.set_ylabel('1 - Cross-Linking Fraction', color='black', fontsize=7, rotation=-90)
ax1.tick_params(axis='both', which='major', labelsize=6, width=0.5, rotation=-90, direction='inout')
ax1.tick_params(axis='both', which='minor', labelsize=4, width=0.5, direction='inout')
ax2.tick_params(axis='both', which='major', labelsize=6, width=0.5, rotation=-90, direction='inout')
ax2.tick_params(axis='both', which='minor', labelsize=4, width=0.5, direction='inout')    
ax1.set_ylim(3,27)
ax2.set_ylim(0, 1.2)
if inv: 
    ax1.set_xlim(df_af.Resid.astype(int).max() + 1, df_af.Resid.astype(int).min() - 1)
    ax2.set_xlim(exp[0].astype(int).max() + 1, exp[0].astype(int).min() - 1)
ax1.legend()
ax2.legend()
plt.savefig("TM1+SensorInterface.svg", transparent=True)

