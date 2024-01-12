import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import uproot
import awkward
import multiprocessing as mp
import matplotlib as mpl
import seaborn as sns
import os
import json
sns.set_style("whitegrid")

def largeSiPMchannel(i):
    if i==2 or i==5 or i==10 or i==13: return False
    else: return True
    
def smallSiPMchannel(i):
    if i==2 or i==5 or i==10 or i==13: return True
    else: return False
    
def read_data_files(init_config, time_included = False):
    Data = {}
    for run in init_config:
        fileName = f"{run}"
        basePath = sorted(Path(init_config[run]).glob(f'**/{fileName}'))
        print("{} files to read in {}".format(len(basePath), init_config[run]))
        dfs = []
        for base in basePath:
            try:
                dfs.append(pd.read_csv(base, header = None, sep = "\t"))
            except:
                pass
        Data[run] = pd.concat(dfs, ignore_index=True)
        df_cols = ["event_number", "hit", "det_id"] + [f"ch_{i}" for i in range(16)]
        if time_included:
            Data[run].columns = ["event_number", "hit", "det_id"] + [f"ch_{i}" for i in range(16)] + [f"time_{i}" for i in range(16)]
        else:
            Data[run].columns = ["event_number", "hit", "det_id"] + [f"ch_{i}" for i in range(16)]
#         Data[run].columns = df_cols
        Data[run]["l"] = (Data[run]["det_id"]%10000)//1000  
        Data[run]["bar"] = (Data[run]["det_id"]%1000)
        Data[run]["s"] = Data[run]["det_id"]//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
        Data[run]["bar_ID"] = Data[run]["s"]*100 + Data[run]["l"]*10 + Data[run]["bar"]
    return Data
def prepare_data(Data_origin, selected = False, wall_origin = None, small_large_good = False):
    data_info = {}
    Data = Data_origin.copy()
    for key in Data:
        data_info[key] = {}
        if selected:
            selected_events = Data[key].groupby("event_number").count().query("hit >= 5").reset_index(0).event_number
            Data[key] = Data[key].loc[Data[key]["event_number"].isin(selected_events)].copy()
        if not wall_origin is None:
            with open(f"scifi_events_{key}.pkl", "rb") as f:
                scifi_info = pd.DataFrame(pickle.load(f))
                if wall_origin != 0:
                    walls_events = scifi_info.query(f"wall == {wall_origin}")["event"]
                    Data[key] = Data[key].loc[Data[key]["event_number"].isin(walls_events)].copy()
                else:
                    walls_events = scifi_info["event"]
                    Data[key] = Data[key].loc[~Data[key]["event_number"].isin(walls_events)].copy()
        if small_large_good:
            Data[key] = Data[key][Data[key] > -10]
        data_info[key]["Data"] = Data[key]
        data_info[key]["Data_small"] = Data[key][["event_number", "hit", "det_id", "l", "bar"] + [f"ch_{i}" for i in list(filter(smallSiPMchannel, range(16)))]].copy()
        data_info[key]["Data_large"] = Data[key].loc[:,~Data[key].columns.isin([f"ch_{i}" for i in list(filter(smallSiPMchannel, range(16)))])].copy()
        data_info[key]["Data_large"]["signal_sum"] = data_info[key]["Data_large"][[f"ch_{i}" for i in list(filter(largeSiPMchannel, range(16)))]][(data_info[key]["Data_large"] > -10) & (data_info[key]["Data_large"] < 120)].sum(axis = 1).copy()
        data_info[key]["Data_small"]["signal_sum"] = data_info[key]["Data_small"][[f"ch_{i}" for i in list(filter(smallSiPMchannel, range(16)))]][(data_info[key]["Data_small"] > -10)].sum(axis = 1)
        data_info[key]["signal_per_event"] = data_info[key]["Data_large"].groupby(["event_number"]).sum()["signal_sum"]
        data_info[key]["fired_bars_per_event"] = data_info[key]["Data_large"].groupby(["event_number"]).count()["signal_sum"]
        data_info[key]["signal_per_plane_per_event"] = data_info[key]["Data_large"].groupby(["event_number", "l"]).sum()["signal_sum"]
        data_info[key]["signal_per_event_small"] = data_info[key]["Data_small"].groupby(["event_number"]).sum()["signal_sum"]
        data_info[key]["fired_bars_per_event_small"] = data_info[key]["Data_small"].groupby(["event_number"]).count()["signal_sum"]
        data_info[key]["signal_per_plane_per_event_small"] = data_info[key]["Data_small"].groupby(["event_number", "l"]).sum()["signal_sum"]
    return data_info
def vis_data(Data, wall = None):
    fig, ax = plt.subplots(2, 3, figsize=(10*len(Data.keys()), 16), dpi = 300)
    energies = [100, 180, 300]
    for i, run in enumerate(Data):
        #print(Data[run])
        h = ax[0][i].hist2d(Data[run]["fired_bars_per_event"], Data[run]["signal_per_event"]/Data[run]["fired_bars_per_event"], 
                     bins = (50, 200), norm=mpl.colors.LogNorm(),
                     cmap = plt.cm.jet)
        ax[0][i].set_xlim([0,50])
#         ax[i].set_ylim([-200,1500])
        ax[0][i].set_xlabel("Number of fired bars per event", fontsize = 16)
        ax[0][0].set_ylabel("(Signal/Number of fired bars) per event [QDC]", fontsize = 16)
        if wall is None:
            ax[0][i].set_title(f"{energies[i]} GeV. Large SiPMs: -10 < signal < 120", fontsize = 18)
        else:
            ax[0][i].set_title(f"{energies[i]} GeV. Wall {wall}. Large SiPMs: -10 < signal < 120", fontsize = 18)
            if wall == 0:
                ax[0][i].set_title(f"{energies[i]} GeV. Not in Scifi. Large SiPMs: -10 < signal < 120", fontsize = 18)          
        try:
            fig.colorbar(h[3], ax=ax[0][i])
        except:
            pass
#     cax,kw = mpl.colorbar.make_axes([axs for axs in ax.flat])
#     plt.colorbar(h, cax=cax, **kw)
    for i, run in enumerate(Data):
        h = ax[1][i].hist2d(Data[run]["fired_bars_per_event_small"], Data[run]["signal_per_event_small"]/Data[run]["fired_bars_per_event_small"], 
                     bins = (50, 200), norm=mpl.colors.LogNorm(),
                     cmap = plt.cm.jet)
        ax[1][i].set_xlim([0,50])
        ax[1][i].set_ylim([-30,150])
        ax[1][i].set_xlabel("Number of fired bars per event", fontsize = 16)
        ax[1][0].set_ylabel("(Signal/Number of fired bars) per event [QDC]", fontsize = 16)
        if wall is None:
            ax[1][i].set_title(f"{energies[i]} GeV. Small SiPMs", fontsize = 18)
        else:
            ax[1][i].set_title(f"{energies[i]} GeV. Wall {wall}. Small SiPMs", fontsize = 18)
            if wall == 0:
                ax[1][i].set_title(f"{energies[i]} GeV. Not in Scifi. Small SiPMs", fontsize = 18)      
        try:
            fig.colorbar(h[3], ax=ax[1][i])
        except:
            pass
    #fig.savefig("signal_fired_bars.pdf")
#     cax,kw = mpl.colorbar.make_axes([axs for axs in ax.flat])
#     plt.colorbar(h, cax=cax, **kw)
    fig.show()
    fig, ax = plt.subplots(1, 2, figsize = (16,8), dpi = 300)
    if wall is None:
        ax[0].set_title("Large SiPMs")
        ax[1].set_title("Small SiPMs")
    else:
        ax[0].set_title(f"Large SiPMs. Wall {wall}")
        ax[1].set_title(f"Small SiPMs. Wall {wall}")    
        if wall == 0:
            ax[0].set_title(f"Large SiPMs. Not in Scifi", fontsize = 18)
            ax[1].set_title(f"Small SiPMs. Not in Scifi", fontsize = 18)             
    for j, run in enumerate(Data):
        signal_l = []
        signal_l_err = []
        signal_l_small = []
        signal_l_small_err = []
        for i in range(5):
            #print(Data[run]["signal_per_plane_per_event"].reset_index("l").query(f"l == {i}").mean()["signal_sum"])
            signal_l.append(Data[run]["signal_per_plane_per_event"].reset_index("l").query(f"l == {i}").mean()["signal_sum"])
            signal_l_small.append(Data[run]["signal_per_plane_per_event_small"].reset_index("l").query(f"l == {i}").mean()["signal_sum"])
            signal_l_err.append(Data[run]["signal_per_plane_per_event"].reset_index("l").query(f"l == {i}").std()["signal_sum"])
            signal_l_small_err.append(Data[run]["signal_per_plane_per_event_small"].reset_index("l").query(f"l == {i}").std()["signal_sum"])
        print(signal_l, signal_l_small)
        for signal, signal_err, axs in zip([signal_l, signal_l_small], [signal_l_err, signal_l_small_err], ax):
            axs.errorbar([plane - 0.1 + j*0.1 for plane in range(1,6)], signal, yerr = signal_err, label = f"{energies[j]} GeV", fmt = "o", capsize=10)
            axs.set_xticks([1, 2, 3, 4, 5])
            axs.set_xticklabels([f"Plane {i}" for i in range(1,6)])
            axs.set_xlabel("Planes", fontsize = 16)
            axs.set_ylabel("Average signal per event [QDC]", fontsize = 16)
            axs.legend(loc = 'upper right')        
    
    
def signal_extrapolate(Data):
    output = {}
    def fun(x, i, bar, slope_info, small_sipm_signal_mean):
        print("check_1")
        if bar % 10 not in [4,5,6]:
            return x
        slope = slope_info.query(f"bar_id == {bar} & ch_other == {i}").apply(np.mean)["slope"]
        slope_err = slope_info.query(f"bar_id == {bar} & ch_other == {i}").apply(np.mean)["slope_err"]

#         print("check_1")
        if x > 140:
            return np.random.normal(slope*small_sipm_signal_mean, slope_err* small_sipm_signal_mean)
        else:
            return x
        
        
    def fun_2(x, slope_info, small_sipm_signal_mean):
#         print(x)
        bar_id = x["bar_ID"]
        if bar_id % 10 not in [4,5,6]:
            return x
        def get_slope(bar_id, slope_info, i):
            try:
                return slope_info.query(f"bar_id == {bar_id} & ch_other == {i}").apply(np.mean)["slope"]
            except:
                return 1
        def get_slope_err(x, slope_info, i):
            try:
                return slope_info.query(f"bar_id == {bar_id} & ch_other == {i}").apply(np.mean)["slope_err"]
            except:
                return 0

        
        x_0 = x.reset_index().copy()
        x_0.columns = ["index", "value"]
        z = x_0.apply(lambda y: np.random.normal(get_slope(bar_id, slope_info, int(y["index"].split("_")[1]))*small_sipm_signal_mean.iloc[x["index_y"]],\
                                               get_slope_err(bar_id, slope_info, int(y["index"].split("_")[1]))* np.abs(small_sipm_signal_mean.iloc[x["index_y"]]))  if (y["value"] > 140 and y["index"] in [f"ch_{i}" for i in range(16)]) else y["value"], axis = 1)
        return z
                                                
                                                 
                                                                                  
    
    for key in Data:
        Data[key] = Data[key][Data[key] != -999]
#         print(Data[key].columns)
        with open(f"large_small_cor_{key}.json", "r") as read_content: 
            slope_info = pd.DataFrame(json.load(read_content))
            slope_info.columns = ["bar_id", "slopes", "slopes_err", "ch_small", "ch_other"]
            Data_large_sipm_only = Data[key][[f"ch_{i}" for i in list(filter(largeSiPMchannel, range(16)))]].copy()
            small_sipm_signal_mean = Data[key].loc[:,Data[key].columns.isin([f"ch_{i}" for i in list(filter(smallSiPMchannel, range(16)))])].copy().mean(axis = 1)
            #output["key"] = Data[key][[f"ch_{i}" for i in range(16)]].apply(lambda col: fun(col, int(col.name.split("_")[1]), Data[key]["bar_ID"], slope_info, small_sipm_signal_mean))
            Data[key] = Data[key].reset_index().copy()
            Data[key].rename(columns = {"index": "index_y"}, inplace=True)
            output["key"] = Data[key].apply(lambda row: fun_2(row, slope_info, small_sipm_signal_mean), axis = 1)
            output["key"][["event_number", "hit", "det_id", "l", "bar", "s", "bar_ID"]] = Data[key][["event_number", "hit", "det_id", "l", "bar", "s", "bar_ID"]]
#             print(slope_info)
    return output        


Data = {}
Data["0"] = pd.read_csv("0", header = None, sep = "\t")
df_cols = ["event_number", "hit", "det_id"] + [f"ch_{i}" for i in range(16)]
if time_included:
    Data[run].columns = ["event_number", "hit", "det_id"] + [f"ch_{i}" for i in range(16)] + [f"time_{i}" for i in range(16)]