import os

for run in ["100631", "100643"]:
    os.system(f"python new_mufi_exp.py -p /eos/experiment/sndlhc/convertedData/commissioning/testbeam_June2023_H8/ -g geofile_full.Ntuple-TGeant4_nom.root --nStart 0 --nEvents 3000000 -r {run}")
