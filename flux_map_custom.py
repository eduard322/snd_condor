#!/usr/bin/env python2
from __future__ import division
import argparse
import numpy as np
import ROOT as r
# Fix https://root-forum.cern.ch/t/pyroot-hijacks-help/15207 :
r.PyConfig.IgnoreCommandLineOptions = True
import shipunit as u
import rootUtils as ut
import logger as log
from array import array


def read_files(chain, eos, filepath, filename, file_num):
    for i in range(1, file_num + 1):
        chain.Add(eos + filepath + "/" + str(i) + "/" + filename)

def main():
    parser = argparse.ArgumentParser()
    #parser.add_argument(
    #    'inputfile',
    #    help='''Simulation results to use as input. '''
    #    '''Supports retrieving files from EOS via the XRootD protocol.''')
    parser.add_argument(
        '-o',
        '--outputfile',
        default='flux_map.root',
        help='''File to write the flux maps to. '''
        '''Will be recreated if it already exists.''')
    # parser.add_argument(
    # '-P',
    # '--Pcut',
    # default= 0.0,
    # help='''set momentum cut''')
    # parser.add_argument(
    # '-E',
    # '--Eloss',
    # default= 0.0,
    # help= '''set Eloss cut''')
    parser.add_argument(
        '--nStart',
        dest="nStart",
        default=0)
    
    parser.add_argument(
        '--nEvents',
        dest="nEvents",
        default=1)
    
        
    args = parser.parse_args()
    f = r.TFile.Open(args.outputfile, 'recreate')
    h = {}
    f.cd()

    eos = "root://eosuser.cern.ch/"
    filepath = "/eos/user/e/ekhaliko/Documents/SND_Data/test_E100-2500_n30k_center_pi-"
    filename = "sndLHC.PG_-211-TGeant4.root"
    file_num = 1
    ch = r.TChain('cbmsim')
    read_files(ch, eos, filepath, filename, file_num)
    # Define histograms
    ut.bookHist(h, 'Mufilter_hits', ';Nhits;', 100, 0, 20000)
    ut.bookHist(h, 'Mufilter_energy', ';Beam energy [GeV];', 100, 0, 2500)
    ut.bookHist(h, 'Mufilter', ';Nhits; Beam energy [GeV]', 100, 0, 20000, 100, 0, 2500)
    
    #ch.Add(args.inputfile)
    n = ch.GetEntries()

    i = 0
    P_dist = []
    mc_hits_dist = []
    for N in range(Nev_st, Nev_en + 1):
#  for event in eventTree:
   #  N+=1
    eventTree.GetEvent(N)
        if i % 1000 == 0:
            print("Event " + str(i))
            pass
        i += 1
        mc_hits = 0
        P = np.sqrt(event.MCTrack[0].GetPx()**2 + event.MCTrack[0].GetPy()**2 + event.MCTrack[0].GetPz()**2)
        #print(P)
        for hit in event.MuFilterPoint:
            if hit:
                if not hit.GetEnergyLoss() > 0:
                    continue
                pid = hit.PdgCode()
                
                # exclude neutrals
                if pid in [22, 2112, 111, 221]: continue


                #P = np.sqrt(hit.GetPx()**2 + hit.GetPy()**2 + hit.GetPz()**2)
                # set momentum cut
                if P < args.Pcut: continue

                Eloss = hit.GetEnergyLoss()
                # set Eloss cut
                if Eloss < args.Eloss: continue
                     
                mc_hits += 1

        h['Mufilter_hits'].Fill(mc_hits)
        h['Mufilter_energy'].Fill(P)
        h['Mufilter'].Fill(mc_hits, P)
        mc_hits_dist.append(mc_hits)
        P_dist.append(P)
    #print(Eloss_total)
    for key in h:
        classname = h[key].Class().GetName()
        if 'TH' in classname or 'TP' in classname:
            h[key].Write()
    mc_hits_parray = r.TVectorD(len(mc_hits_dist), array('d', mc_hits_dist))
    P_parray = r.TVectorD(len(P_dist), array('d', P_dist))
    mc_hits_parray.Write('mc_hits')
    P_parray.Write('P')

    f.Close()
    #np.savetxt('test_out', np.array(Eloss_total), fmt='%f') 
    #f = open(args.norm + "/flux", "w")
    #f.write(str(B_ids_unw) + "\t" + str(B_ids))
    #f.close()


if __name__ == '__main__':
    r.gErrorIgnoreLevel = r.kWarning
    r.gROOT.SetBatch(True)
    main()
