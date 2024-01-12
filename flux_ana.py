#!/usr/bin/env python2
from __future__ import division
import argparse
import numpy as np
import ROOT as r
# Fix https://root-forum.cern.ch/t/pyroot-hijacks-help/15207 :
r.PyConfig.IgnoreCommandLineOptions = True
#import shipunit as u
#import rootUtils as ut
#import logger as log


def main():
    parser = argparse.ArgumentParser(description='Script to create flux maps.')
    parser.add_argument(
        'inputfile',
        help='''Simulation results to use as input. '''
        '''Supports retrieving files from EOS via the XRootD protocol.''')
    parser.add_argument(
        '-o',
        '--outputfile',
        default='flux_map.root',
        help='''File to write the flux maps to. '''
        '''Will be recreated if it already exists.''')
    parser.add_argument(
        '-n',
        '--norm',
        default='flux',
        help='''File to write the flux maps to. '''
        '''Will be recreated if it already exists.''')
    
    
    args = parser.parse_args()
    f = r.TFile.Open(args.outputfile, 'recreate')
    f.cd()
    #maxpt = 10. * u.GeV
    #maxp = 360. * u.GeV
    h = {}

    # Define histograms
    #ut.bookHist(h, 'Mufilter', '#mu#pm;Eloss[GeV];', 100, 0, 100)
    ch = r.TChain('cbmsim')
    ch.Add(args.inputfile)
    n = ch.GetEntries()
    B_ids = 0
    B_ids_unw = 0
    SUM = 0
    Eloss_total = {20:[], 21:[], 22:[], 23:[], 24:[], 30:[], 31:[], 32:[]}
    #log.info(n)
    i = 0
    for event in ch:
        ids = {}
        if i % 10000 == 0:
            #log.info('{}/{}'.format(i, n))
            pass
        i += 1
        Eloss = 0
        weight = event.MCTrack[0].GetWeight()
        k = 1
        for hit in event.MuFilterPoint:
            if hit:
                if not hit.GetEnergyLoss() > 0:
                    continue
                pid = hit.PdgCode()
                if abs(pid) == 13:
                    P = np.sqrt(hit.GetPx()**2 + hit.GetPy()**2 + hit.GetPz()**2) 
                    trid = hit.GetTrackID()
                    if k:
                        x = hit.GetDetectorID() // 1000
                        k = 0
                    detector_ID = hit.GetDetectorID()
                    station = detector_ID // 1000
                    if station != x:
                        if x not in Eloss_total.keys(): continue
                        Eloss_total[x].append(Eloss) 
                        Eloss = 0
                    x = station
                    Eloss += hit.GetEnergyLoss()
                    if station <= 33: continue
                    print(i, hit.GetZ(), P, hit.GetEnergyLoss(), station)
        if Eloss == 0: continue
        #Eloss_total.append(Eloss)
        #h['Mufilter'].Fill(Eloss, weight)
    print(Eloss_total.values())
#    for key in h:
#        classname = h[key].Class().GetName()
#        if 'TH' in classname or 'TP' in classname:
#            h[key].Write()
    f.Close()
    Array =  np.array(list(Eloss_total.values()))
    print(Array.shape)
    for name, row in zip(Eloss_total.keys(), Array):
         f = open("output_" + str(name), "w")
         for el in row:
             f.write(str(el) + '\n')
         #f.write('\n')
         f.close()
    #f = open(args.norm + "/flux", "w")
    #f.write(str(B_ids_unw) + "\t" + str(B_ids))
    #f.close()


if __name__ == '__main__':
    r.gErrorIgnoreLevel = r.kWarning
    r.gROOT.SetBatch(True)
    main()
