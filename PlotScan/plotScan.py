import os
import sys
import yaml
import glob
import ROOT
import argparse
import enlighten
import itertools
from pprint import pprint
from copy import copy,deepcopy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm

from data_helper import Data, Observable

class PlotScan:
    def __init__(self,path,suffix,files,hists,parameters,efficiency,**kwargs):
        self.path           = path
        self.files          = files
        self.suffix         = suffix
        self.hists          = hists
        self.parameters     = parameters
        self.efficiency     = efficiency
        self.content        = []
        #self.cache          = os.path.join(os.path.abspath(os.path.dirname(__file__)),'cache','cache_{}.npy')
        self.cache          = os.path.join(os.path.abspath(os.path.dirname(__file__)),'cache','cache_{}.pkl'.format(self.suffix))
        self.observables = ['Efficiency','Out-of-time fraction','Mean','Next BX contamination','Previous BX contamination']
        self.limits = {'Efficiency'             : (0.,1.),
                       'Out-of-time fraction'   : (0.,1.),
                       'Mean'                   : (-1.,1.),
                       'Next BX contamination'  : (0.,1.),
                       'Previous BX contamination'  : (0.,1.)}

        self._getRootFiles()

    def _getRootFiles(self):

        if not os.path.exists(self.cache):
            print('Not all cache files are present, will process root files')
            for f in self.files.keys():
                if not os.path.exists(os.path.join(self.path,f)):
                    print ('Could not find file : {}'.format(os.path.join(self.path,f)))
            pbar = enlighten.Counter(total=len(self.files.keys()), desc='Progress', unit='files')
            i = 0
            for f,params in self.files.items():
                pbar.update()
                hDict = self._getHistograms(os.path.join(self.path,f))
                effDict,contNextDict,contPrevDict = self._getEfficiency(os.path.join(self.path,f))

                for h,values in hDict.items():
                    if values['type'] == 'TH1':
                        entry = {k:v for k,v in zip(self.parameters,params)}
                        entry.update(values)
                        entry['Efficiency'] = effDict[values['mode']][values['delay']]
                        entry['Next BX contamination'] = contNextDict[values['mode']][values['delay']]
                        entry['Previous BX contamination'] = contPrevDict[values['mode']][values['delay']]
                        fo,fo_err = self._fracOutHisto(h)
                        entry['Out-of-time fraction'] = fo
                        entry['Out-of-time fraction error'] = fo_err
                        m,m_err = self._meanHisto(h)
                        entry['Mean'] = m
                        entry['Mean error'] = m_err
                        
                        self.content.append(entry)
            df = pd.DataFrame(self.content)
            df.to_pickle(self.cache)
            print ('Saved cache in {}'.format(self.cache))

        else:
            df = pd.read_pickle(self.cache)
            print ('Loaded cache from {}'.format(self.cache))
        self.data = Data(df)
        self.data.SetParameters(self.parameters + ['delay'])

    def _getHistograms(self,f):
        F = ROOT.TFile(f)
        #ROOT.SetOwnership(F,False)
        histDict = {}
        for name,values in self.hists.items():
            h = deepcopy(F.Get(values['dir']+"/"+name))
            histDict[h] = {k:v for k,v in values.items() if k!='dir'}
        F.Close()
        return histDict

    def _getEfficiency(self,f):
        F = ROOT.TFile(f)
        effDict = {mode:{} for mode in self.efficiency.keys()}
        contNextDict = {mode:{} for mode in self.efficiency.keys()}
        contPrevDict = {mode:{} for mode in self.efficiency.keys()}
        for mode,d in self.efficiency.items():
            trueName = d['dir'] + '/' + d['truth']
            recoName = d['dir'] + '/' + d['reco']
            h_true = F.Get(trueName)
            h_reco = F.Get(recoName)
            h_eff = h_reco.Clone("Eff")
            h_eff.Divide(h_true)
            bin_center  = h_eff.GetXaxis().FindBin(0)
            bin_next    = h_eff.GetXaxis().FindBin(1)
            bin_prev    = h_eff.GetXaxis().FindBin(-1)
            for y in range(1,h_eff.GetNbinsY()+1):
                delay = round(h_eff.GetYaxis().GetBinCenter(y),1)
                effDict[mode][delay] = h_eff.GetBinContent(bin_center,y)
                contNextDict[mode][delay] = h_eff.GetBinContent(bin_next,y)
                contPrevDict[mode][delay] = h_eff.GetBinContent(bin_prev,y)

        F.Close()
        return effDict,contNextDict,contPrevDict
        

    @staticmethod
    def _getHistContent(h):
        e = [h.GetXaxis().GetBinUpEdge(0)]
        w = []
        s = []
        for i in range(1,h.GetNbinsX()+1):
            e.append(h.GetXaxis().GetBinUpEdge(i))
            w.append(h.GetBinContent(i))
            s.append(h.GetBinError(i))
        return np.array(e),np.array(w),np.array(s)
 
    def _fracOutHisto(self,h): 
        e,w,s = self._getHistContent(h)
        c = (e[:-1]+e[1:])/2
        idx = np.delete(np.arange(w.shape[0]),np.abs(c).argmin())
        return w[idx].sum()/w.sum(),np.sqrt((s[idx]**2).sum())/w.sum()

    def _meanHisto(self,h): 
        e,w,s = self._getHistContent(h)
        c = (e[:-1]+e[1:])/2
        return (w*c).sum()/w.sum(),np.sqrt(((s*c)**2).sum())/w.sum()
    
    def _makeArray(self,var):
        axNames = self.parameters + ['delay']
        df = self.df[axNames+[var]]
        grouped = df.groupby(axNames)[var].mean()
        shape = tuple(map(len, grouped.index.levels))
        arr = np.full(shape, np.nan)
        arr[tuple(grouped.index.codes)] = grouped.values.flat
        labels = {level.name:level.values for level in grouped.index.levels}
        print ('Produced array for variable {}'.format(var))
        return Data(arr,labels)

    def Plots(self):
        from IPython import embed
        for obsName in self.observables:
            observable = self.data.GetObservable(obsName)
            for param in self.parameters:
                print ("Plotting {} curves of {} parameter".format(obsName,param))
                with PdfPages('pdf/curves_{}_{}_{}.pdf'.format(self.suffix,obsName,param).replace(' ','')) as pdf:
                    labels_to_vary = {key:observable.GetLabels()[key] for key in observable.GetLabels().keys() if key != param and key != 'delay'}

                    for comb in itertools.product(*list(labels_to_vary.values())):
                        varied_labels = {k:c for k,c in zip(labels_to_vary.keys(),comb)}
                        obs2d = observable.GetSlice(**varied_labels)

                        # Plot 1D #
                        fig, ax = plt.subplots(figsize=(8,7))
                        plt.subplots_adjust(left=0.15, right=0.90, top=0.85, bottom=0.12)
                        paramValues = obs2d.GetLabels()[param]
                        colors = cm.jet(np.linspace(0,1,paramValues.shape[0]))
                        for paramVal,color in zip(paramValues,colors):
                            obs1d = obs2d.GetSlice(**{param:paramVal})
                            obs1d.Pyplot1D(ax,color=color)
                        sm = cm.ScalarMappable(cmap=cm.rainbow, norm=plt.Normalize(vmin=paramValues.min(), vmax=paramValues.max()))
                        cbar = plt.colorbar(sm)
                        cbar.set_label(param,fontsize=18,labelpad=20)
                        cbar.ax.tick_params(labelsize=14)
                        plt.xlabel('Delay [ns]',fontsize=18,labelpad=20)
                        plt.ylabel(obsName,fontsize=18,labelpad=10)
                        plt.ylim(self.limits[obsName])
                        title = '{} curves\n('.format(param)+','.join(['{} = {}'.format(p,v) for p,v in varied_labels.items()])+')'
                        plt.title(title,fontsize=20,pad=25)
                        plt.tick_params(axis='both', which='major', labelsize=14)
                        print ('\t'+title)
                        pdf.savefig()
                        plt.close()

                        # Plot 2D #
                        fig,ax = plt.subplots(figsize=(8,7))
                        plt.subplots_adjust(left=0.15, right=0.95, top=0.85, bottom=0.1)
                        obs2d.Pyplot2D(x='delay',y=param, ax=ax, vmin=self.limits[obsName][0], vmax=self.limits[obsName][1], shading='auto', linewidth=0,rasterized=True)
                        cbar = fig.colorbar(ax.collections[0])
                        cbar.set_label(obsName,fontsize=18,labelpad=20)
                        cbar.ax.tick_params(labelsize=14)
                        ax.set_xlabel('Delay [ns]',fontsize=18,labelpad=10)
                        ax.set_ylabel(param,fontsize=18,labelpad=20)
                        ax.set_title(title,fontsize=20,pad=25)
                        ax.tick_params(axis='both', which='major', labelsize=14)
                        pdf.savefig()
                        plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Produce datacards')
    parser.add_argument('--yaml', action='store', required=True, type=str,
                        help='Yaml containing parameters')
    args = parser.parse_args()


    if args.yaml is None:
        raise RuntimeError("Must provide the YAML file")
    with open(args.yaml,'r') as handle:
        f = yaml.load(handle,Loader=yaml.FullLoader)

    instance = PlotScan(**f)
    instance.Plots()


