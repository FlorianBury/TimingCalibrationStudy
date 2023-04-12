import os
import sys
import re
import yaml
import math
import glob
import ROOT
import argparse
import enlighten
import itertools
import datetime
from pprint import pprint
from copy import copy,deepcopy
from IPython import embed
import numpy as np

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
ROOT.TH1.AddDirectory(False)

DELTA_COARSE = 1.
DELTA_FINE = 0.01

class UseEvaluator:
    def __init__(self,filepath):
        if not os.path.exists(filepath):
            raise RuntimeError(f"File {filepath} does not exist")
        self.loadFile(filepath)

    @property
    def name(self):
        raise NotImplementedError

    def loadFile(self,filepath):
        if not os.path.exists(filepath):
            raise RuntimeError("Could not find file {}".format(filepath))
        F = ROOT.TFile(filepath)
        self.params = self._getParameters(F)
        return F

    @staticmethod
    def _getParameters(F):
        params = {}
        for key in F.GetListOfKeys():
            if 'TParameter' in key.GetClassName():
                params[key.GetName()] = F.Get(key.GetName()).GetVal()
            if 'TNamed' in key.GetClassName():
                params[key.GetName()] = F.Get(key.GetName()).GetTitle()
        return params

    def __call__(self,hist,return_graphs=False,verbose=False):
        hist.Scale(1/hist.Integral()) 

        decimals = lambda x : str(x)[::-1].find('.')

        # Coarse scan #
        if verbose:
            print (f"Starting coarse scan : [0.,50.] with step {DELTA_COARSE:0.5f}")
        bxCoarse = np.arange(0.,50.,DELTA_COARSE).round(decimals(DELTA_COARSE))
        chi2Coarse = self.scan(hist,bxCoarse)
        chi2Coarse -= chi2Coarse.min()
        # TODO -> scan between 1 and 1 limits in chi2

        # Get two minimums #
        idxMin = min(bxCoarse.shape[0]-2,max(1,chi2Coarse.argmin()))
        #idxTwoMins = np.sort(np.argpartition(chi2Coarse,2)[:2])
        x1 = bxCoarse[idxMin-1]
        x2 = bxCoarse[idxMin]
        x3 = bxCoarse[idxMin+1]
        y1 = chi2Coarse[idxMin-1]
        y2 = chi2Coarse[idxMin]
        y3 = chi2Coarse[idxMin+1]
        denom = (x1 - x2)*(x1 - x3)*(x2 - x3)
        a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
        b = (x3**2 * (y1 - y2) + x2**2 * (y3 - y1) + x1**2 * (y2 - y3)) / denom
        c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom
        if a <= 0.: # parabola should be pointed upwards
            xmin = bxCoarse[chi2Coarse.argmin()]
        else:
            xmin = min(50.,max(0.,-b/(2*a)))

        if verbose:
            print (f"Found three points around minimum : [{x1:0.5f},{x2:0.5f},{x3:0.5f}] with values [{y1:0.5f},{y2:0.5f},{y3:0.5f}]")
            print (f"Coarse estimated minimum : {xmin:0.5f}")

        # Fine scan #
        bxStart = max(0.,xmin-DELTA_COARSE)
        bxStop  = min(50.,xmin+DELTA_COARSE)
        if verbose:
            print (f"Starting fine scan : [{bxStart:0.5f},{bxStop:0.5f}] with step {DELTA_FINE:0.5f}")
        bxFine = np.arange(bxStart,bxStop,DELTA_FINE).round(decimals(DELTA_FINE))
        chi2Fine = self.scan(hist,bxFine)
        chi2Fine -= chi2Fine.min()

        gCoarse = ROOT.TGraph(bxCoarse.shape[0],bxCoarse,chi2Coarse)
        gFine = ROOT.TGraph(bxFine.shape[0],bxFine,chi2Fine)
            
        gCoarse.SetTitle(";Test delay [ns];#chi^{2} value")
        gFine.SetTitle(";Test delay [ns];#chi^{2} value")
        numerical_min = bxFine[np.argmin(chi2Fine)]
        sigma = self._findSigma(gFine)
        if verbose:
            print ("Numerical = {}, fit = {} +/- {}".format(numerical_min,fit_min,sigma))
        if return_graphs:
            return numerical_min,sigma,gCoarse,gFine
        else:
            return numerical_min,sigma


    @staticmethod
    def _findSigma(graph):
        # TODO : find range of X in which chi2 <= 1
        raise NotImplementedError
#        graph.Fit("pol2","SQR","",x1,x2)
#        fit = graph.GetListOfFunctions().FindObject("pol2")
#        c = fit.GetParameter(0)
#        b = fit.GetParameter(1)
#        a = fit.GetParameter(2)
#        f = lambda x : a*x**2+b*x+c
#        if a <= 0. : 
#            return None,50.
#        xmin = -b/(2*a)
#        ymin = f(xmin)
#        x1 = (-b-math.sqrt(4*a))/(2*a)
#        x2 = (-b+math.sqrt(4*a))/(2*a)
#        assert f(x1)-ymin-1 < 1e-9
#        assert f(x2)-ymin-1 < 1e-9
#        sigma = x2-x1
#        return xmin,sigma

    def scan(self,recoHist,bxRange):
        raise NotImplementedError


class MorphingEvaluator(UseEvaluator):
    def loadFile(self,filepath):
        F = super(MorphingEvaluator,self).loadFile(filepath)
        self.workspace = deepcopy(F.Get("MorphWorkspace"))
        self.delayVar = self.workspace.var("delayVar")
        self.bxidVar = self.workspace.var("bxidVar")
        self.morphing = self.workspace.pdf("morph")
        self.shape = deepcopy(F.Get("Shape2D"))
        self.morphing.useHorizontalMorphing(True)
        F.Close()

        self.cache_hist = {}

    @property
    def name(self):
        return "Morphing method"

    def scan(self,recoHist,bxRange):
        chi2 = []
        count = 0
        for i in range(bxRange.shape[0]):
            # Check cache #
            if bxRange[i] not in self.cache_hist.keys():
                self.delayVar.setVal(bxRange[i])
                trueHist = self.morphing.createHistogram("trueHist{:.5f}".format(bxRange[i]),self.bxidVar)
                self.cache_hist[bxRange[i]] = trueHist
                count += 1
            else:
                trueHist = self.cache_hist[bxRange[i]]
            chi2.append(recoHist.Chi2Test(trueHist,"WW NORM CHI2"))
        print ('In dict {}, new elements {} / {}'.format(len(self.cache_hist),count,bxRange.shape[0]))
        size = sys.getsizeof(self.cache_hist)
        for k,v in self.cache_hist.items():
            size += sys.getsizeof(k)
            size += sys.getsizeof(v)
        print ('Size = ',size/1e6)
        return np.array(chi2)

class MeanEvaluator(UseEvaluator):
    def loadFile(self,filepath):
        F = super().loadFile(filepath)
        self.mGraph = deepcopy(F.Get("Mean"))
        F.Close()

    @property
    def name(self):
        return "Mean method"

    def scan(self,recoHist,bxRange):
        chi2 = []
        mean_reco = recoHist.GetMean()
        mean_reco_err = recoHist.GetMeanError()
        if mean_reco_err == 0.:
            mean_reco_err = 1.
        for i in range(bxRange.shape[0]):
            mean_true = self.mGraph.Eval(bxRange[i])
            chi2.append((mean_reco-mean_true)**2/mean_reco_err)
        return np.array(chi2)

        
