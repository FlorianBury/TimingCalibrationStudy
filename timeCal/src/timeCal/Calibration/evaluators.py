import os
import sys
import yaml
import glob
import math
import ROOT
import logging
import argparse
import enlighten
import itertools
import datetime
import inspect
from copy import copy,deepcopy
from IPython import embed
import numpy as np

from ..utils.environment import getEnv
from ..utils.context import TFileOpen

from .interpolation import Interpolation
from ..PlotScan.data_helper import Observable

#ROOT.gInterpreter.ProcessLine(f'#include "{os.path.join(path,"morhing.h")}"')

DELTA_COARSE = 1.
DELTA_FINE   = 0.01

CHI2_1SIGMA_PER_NDOF = {1: 1.0, 2: 2.3, 3: 3.53, 4: 4.72, 5: 5.89, 6: 7.04, 7: 8.18, 8: 9.3, 9: 10.42, 10: 11.54}


#*******************************************************#
#                    Base evaluator                     #
#*******************************************************#
class BaseEvaluator:
    #######################################
    #               GLOBAL                #
    #######################################
    _evaluatorName = 'Base'
    #######################################
    #              USE PART               #
    #######################################
    def __init__(self,params,objDict):
        self.params     = params
        self.objDict    = objDict

    def __call__(self,hist,return_graphs=False):

        decimals = lambda x : str(x)[::-1].find('.')

        oneSigNLL = self.getOneSigmaNLL()

        # Coarse scan #
        logging.debug(f"Starting coarse scan : [0.,50.] with step {DELTA_COARSE:0.5f}")
        delayCoarse = np.arange(0.,50.,DELTA_COARSE).round(decimals(DELTA_COARSE))
        NLLCoarse = self.scan(hist,delayCoarse)
        NLLCoarse -= NLLCoarse.min()

        # Find limits of fine scan #
        if len(delayCoarse[NLLCoarse<=oneSigNLL]) > 0:
            delayStart = max(0.,delayCoarse[NLLCoarse<=oneSigNLL].min()-DELTA_COARSE)
            delayStop  = min(50.,delayCoarse[NLLCoarse<=oneSigNLL].max()+DELTA_COARSE)
        else:
            idxMins = np.argpartition(delayCoarse,1)[0:2]
            delayStart = min(delayCoarse[idxMins])
            delayStop  = max(delayCoarse[idxMins])
        if delayStart == delayStop:
            delayStart = delayCoarse[NLLCoarse.argmin()]-DELTA_COARSE
            delayStop  = delayCoarse[NLLCoarse.argmin()]+DELTA_COARSE

        # Fine scan #
        logging.debug(f"Starting fine scan : [{delayStart:0.5f},{delayStop:0.5f}] with step {DELTA_FINE:0.5f}")
        delayFine = np.arange(delayStart,delayStop,DELTA_FINE).round(decimals(DELTA_FINE))
        NLLFine = self.scan(hist,delayFine)
        NLLFine -= NLLFine.min()

        gCoarse = ROOT.TGraph(delayCoarse.shape[0],delayCoarse,NLLCoarse)
        gFine = ROOT.TGraph(delayFine.shape[0],delayFine,NLLFine)
            
        gCoarse.SetTitle(";Test delay [ns];#chi^{2} value")
        gFine.SetTitle(";Test delay [ns];#chi^{2} value")
        numerical_min = delayFine[np.argmin(NLLFine)]
        logging.debug(f"Found minimum at {numerical_min}")
        x_minus, x_plus = self._findSigma(gFine,oneSigNLL)
        logging.debug(f"Best value : {numerical_min} - {x_minus} + {x_plus}")
        if return_graphs:
            return numerical_min,x_minus,x_plus,gCoarse,gFine
        else:
            return numerical_min,x_minus,x_plus


    @staticmethod
    def _findSigma(graph,oneSigNLL):
        x = np.linspace(graph.GetX()[0],graph.GetX()[graph.GetN()-1],10000)
        y = np.array([graph.Eval(xi) for xi in x])
        left_idx  = np.arange(0,np.where(y==y.min())[0][0])
        right_idx = np.arange(np.where(y==y.min())[0][-1],y.shape[0])

        x_left  = x[np.argmin(abs(y[left_idx]-oneSigNLL))] if left_idx.shape[0]>0 else 0.
        x_right = x[np.argmin(abs(y[right_idx]-oneSigNLL)) + right_idx[0]] if right_idx.shape[0]>0 else 0.

        x_min = x[y.argmin()]

        return x_min-x_left,x_right-x_min


    def scan(self,recoHist,delayRange):
        raise NotImplemented

    def getOneSigmaNLL(self):
        raise NotImplemented

    #######################################
    #             SAVE PART               #
    #######################################
    def save(self,subdir):
        #Open root file #
        f = self.__class__.get_evaluator_path(subdir,self.params)
        path = os.path.dirname(f)
        if not os.path.exists(path):
            logging.warning(f'Output path {path} absent, will create it')
            os.makedirs(path)

        # Check in case another run is saving the evaluator #
        if not self.__class__.check_evaluator_file(subdir,self.params):
            with TFileOpen(f,'w') as F:
                # Save parameters in the file #
                self._saveParameters(F)
                # Save the objects to the file #
                self.saveObjects(F)
                # Close file #
                F.Close()

    def _saveParameters(self,F):
        for pName,pVal in self.params.items():
            if isinstance(pVal,str):
                p = ROOT.TNamed(pName,pVal)
            else:
                p = ROOT.TParameter(type(paramVal))(pName,pVal)
            p.Write()

    def saveObjects(self,F):
        raise NotImplementedError

    #######################################
    #             BUILD PART              #
    #######################################
    @classmethod
    def build(cls,F,hists,paramNames):
        # Get file parameters #
        params = cls.getParameters(F,paramNames)
        # Get dict of histograms #
        delayDict = deepcopy(cls.getHistograms(F,hists))

        # Pass to method to create evaluator #
        objDict = cls.makeEvaluator(delayDict)
        # Make object instance #
        instance = cls(params,objDict)
        ## Save #
        #instance.save()

        # Return 
        return instance

    @staticmethod
    def makeEvaluator(delayDict):
        """
            delayDict [dict] = {delay value [float] : BX ID histogram [TH1F]}
        """
        raise NotImplementedError

    @staticmethod
    def getParameters(F,parameters):
        params = {}
        for pFile,pName in parameters.items():
            key = F.GetListOfKeys().FindObject(pFile)
            if not key:
                raise RuntimeError(f'Could not find parameter {pFile} in file {F}')
            if 'TParameter' in key.GetClassName():
                params[pName] = F.Get(key.GetName()).GetVal()
            if 'TNamed' in key.GetClassName():
                params[pName] = F.Get(key.GetName()).GetTitle()
        return params

    @classmethod
    def getHistograms(cls,F,hists):
        histDict = {}
        for name,values in hists.items():
            histDict[float(values['delay'])] = cls._getHistogram(F,values['dir'],name)
        return histDict

    @staticmethod
    def _getHistogram(F,dirPath,name):
        # Get through subdirectories by recurrence #
        subdir = F
        for sub in dirPath.split('/'):
            if not hasattr(subdir,sub):
                raise RuntimeError(f'Could not find subdirectory {sub.GetName()} from {subdir.GeName()} in file {F.GetName()}')
            subdir = getattr(subdir,sub)
        if not subdir.GetListOfKeys().FindObject(name):
            raise RuntimeError(f'Could not find histogram {name} from subdir {dirPath} in file {F.GetName()}')
        return deepcopy(subdir.Get(name))



    #######################################
    #              LOAD PART              #
    #######################################
    @classmethod
    def load(cls,subdir,params):
        f = cls.get_evaluator_path(subdir,params)
        # Open file #
        if cls.check_evaluator_file(subdir,params):
            with TFileOpen(f,'r') as F:
                ## Get parameters #
                #params = cls.getParameters(F,self.parameters)
                # Load the objects # 
                objDict = deepcopy(cls.loadObjects(F))
            # Return object #
            return cls(params,objDict)
        else:
            raise RuntimeError(f'File {f} does not exist')
                

    @staticmethod
    def loadObjects(F):
        raise NotImplementedError


    #######################################
    #                UTILS                #
    #######################################

    @classmethod
    def get_evaluator_path(cls,subdir,params):
        return os.path.join(getEnv()['paths']['calibration'],
                            'evaluators',
                            cls._evaluatorName,
                            subdir,
                            "_".join([f"{pName}_{pVal}" for pName,pVal in params.items()]).replace('.','p').replace(' ','_')+'.root')
                            

    @classmethod
    def check_evaluator_file(cls,subdir,params):
        path = cls.get_evaluator_path(subdir,params)
        return os.path.exists(path)


    def makeIllustrationPlot(self):
        """
            Produce an illustration of the method
            Returns : illustration object
        """
        raise NotImplementedError

    @classmethod
    def makeIllustrationPlotStandalone(cls,delayDict):
        """
            Produce an illustration of the method from standalone (only the delay dict)
            -> can be obtained even without reconstructing the evaluator
            Returns : illustration object
        """
        raise NotImplementedError


#*******************************************************#
#                    Mean evaluator                     #
#*******************************************************#
            
class MeanEvaluator(BaseEvaluator):
    _evaluatorName = 'Mean'

    @staticmethod
    def makeEvaluator(delayDict):
        """
            delayDict [dict] = {delay value [float] : BX ID histogram [TH1F]}
        """
        mGraph_nom   = ROOT.TGraph(len(delayDict))
        mGraph_error = ROOT.TGraph(len(delayDict))
        for i,(delay,hist) in enumerate(sorted(delayDict.items())):
            mean = hist.GetMean() + 0.5 # shift of 0.5 to center the bins
            error = hist.GetMeanError()
            mGraph_nom.SetPoint(i,delay,mean)
            mGraph_error.SetPoint(i,delay,hist.GetMeanError())

        objDict = {"mean": mGraph_nom,"mean_error":mGraph_error}

        return objDict

    def saveObjects(self,F):
        for name,obj in self.objDict.items():
            self.objDict[name].Write(name)

    @staticmethod
    def loadObjects(F):
        return {'mean': deepcopy(F.Get('mean')),'mean_error': deepcopy(F.Get('mean_error'))}

    def getOneSigmaNLL(self):
        return 1.


    def scan(self,recoHist,delayRange):
        NLLs = []
        mean_reco = recoHist.GetMean() + 0.5
        mean_reco_err = recoHist.GetMeanError()
        for i in range(delayRange.shape[0]):
            mean_true = self.objDict['mean'].Eval(delayRange[i])
            mean_true_error = self.objDict['mean_error'].Eval(delayRange[i])
            NLLs.append(abs(mean_reco-mean_true)/math.sqrt(mean_reco_err**2+mean_true_error**2+1e-20))
        return np.array(NLLs)
    
    @classmethod
    def _editIllustrationPlot(cls,graph):
        graph.GetHistogram().SetTitle('')
        graph.GetHistogram().GetXaxis().SetTitle('Delay [ns]')
        graph.GetHistogram().GetYaxis().SetTitle('Mean of BX histogram')
        graph.SetLineWidth(2)
        return graph
    
    def makeIllustrationPlot(self):
        return self.__class__._editIllustrationPlot(self.objDict['mean'])

    @classmethod
    def makeIllustrationPlotStandalone(cls,delayDict):
        objDict = cls.makeEvaluator(delayDict)
        graph = objDict['mean']
        return cls._editIllustrationPlot(graph)
        

#*******************************************************#
#           Linear interpolation evaluator              #
#*******************************************************#
            
class LinearInterpolationEvaluator(BaseEvaluator):
    _evaluatorName = 'LinearInterpolation'

    @staticmethod
    def makeEvaluator(delayDict):
        """
            delayDict [dict] = {delay value [float] : BX ID histogram [TH1F]}
        """
        return {'delayDict':delayDict}

    def saveObjects(self,F):
        for delay,hist in self.objDict['delayDict'].items():
            hist.Write(f'hist_{delay}')

    @staticmethod
    def loadObjects(F):
        objDict = {'delayDict':{}}
        for key in F.GetListOfKeys():
            if key.GetName().startswith('hist'):
                name = key.GetName()
                h = F.Get(key.GetName())
                objDict['delayDict'][float(name.replace('hist_',''))] = deepcopy(h)
        return objDict

    @staticmethod
    def _getNLL(f,d):
        """
            f = expected counts histogram
            d = measured counts histogram
        """
        NLL_Poisson = 0.
        NLL_Gaussian = 0.
        epsilon = 1e-9
        for i in range(1,f.GetNbinsX()+1):
            fi = f.GetBinContent(i)
            di = d.GetBinContent(i)
            #if fi > 0:# and di > 0:
            NLL_Poisson += 2 * (fi-di + di* math.log((di+epsilon)/(fi+epsilon)))
            NLL_Gaussian += (di-fi)**2/(fi+epsilon)
        return NLL_Poisson,NLL_Gaussian
                

    def getOneSigmaNLL(self):
        return 1.
        #ndof = sum([h.GetBinContent(i)>0 for i in range(1,h.GetNbinsX()+1)]) - 1 
        #if ndof not in CHI2_1SIGMA_PER_NDOF.keys():
        #    raise RuntimeError(f'Ndof {ndof} not computed')
        #return CHI2_1SIGMA_PER_NDOF[ndof]

    def scan(self,recoHist,delayRange):
        NLLs = []
        # Loop over the delay range required #
        for delay in delayRange:
            if delay in self.objDict['delayDict'].keys():
                # Hist already in content #
                trueHist = self.objDict['delayDict'][delay]
            else:
                # Find two closest point #
                #delays = np.sort(np.array(list(self.objDict['delayDict'].keys())))
                delays = np.array(list(self.objDict['delayDict'].keys()))
                closest_delay = delays[np.argmin(abs(delays-delay))]
                if closest_delay < delay:
                    # Search on the right #
                    remaining_delays = delays[delays>closest_delay]
                    next_closest_delay = remaining_delays[np.argmin(abs(remaining_delays-delay))]
                    delayLeft  = closest_delay
                    delayRight = next_closest_delay
                else:
                    # Search on the left #
                    remaining_delays = delays[delays<closest_delay]
                    next_closest_delay = remaining_delays[np.argmin(abs(remaining_delays-delay))]
                    delayLeft  = next_closest_delay
                    delayRight = closest_delay
                assert delayLeft <= delay and delay <= delayRight
                interpolation = Interpolation(delayLeft,delayRight,delay)
                trueHist = interpolation(self.objDict['delayDict'][delayLeft],self.objDict['delayDict'][delayRight],f'hist_{delay}')
                self.objDict['delayDict'][delay] = trueHist

            #print ('------')
            #print ([trueHist.GetBinContent(i) for i in range(1,12)])
            #print (delay,recoHist.Chi2Test(trueHist,"UU CHI2"))
            #NLLs.append(recoHist.Chi2Test(trueHist,"UU CHI2/NDF"))
            NLL_Poisson,NLL_Gaussian = self._getNLL(trueHist,recoHist)
#            if 2*abs(NLL_Poisson-NLL_Gaussian) / (NLL_Poisson+NLL_Gaussian + 1e-9) > 0.01:
#                logging.warning(f'In {self.__class__.__name__} scan, NLL(Poisson) = {NLL_Poisson} != NLL(Gaussian) = {NLL_Gaussian} -> Approximation might not stand')
            #print (delay,NLL_Poisson,NLL_Gaussian,2*recoHist.Chi2Test(trueHist,"UU CHI2"))
            NLLs.append(NLL_Poisson)
            #NLLs.append(NLL_Gaussian)
            #NLLs.append(2*recoHist.Chi2Test(trueHist,"UU CHI2/NDF"))
        #embed()

        return np.array(NLLs)

    @classmethod
    def _make2DPlot(cls,delayDict):
        delays = np.array(list(delayDict.keys()))
        delays.sort()
        h_dummy = delayDict[delays[0]]
        x_edges = np.array([h_dummy.GetXaxis().GetBinLowEdge(i) for i in range(1,h_dummy.GetNbinsX()+2)],dtype=np.float32)
        diffs = np.diff(delays)/2
        y_edges = Observable._getEdgesFromCenters(delays).astype(np.float32)
        h2D = ROOT.TH2F('h2D','h2D',
                        len(x_edges)-1,
                        x_edges,
                        len(y_edges)-1,
                        y_edges)
        for iy,delay in enumerate(delays,1):
            h = delayDict[delay]
            for ix in range(1,h.GetNbinsX()+1):
                h2D.SetBinContent(ix,iy,h.GetBinContent(ix))
        h2D.GetXaxis().SetTitle('BX')
        h2D.GetYaxis().SetTitle('Delay [ns]')
        h2D.GetZaxis().SetTitle('Events')
        return h2D
    
    def makeIllustrationPlot(self):
        return self.__class__._make2DPlot(self.objDict['delayDict'])

    @classmethod
    def makeIllustrationPlotStandalone(cls,delayDict):
        return cls._make2DPlot(delayDict)
        


