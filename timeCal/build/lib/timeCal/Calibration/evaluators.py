import os
import sys
import yaml
import glob
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
        hist.Scale(1/hist.Integral()) 

        decimals = lambda x : str(x)[::-1].find('.')

        # Coarse scan #
        logging.debug(f"Starting coarse scan : [0.,50.] with step {DELTA_COARSE:0.5f}")
        delayCoarse = np.arange(0.,50.,DELTA_COARSE).round(decimals(DELTA_COARSE))
        chi2Coarse = self.scan(hist,delayCoarse)
        chi2Coarse -= chi2Coarse.min()
        # TODO -> scan between 1 and 1 limits in chi2

        # Get two minimums #
        idxMin = min(delayCoarse.shape[0]-2,max(1,chi2Coarse.argmin()))
        #idxTwoMins = np.sort(np.argpartition(chi2Coarse,2)[:2])
        x1 = delayCoarse[idxMin-1]
        x2 = delayCoarse[idxMin]
        x3 = delayCoarse[idxMin+1]
        y1 = chi2Coarse[idxMin-1]
        y2 = chi2Coarse[idxMin]
        y3 = chi2Coarse[idxMin+1]
        denom = (x1 - x2)*(x1 - x3)*(x2 - x3)
        a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
        b = (x3**2 * (y1 - y2) + x2**2 * (y3 - y1) + x1**2 * (y2 - y3)) / denom
        c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom
        if a <= 0.: # parabola should be pointed upwards
            xmin = delayCoarse[chi2Coarse.argmin()]
        else:
            xmin = min(50.,max(0.,-b/(2*a)))

        logging.debug(f"Found three points around minimum : [{x1:0.5f},{x2:0.5f},{x3:0.5f}] with values [{y1:0.5f},{y2:0.5f},{y3:0.5f}]")
        logging.info(f"Coarse estimated minimum : {xmin:0.5f}")

        # Fine scan #
        delayStart = max(0.,xmin-DELTA_COARSE)
        delayStop  = min(50.,xmin+DELTA_COARSE)
        logging.debug(f"Starting fine scan : [{delayStart:0.5f},{delayStop:0.5f}] with step {DELTA_FINE:0.5f}")
        delayFine = np.arange(delayStart,delayStop,DELTA_FINE).round(decimals(DELTA_FINE))
        chi2Fine = self.scan(hist,delayFine)
        chi2Fine -= chi2Fine.min()

        gCoarse = ROOT.TGraph(delayCoarse.shape[0],delayCoarse,chi2Coarse)
        gFine = ROOT.TGraph(delayFine.shape[0],delayFine,chi2Fine)
            
        gCoarse.SetTitle(";Test delay [ns];#chi^{2} value")
        gFine.SetTitle(";Test delay [ns];#chi^{2} value")
        numerical_min = delayFine[np.argmin(chi2Fine)]
        sigma = self._findSigma(gFine)
        logging.info("Numerical = {}, fit = {} +/- {}".format(numerical_min,fit_min,sigma))
        if return_graphs:
            return numerical_min,sigma,gCoarse,gFine
        else:
            return numerical_min,sigma


    @staticmethod
    def _findSigma(graph):
        # TODO : find range of X in which chi2 <= 1
        raise NotImplementedError

    def scan(self,recoHist,delayRange):
        raise NotImplementedError

    #######################################
    #             SAVE PART               #
    #######################################
    def save(self):
        #Open root file #
        f = self.__class__.get_evaluator_path(self.params)
        path = os.path.dirname(f)
        if not os.path.exists(path):
            logging.warning(f'Output path {path} absent, will create it')
            os.makedirs(path)

        # Check in case another run is saving the evaluator #
        if not self.__class__.check_evaluator_file(self.params):
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
        # Save #
        instance.save()

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
    def load(cls):
        f = cls.get_evaluator_path(self.params)
        # Open file #
        if cls.check_evaluator_file():
            with TFileOpen(f,'r') as F:
                # Get parameters #
                params = cls.getParameters(F,self.parameters)
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
    def get_evaluator_path(cls,params):
        return os.path.join(getEnv()['paths']['calibration'],
                            'evaluators',
                            cls._evaluatorName,
                            "_".join([f"{pName}_{pVal}" for pName,pVal in params.items()]).replace('.','p'))

    @classmethod
    def check_evaluator_file(cls,params):
        path = cls.get_evaluator_path(params)
        return os.path.exists(path)




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
        mGraph = ROOT.TGraph(len(delayDict))
        for i,(delay,hist) in enumerate(sorted(delayDict.items())):
            mGraph.SetPoint(i,delay,hist.GetMean())

        objDict = {"mean": mGraph}

        return objDict

    def saveObjects(self,F):
        self.objDict['mean'].Write('mean')

    @staticmethod
    def loadObjects(F):
        return {'mean': F.Get('mean')}


    def scan(self,recoHist,delayRange):
        mGraph = self.objDict['mean']
        chi2 = []
        mean_reco = recoHist.GetMean()
        mean_reco_err = recoHist.GetMeanError()
        if mean_reco_err == 0.:
            mean_reco_err = 1.
        for i in range(delayRange.shape[0]):
            mean_true = mGraph.Eval(delayRange[i])
            chi2.append((mean_reco-mean_true)**2/mean_reco_err)
        return np.array(chi2)


#*******************************************************#
#                  Morphing evaluator                   #
#*******************************************************#
class MorphingEvaluator(BaseEvaluator):
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

    _evaluatorName = 'Morphing'

    def __init__(self,params,objDict):
        super().__init__(params,objDict)
        workspace       = self.objDict['workspace']
        self.delayVar   = workspace.var("delayVar")
        self.bxidVar    = workspace.var("bxidVar")
        self.morphing   = workspace.pdf("morph")

        self.cache_hist = {}

        self.morphing.useHorizontalMorphing(True)

        # Get Bx binning info for the histogram creation #
        bxIdBinning = self.bxidVar.GetBinning()
        self.binCenters = []
        self.binEdges = []
        for ibin in range(bxIdBinning.numBins()):
            self.binCenters.append(bxidBinning.binCenter(ibin))
            self.binEdges.append(bxidBinning.binLow(ibin))
        self.binEdges.append(bxidBinning.binHigh(bxidBinning.numBins() - 1))

        self.normSet = ROOT.RooArgSet(self.bxidVar)
            
    @staticmethod
    def makeEvaluator(delayDict):
        """
            delayDict [dict] = {delay value [float] : BX ID histogram [TH1F]}
        """

        bxidVar = ROOT.RooRealVar("bxidVar","bxidVar",-5.5,5.5)
        delayVar = ROOT.RooRealVar("delayVar","delayVar",min(list(delayDict.keys())),max(list(delayDict.keys())))
        bxidVar.setBins(11)
        delayVar.setBins(len(delayDict))
        paramVec = ROOT.TVectorD(len(delayDict))

        listOfMorphs = ROOT.RooArgList("listOfMorphs")
        listPdfs = []
        listHist = []

        for i,delay in enumerate(sorted(list(delayDict.keys()))):
            h = delayDict[delay]
            paramVec[i] = delay 
            # Delta for convergence #
            h.Scale(1./h.Integral())
            delta = 1e-3
            idx = np.array([h.GetBinContent(i) for i in range(1,h.GetNbinsX()+1)]).argmax()+1
            if idx > 0:
                h.SetBinContent(int(idx-1),h.GetBinContent(int(idx-1))+delta)
            if idx > 1:
                h.SetBinContent(int(idx-2),h.GetBinContent(int(idx-2))+delta)
            if idx < h.GetNbinsX():
                h.SetBinContent(int(idx+1),h.GetBinContent(int(idx+1))+delta)
            if idx < h.GetNbinsX()-1:
                h.SetBinContent(int(idx+2),h.GetBinContent(int(idx+2))+delta)

            hD = ROOT.RooDataHist(h.GetName()+"DataHist",h.GetName()+"DataHist",ROOT.RooArgList(bxidVar),h)
            listHist.append(h)
            hPdf = ROOT.RooHistPdf(h.GetName()+"Pdf",h.GetName()+"Pdf",ROOT.RooArgSet(bxidVar),hD)
            listPdfs.append(deepcopy(hPdf))

        for pdf in listPdfs:
            listOfMorphs.add(pdf)

        morph = ROOT.RooMomentMorph('morph','morph',
                                    delayVar,
                                    ROOT.RooArgList(bxidVar),
                                    listOfMorphs,
                                    paramVec,
                                    ROOT.RooMomentMorph.Linear)
                                    #ROOT.RooMomentMorph.SineLinear)

        h2D = morph.createHistogram("test", 
                                    bxidVar, 
                                    ROOT.RooFit.Binning(11),
                                    ROOT.RooFit.YVar(delayVar))

        w = ROOT.RooWorkspace("MorphWorkspace","MorphWorkspace")
        getattr(w,'import')(morph)

        objDict = {"Shape2D": h2D,'workspace':w}

        return objDict
            
    def saveObjects(self,F):
        self.objDict['Shape2D'].Write('Shape2D')
        self.objDict['workspace'].writeToFile(F.GetName(),True) # True = recreate

    @staticmethod
    def loadObjects(F):
        return {'workspace': F.Get("MorphWorkspace")}


    def scan(self,recoHist,delayRange):
        chi2 = []
        count = 0
        # Loop over the delay range required #
        for i in range(delayRange.shape[0]):
            # Check cache #
            if delayRange[i] not in self.cache_hist.keys():
                trueHist = self._makeBXhistogram(delayRange[i]) #self.morphing.createHistogram("trueHist{:.5f}".format(delayRange[i]),self.bxidVar)
                self.cache_hist[delayRange[i]] = trueHist
                count += 1
            else:
                trueHist = self.cache_hist[delayRange[i]]
            chi2.append(recoHist.Chi2Test(trueHist,"WW NORM CHI2"))

        # Debug printout of the cache size #
        logging.debug('In dict {}, new elements {} / {}'.format(len(self.cache_hist),count,delayRange.shape[0]))
        size = sys.getsizeof(self.cache_hist)
        for k,v in self.cache_hist.items():
            size += sys.getsizeof(k)
            size += sys.getsizeof(v)
        logging.debug('Size = ',size/1e6)
        return np.array(chi2)

    def _makeBXhistogram(self,delay):
        self.delayVar.setVal(delayRange[i])
        trueHist = ROOT.TH1D(f"trueHist{delay:.5f}",f"trueHist{delay:.5f}",len(self.binCenters),np.array(self.binEdges))
        for ibi,binCenter in enumerate(self.binCenters,1):
            self.bxidVar.SetVal(binCenter)
            trueHist.SetBinContent(ibin,morphing.getVal(self.normSet))
        
        #trueHist = self.morphing.createHistogram(f"trueHist{delay:.5f}",self.bxidVar) # Slow 
        return trueHist
  

#if __name__ == "__main__":
#        parser = argparse.ArgumentParser(description='Produce template morphings')
#        parser.add_argument('--yaml', action='store', required=True, type=str, 
#                            help='Yaml containing parameters')
#        parser.add_argument('--submit', action='store', required=False, type=str, default=None,
#                            help='Name for submission')
#        parser.add_argument('--split', action='store', required=False, type=int, default=None,
#                            help='Number of files per job')
#        parser.add_argument('--debug', action='store_true', required=False, default=False,
#                            help='Debug for submit jobs')
#        args = parser.parse_args()
#
#        if args.submit is not None and args.split is not None:
#            from utils import submit
#            submit(args.yaml,args.submit,args.split,['*.root'],['files'],args.debug)
#        else:
#            clsmembers = {clsName:cls for clsName,cls in inspect.getmembers(sys.modules[__name__], inspect.isclass)}
#            with open(args.yaml,'r') as handle:
#                f = yaml.load(handle,Loader=yaml.FullLoader)
#            evaluator = f['evaluator']
#            if evaluator not in clsmembers.keys():
#                raise RuntimeError(f"Evaluator {evaluator} class not present in this file")
#            cls = clsmembers[evaluator]
#            instance = cls(**f)
#

