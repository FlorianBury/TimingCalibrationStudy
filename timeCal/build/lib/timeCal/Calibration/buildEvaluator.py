import os
import sys
import yaml
import glob
import ROOT
import argparse
import enlighten
import itertools
import datetime
import inspect
from copy import copy,deepcopy
from IPython import embed
import numpy as np

class BuildEvaluator:
    def __init__(self,path,files,hists,parameters,**kwargs):
        self.path           = path
        self.files          = files
        self.hists          = hists
        self.parameters     = parameters
        
        self._loopOverRootFiles()

    def _loopOverRootFiles(self):
        for f in self.files:
            if not os.path.exists(os.path.join(self.path,f)):
                print ('Could not find file : {}'.format(f))
        pbar = enlighten.Counter(total=len(self.files), desc='Progress', unit='files')
        for f in self.files:
            F = ROOT.TFile(os.path.join(self.path,f))
            histDict = self._getHistograms(F)
            params = self._getParameters(F)
            F.Close()
            for mode, delayDict in histDict.items():
                paramDict = {'mode':mode,**params}
                objDict = self.makeEvaluator(delayDict)
                filename = self._makeFileName(paramDict)
                self._fileWriter(filename,paramDict,objDict)
            pbar.update()

    def _getHistograms(self,F):
        histDict = {}
        for name,values in self.hists.items():
            h = deepcopy(F.Get(values['dir']+"/"+name))
            if values['mode'] not in histDict.keys():
                histDict[values['mode']] = {}
            histDict[values['mode']][values['delay']] = h
        return histDict

    def _getParameters(self,F):
        params = {}
        for key in F.GetListOfKeys():
            if key.GetName() not in self.parameters:
                continue
            if 'TParameter' in key.GetClassName():
                params[key.GetName()] = F.Get(key.GetName()).GetVal()
            if 'TNamed' in key.GetClassName():
                params[key.GetName()] = F.Get(key.GetName()).GetTitle()
        return {paramName:params[paramName] for paramName in self.parameters} # Reorder them

    @property
    def evaluatorName(self):
        raise NotImplementedError
        
    def _makeFileName(self,paramDict):
        """
            paramDict [dict] = {parameter name [str] : parameter value [float]}
        """
        name = self.evaluatorName
        for k,v in paramDict.items():
            if isinstance(v,float):
                name += ("_%s_%0.1f"%(k,v)).replace('.','p')
            elif isinstance(v,int):
                name += "_%s_%d"%(k,v)
            elif isinstance(v,str):
                name += ("_%s_%s"%(k,v)).replace('.','p')
            else:
                raise RuntimeError("Param value type not understood")
        name += ".root"
        return name

    def _fileWriter(self,filename,paramDict,objDict):
        """
            filename [str]   = name of the root file to be created
            paramDict [dict] = {parameter name [str] : parameter value [float]}
            objDict [dict]   = {object name [str] : object [ROOT.*]}
        """
        if 'workspace' in objDict.keys():
            objDict['workspace'].writeToFile(filename,True) # True = recreate 
            
        F = ROOT.TFile(filename,"UPDATE")
        for paramName,paramVal in paramDict.items():
            if isinstance(paramVal,str):
                p = ROOT.TNamed(paramName,paramVal)
            else:
                p = ROOT.TParameter(type(paramVal))(paramName,paramVal)
            p.Write()
        for objName, obj in objDict.items():
            if objName == 'workspace':
                continue
            obj.Write(objName)
        F.Close()


    @staticmethod
    def makeEvaluator(delayDict):
        """
            delayDict [dict] = {delay value [float] : BX ID histogram [TH1F]}
        """
        raise NotImplementedError
   

class MorphingEvaluator(BuildEvaluator):
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

    @property
    def evaluatorName(self):
        return "Morphing"

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
            
class MeanEvaluator(BuildEvaluator):
    @property
    def evaluatorName(self):
        return "Mean"

    @staticmethod
    def makeEvaluator(delayDict):
        """
            delayDict [dict] = {delay value [float] : BX ID histogram [TH1F]}
        """
        mGraph = ROOT.TGraph(len(delayDict))
        for i,(delay,hist) in enumerate(sorted(delayDict.items())):
            mGraph.SetPoint(i,delay,hist.GetMean())

        objDict = {"Mean": mGraph}

        return objDict
   

if __name__ == "__main__":
        parser = argparse.ArgumentParser(description='Produce template morphings')
        parser.add_argument('--yaml', action='store', required=True, type=str, 
                            help='Yaml containing parameters')
        parser.add_argument('--submit', action='store', required=False, type=str, default=None,
                            help='Name for submission')
        parser.add_argument('--split', action='store', required=False, type=int, default=None,
                            help='Number of files per job')
        parser.add_argument('--debug', action='store_true', required=False, default=False,
                            help='Debug for submit jobs')
        args = parser.parse_args()

        if args.submit is not None and args.split is not None:
            from utils import submit
            submit(args.yaml,args.submit,args.split,['*.root'],['files'],args.debug)
        else:
            clsmembers = {clsName:cls for clsName,cls in inspect.getmembers(sys.modules[__name__], inspect.isclass)}
            with open(args.yaml,'r') as handle:
                f = yaml.load(handle,Loader=yaml.FullLoader)
            evaluator = f['evaluator']
            if evaluator not in clsmembers.keys():
                raise RuntimeError(f"Evaluator {evaluator} class not present in this file")
            cls = clsmembers[evaluator]
            instance = cls(**f)


