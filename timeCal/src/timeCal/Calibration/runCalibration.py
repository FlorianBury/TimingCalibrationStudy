import os
import yaml
import glob
import ROOT
import logging
import inspect
import argparse
import enlighten
import numpy as np
import multiprocessing as mp
from IPython import embed
from copy import deepcopy

from . import evaluators as Evaluators

from ..utils.environment import getEnv
from ..utils.yamlLoader import parseYaml
from ..utils.context import TFileOpen

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


from dataclasses import dataclass,field

@dataclass
class Scenario:
    paramNames: dict = None                             # Link between param Names and names for plots
    testParams: dict = None                             # Values of the parameters of the test
    recoParams: dict = None                             # Values of the parameters of the evaluator
    recoConfig: dict = None                             # Evaluator config
    testHists : dict = None                             # Test Histograms {delay:h}
    recoIllust: 'typing.Any' = object()                 # Illustration object of the evaluator method
    testIllust: 'typing.Any' = object()                 # Illustration object of the tests as if there were used to build the evaluator 
    testDelays: list = field(default_factory=list)      # True delay values
    recoDelays: list = field(default_factory=list)      # Reconstructed decays values
    recoErrors: list = field(default_factory=list)      # Errors from reconstruction
    coarseScan: dict = field(default_factory=dict)      # Coarse scan delay value tested
    fineScan  : dict = field(default_factory=dict)      # Fine scan delay value tested
    evaluator : 'typing.Any' = object()                 # Evaluator object

def RunEvaluator(scObj):
    # Obtain class and parameters #
    f = os.path.join(scObj.recoConfig['dir'],scObj.recoConfig['file'])
    if scObj.recoConfig['method'] not in  [module[0] for module in inspect.getmembers(Evaluators)]:
        raise RuntimeError(f"Evaluator {scObj.recoConfig['method']} class not found")

    evaluatorCls = getattr(Evaluators,scObj.recoConfig['method'])
    with TFileOpen(f,'r') as F:
        scObj.recoParams = deepcopy(evaluatorCls.getParameters(F,scObj.paramNames))
        # If evaluator exists, load it #
        if evaluatorCls.check_evaluator_file(scObj.recoParams):
            logging.info(f'Found evaluator {evaluatorCls.get_evaluator_path(scObj.recoParams)}, will load it')
            scObj.evaluator = evaluatorCls.load(scObj.recoParams)
            logging.info('... done')
        # Else, build it #
        else:
            logging.info(f'Evaluator {evaluatorCls.get_evaluator_path(scObj.recoParams)} not found, will build it')
            scObj.evaluator = evaluatorCls.build(F,scObj.recoConfig['hists'],scObj.paramNames)
            scObj.evaluator.save()
            logging.info('... done')

    for testDelay,hist in scObj.testHists.items():
        # Use method to estimate reco delay #
        logging.debug(f'True delay = {testDelay}')
        recoDelay, recoErrorMinus, recoErrorPlus, gCoarse, gFine = scObj.evaluator(hist,return_graphs=True)

        # Save values #
        scObj.testDelays.append(testDelay)
        scObj.recoDelays.append(recoDelay)
        scObj.recoErrors.append((recoErrorMinus,recoErrorPlus))

        # Save scan graphs #
        scObj.coarseScan[testDelay] = gCoarse
        scObj.fineScan[testDelay]   = gFine

    # Create arrays #
    x = np.array(scObj.testDelays,dtype=np.float32)
    y = np.array(scObj.recoDelays,dtype=np.float32)
    e = np.array(scObj.recoErrors,dtype=np.float32)
    e_down   = e[:,0]
    e_up     = e[:,1]
        
    # Create graphs #
    N = x.shape[0]
    scObj.graph = ROOT.TGraphAsymmErrors(N,
                                         x.astype(np.float32),
                                         y.astype(np.float32),
                                         np.zeros(N,dtype=np.float32),
                                         np.zeros(N,dtype=np.float32),
                                         e_down.astype(np.float32),
                                         e_up.astype(np.float32))

    scObj.diff = ROOT.TGraph(N,x,abs(x-y))

    # Create illustration #
    try:
        scObj.recoIllust = scObj.evaluator.makeIllustrationPlot()
        scObj.testIllust = evaluatorCls.makeIllustrationPlotStandalone(scObj.testHists)
    except Exception as e:
        logging.warning(f'Could not produce illustration due to `{e}`')

    return scObj

class RunCalibration:
    def __init__(self,suffix,parameters,scenarios,jobs=None,**kwargs):
        self.suffix = suffix
        self.paramNames = parameters
        self.scenarios = scenarios
        self.jobs = jobs
        self.scObjs = []

        self.run()

        for scObj in self.scObjs:
            self.producePlots(scObj)


    def run(self):
        """
            scenario [dict] :
                 'test' : {
                            'file'    : root file [str],
                            'dir'     : path to the file [str],
                            'params'  : parameters of the file [list[float]],
                            'hists'   : histogram config [dict]
                                         key   : name of the histogram
                                         value : config [dict] 
                                                {'delay': delay [float],'dir': directory [srt]}
                           }
                 'evaluators' [list] : list of evaluator configs [dict]
                  {
                                'method' : name of the class in useEvaluator.py [str],
                                'file'   : file containing evaluator [str],
                                'dir'    : path to the file [str],
                                'params' : parameters of the evaluator [list[float]],                              }
        """

        for scenario in self.scenarios:
            # Get the test content 
            f = os.path.join(scenario['test']['dir'],scenario['test']['file'])
            with TFileOpen(f,'r') as F:
                histDict = deepcopy(Evaluators.BaseEvaluator.getHistograms(F,scenario['test']['hists']))
                params = deepcopy(Evaluators.BaseEvaluator.getParameters(F,self.paramNames))
            # Instantiate the evaluator runs #
            for recoConfig in scenario['evaluators']:
                scObj = Scenario(paramNames   = self.paramNames,
                                   testParams   = params,
                                   recoConfig   = recoConfig,
                                   testHists    = histDict)
                self.scObjs.append(scObj)

        # Start the scenario running #
        pbar = enlighten.Counter(total=len(self.scObjs), desc='Progress', unit='scenarios')
        if self.jobs is None:
            for scObj in self.scObjs:
                RunEvaluator(scObj)
                pbar.update()
        else:
            results = []
            pool = mp.Pool(self.jobs)
            for result in pool.imap(RunEvaluator,self.scObjs):
                results.append(result)
                pbar.update()
            self.scObjs = results

    def producePlots(self,scObj):
        # Make PDF name #
        dirName = os.path.join(getEnv()['paths']['calibration'],'plots',self.suffix)
        if not os.path.exists(dirName):
            os.makedirs(dirName)

        pdfName = os.path.join(dirName,
                               ('Test_' + \
                                    '_'.join([f'{pName}_{pVal}' 
                                            for pName,pVal in scObj.testParams.items()]) + \
                                    f'_{scObj.recoConfig["method"]}_' + \
                                    '_'.join([f'{pName}_{pVal}' 
                                            for pName,pVal in scObj.recoParams.items()])
                               ).replace(' ','_').replace('.','p') + '.pdf')


        # Make canvas #
        C = ROOT.TCanvas('C','C',1000,600)
        C.SetLeftMargin(0.1)
        C.SetRightMargin(0.2)
        C.SetTopMargin(0.05)
        C.SetBottomMargin(0.10)
        C.SetGridx()
        C.SetGridy()
        C.SetTickx()
        C.SetTicky()    
        C.Print(pdfName+'[')

        # Esthetics #
        scObj.graph.SetFillColorAlpha(ROOT.kGreen+2, 0.5)
        scObj.graph.SetLineStyle(1)
        scObj.graph.SetLineWidth(2)
        scObj.graph.SetLineColor(9)

        # Print reco versus true delay #
        C.Clear()
        scObj.graph.GetHistogram().GetXaxis().SetRangeUser(-1,51)
        scObj.graph.GetHistogram().GetYaxis().SetRangeUser(-1,51)
        scObj.graph.GetHistogram().SetTitle("")
        scObj.graph.GetHistogram().GetXaxis().SetTitle("True delay [ns]")
        scObj.graph.GetHistogram().GetYaxis().SetTitle("Reco delay [ns]")

        scObj.graph.Draw("ALE4")

        line = ROOT.TLine(0.,0.,50.,50.)
        line.SetLineWidth(2)
        line.SetLineStyle(2)
        line.Draw()

        # Add the parameters #
        text = ['Histogram parameters :']
        for paramName,paramVal in scObj.testParams.items():
            text.append(f'  {paramName} = {paramVal}')
        text.append('')
        text.append(f'Reco : {scObj.recoConfig["method"]}')
        for paramName,paramVal in scObj.evaluator.params.items():
            text.append(f'  {paramName} = {paramVal}')

        pt = ROOT.TPaveText(.81,.5,.99,.9,"NDC")
        for txtline in text:
            pt.AddText(txtline)    
        pt.SetTextAlign(11)
        pt.SetFillStyle(4001)
        pt.SetBorderSize(0)

        pt.Draw()

        C.Print(pdfName,'Title:True-vs-Reco')

        # Print diff #
        C.Clear()
        scObj.diff.GetHistogram().GetYaxis().SetRangeUser(0.,np.array(scObj.diff.GetY()).max()*1.1)
        scObj.diff.GetHistogram().SetTitle("")
        scObj.diff.GetHistogram().GetXaxis().SetTitle("True delay [ns]")
        scObj.diff.GetHistogram().GetYaxis().SetTitle("|True-reco| delay [ns]")

        scObj.diff.Draw("APL")
        pt.Draw()
        C.Print(pdfName,'Title:True-Reco')

        # Print reco and test illustrations #

        C.Clear()
        if scObj.testIllust.__class__ != object and scObj.recoIllust.__class__ != object:
            recoType = scObj.recoIllust.__class__.__name__
            testType = scObj.testIllust.__class__.__name__
            if recoType != testType:
                raise RuntimeError(f'Illustration of test is of type {testType} while reco is of type {recoType}')
            if recoType == 'TGraph':
                legend = ROOT.TLegend(0.81,0.2,0.99,0.4) 
                legend.AddEntry(scObj.testIllust,'Test evaluator')
                legend.AddEntry(scObj.recoIllust,'Reco evaluator')
                legend.SetFillStyle(4001)
                legend.SetBorderSize(0)
                scObj.testIllust.GetHistogram().GetYaxis().SetRangeUser(
                            min(list(scObj.testIllust.GetY())+list(scObj.recoIllust.GetY())),
                            max(list(scObj.testIllust.GetY())+list(scObj.recoIllust.GetY())))
                scObj.testIllust.SetLineColor(ROOT.kGreen+2)
                scObj.recoIllust.SetLineColor(ROOT.kBlue+2)
                scObj.testIllust.Draw()
                scObj.recoIllust.Draw('same')    
                legend.Draw()
            elif recoType.startswith('TH2'):
                mainpad = ROOT.TPad('mainpad','mainpad',0,0,0.75,1)
                mainpad.Divide(2)
                mainpad.Draw()
                pad1 = mainpad.cd(1)
                pad1.SetLogz()
                scObj.testIllust.SetTitle('Test evaluator')
                scObj.testIllust.Draw('colz')
                pad2 = mainpad.cd(2)
                pad2.SetLogz()
                scObj.recoIllust.SetTitle('Reco evaluator')
                scObj.recoIllust.Draw('colz')
                C.cd()
            else:
                raise RuntimeError(f'Type {recoType} not understood')

            pt.Draw()

            C.Print(pdfName,'Title:Method illustration')
        else:
            logging.warning(f'Will not produce comparative illustration plot as they were not produced ')


        # Print scans #
        for idx,testDelay in enumerate(scObj.testDelays):
            # Get associated reco values #
            recoDelay  = scObj.recoDelays[idx]
            recoErrors = scObj.recoErrors[idx]
            # Canvas #
            C.Clear()
            C.Divide(2)
            # Coarse scan #
            c1 = C.cd(1)
            c1.SetLogy(True)
            scObj.coarseScan[testDelay].SetTitle(f"Coarse scan (true delay = {testDelay});Delay [ns];-2#Deltaln(L)")
            scObj.coarseScan[testDelay].GetHistogram().GetYaxis().SetRangeUser(
                    min(list(scObj.coarseScan[testDelay].GetY())) * 0.1,
                    max(list(scObj.coarseScan[testDelay].GetY())) * 10)
            scObj.coarseScan[testDelay].Draw()
            # Fine scan #
            c2 = C.cd(2)
            scObj.fineScan[testDelay].SetTitle(f"Fine scan (true delay = {testDelay});Delay [ns];-2#Deltaln(L)")
            fineX = np.array(scObj.fineScan[testDelay].GetX())
            fineY = np.array(scObj.fineScan[testDelay].GetY())
            scObj.fineScan[testDelay].GetHistogram().GetXaxis().SetRangeUser(
                        fineX[fineY<=2.].min(),
                        fineX[fineY<=2.].max())
            scObj.fineScan[testDelay].GetHistogram().GetYaxis().SetRangeUser(0.,2.)
            scObj.fineScan[testDelay].Draw()
            # Draw 1 sigma line #
            lines = [
                ROOT.TLine(recoDelay-recoErrors[0],1.,recoDelay+recoErrors[1],1.),
                ROOT.TLine(recoDelay,0.,recoDelay,1.),
                ROOT.TLine(recoDelay-recoErrors[0],0.,recoDelay-recoErrors[0],1.),
                ROOT.TLine(recoDelay+recoErrors[1],0.,recoDelay+recoErrors[1],1.)]
            for line in lines:
                line.SetLineColor(ROOT.kRed+1)
                line.Draw()
            # Print values #
            pt = ROOT.TPaveText(.5,.8,.9,.9,"NDC")
            pt.AddText(f"delay = {recoDelay:.3f}_{{-{recoErrors[0]:.3f}}}^{{+{recoErrors[1]:.3f}}}")
            pt.SetTextAlign(11)
            pt.SetFillStyle(4001)
            pt.SetBorderSize(0)
            pt.Draw()
            # Print #
            C.Print(pdfName,f'Title:Scan_{testDelay}')

        # Final save #
        C.Print(pdfName+']')


        # Printout #
        logging.info(f'Produced plot in {pdfName}')
        
def main():
    parser = argparse.ArgumentParser(description='Run calibration')
    parser.add_argument('--yaml', action='store', required=True, type=str, 
                        help='Yaml containing parameters')
    parser.add_argument('--custom', action='store', required=False, default=None, nargs='*',
                        help='Format the yaml file')
    parser.add_argument('--submit', action='store', required=False, type=str, default=None,
                        help='Name for submission')
    parser.add_argument('--split', action='store', required=False, type=int, default=None,
                        help='Number of files per job')
    parser.add_argument('-v','--verbose', action='store_true', required=False, default=False,
                        help='Verbose mode')
    parser.add_argument('-j','--jobs', action='store', required=False, default=None, type=int,
                        help='Number of jobs for multiprocessing')

    args = parser.parse_args()

    logging.basicConfig(level   = logging.DEBUG,
                        format  = '%(asctime)s - %(levelname)s - %(message)s',
                        datefmt = '%m/%d/%Y %H:%M:%S')

    # Verbose level #
    if not args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    
    if args.submit is not None and args.split is not None:
        from utils import submit
        submit(args.yaml,args.submit,args.split,['*.pdf'],['scenarios'],args.jobs)
    else:
        f = parseYaml(args.yaml,args.custom)
        RunCalibration(**f,jobs=args.jobs)

if __name__ == "__main__":
    main()
