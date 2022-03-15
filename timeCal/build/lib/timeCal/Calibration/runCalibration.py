import os
import yaml
import glob
import ROOT
import logging
import inspect
import argparse
import enlighten
import numpy as np
from IPython import embed
from copy import deepcopy

from . import evaluators as Evaluators

from ..utils.environment import getEnv
from ..utils.yamlLoader import parseYaml
from ..utils.context import TFileOpen

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


class RunEvaluator:
    def __init__(self,paramNames,testParams,evalConfig,testDict):
        # Save the test content #
        self.testParams = testParams
        self.testDict   = testDict

        # Obtain class and parameters #
        f = os.path.join(evalConfig['dir'],evalConfig['file'])
        if evalConfig['method'] not in  [module[0] for module in inspect.getmembers(Evaluators)]:
            raise RuntimeError(f"Evaluator {evalConfig['method']} class not found")

        evaluatorCls = getattr(Evaluators,evalConfig['method'])
        with TFileOpen(f,'r') as F:
            self.evalParams = deepcopy(evaluatorCls.getParameters(F,paramNames))

            # If evaluator exists, load it #
            if evaluatorCls.check_evaluator_file(self.evalParams):
                logging.info(f'Found evaluator {evaluatorCls.get_evaluator_path(self.evalParams)}, will load it')
                self.evaluator = evaluatorCls.load()
                logging.info('... done')
            # Else, build it #
            else:
                logging.info(f'Evaluator {evaluatorCls.get_evaluator_path(self.evalParams)} not found, will build it')
                self.evaluator = evaluatorCls.build(F,evalConfig['hists'],paramNames)
                logging.info('... done')


            
    def run(self):
        self.trueDelays = [] 
        self.recoDelays = []
        self.recoErrors = []
        for trueDelay,hist in self.testDict.items():
            # Use method to estimate reco delay #
            recoDelay, recoErrorPlus, recoErrorMinus = self.evaluator(hist)

            # Save values #
            self.trueDelays.append(trueDelay)
            self.recoDelays.append(recoDelay)
            self.recoErrors.append((recoErrorPlus, recoErrorMinus))

        self.produceRecoVsTrueDelayPlot()

    def produceRecoVsTrueDelayPlot(self):
        # Create arrays #
        x = np.array(self.trueDelays)
        y = np.array(self.recoDelays)
        e = np.array(self.recoErrors)
        e_up   = e[:,0]
        e_down = e[:,1]
            
        # Create bands #
        band = np.concatenate((e_up,e_down[::-1]))
        x_all = np.concatenate((x,x[::-1]))

        # Create graphs #
        self.g_central = ROOT.TGraph(x.shape[0], x, central)
        self.g_band = ROOT.TGraph(x_all.shape[0], x_all, band)

        self.g_band.SetFillColorAlpha(ROOT.kGreen+2, 0.5)
        self.g_band.SetLineStyle(1)
        self.g_central.SetLineWidth(2)
        self.g_central.SetLineColor(9)

        # Plot #
#        minx = x.min()-1
#        maxx = x.max()+1
        background = ROOT.TH1F("b","", x.shape[0]*2, -1.,51.)
#        background.SetTitle(" ")
#        background.GetXaxis().SetTitle("True delay [ns]")
#        background.GetYaxis().SetRangeUser(minx,maxx)
#        background.GetYaxis().SetTitle("Reco delay [ns]")
#        background.SetStats(0)
        
    def printPlots(self,pdfName):
        C = ROOT.TCanvas('C','C',1000,600)
        C.SetLeftMargin(0.1)
        C.SetRightMargin(0.2)
        C.SetTopMargin(0.05)
        C.SetBottomMargin(0.10)
        C.SetGridx()
        C.SetGridy()
        C.SetTickx()
        C.SetTicky()    

        # Print reco versus true delay #
        background = ROOT.TH1F("b","", len(self.trueDelays)*2, -1.,51.)
        background.GetXaxis().SetRangeUser(-1,51)
        background.GetYaxis().SetRangeUser(-1,51)
        background.GetXaxis().SetTitle("True delay [ns]")
        background.GetYaxis().SetTitle("Reco delay [ns]")

        background.Draw() 
        self.g_central.Draw("fe3same")
        self.g_central.Draw("lpsame")

        # Add the parameters #
        pt = ROOT.TPaveText(.81,.5,.99,.9,"NDC")
        for line in text:
            pt.AddText(line)    
        pt.SetTextAlign(11)
        pt.SetFillStyle(4001)
        pt.SetBorderSize(0)

        text = ['Histogram parameters :']
        for paramName,paramVal in self.testParams.items():
            text.append(f'  {paramName} = {paramVal}')
        text.append('')
        text.append(f'Reco : {evaluator.name}')
        for paramName,paramVal in self.evaluator.params.items():
            text.append(f'  {paramName} = {paramVal}')
        pt.Draw()


class RunCalibration:
    def __init__(self,suffix,parameters,scenarios,**kwargs):
        self.suffix = suffix
        self.paramNames = parameters
        self.scenarios = scenarios

        self.run()


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

        self.instances = []
        for scenario in self.scenarios:
            # Get the test content 
            f = os.path.join(scenario['test']['dir'],scenario['test']['file'])
            with TFileOpen(f,'r') as F:
                histDict = deepcopy(Evaluators.BaseEvaluator.getHistograms(F,scenario['test']['hists']))
                params = deepcopy(Evaluators.BaseEvaluator.getParameters(F,self.paramNames))
            # Instantiate the evaluator runs #
            for evalConfig in scenario['evaluators']:
                self.instances.append(RunEvaluator(self.paramNames,
                                                   params,
                                                   evalConfig,
                                                   histDict))

        # Start the scenario running #
        pbar = enlighten.Counter(total=len(self.instances), desc='Progress', unit='scenarios')
        for instance in self.instances:
            instance.run()
            pbar.update()
        
#    def producePlots(self,pdfName):
#        for instance in self.instances:
            


        

#class RunEvaluator:
#    def __init__(self,parameters,evalConfig,histDict,**kwargs):
#        self.parameters = parameters
#        self.evalConfig = evalConfig
#        self.histDict = histDict
##        for i,scenario in enumerate(scenarios,1):
##            print (f'Running scenario {i}/{len(scenarios)} [{i*100./len(scenarios):3.2f}%]')
##            self._runScenario(scenario)
#        
#    def run(self):
#        # Initialize evaluator #
#
#
#        evaluators = []
#        for evalConfig in scenario['evaluators']:
#            if evalConfig['method'] not in clsmembers.keys():
#            evaluator = getattr(useEvaluator,evalConfig['method'])(os.path.join(evalConfig['dir'],evalConfig['file']))
#            evaluators.append(evaluator)
#
#        # Loop over histograms #
#        F = ROOT.TFile(os.path.join(scenario['test']['dir'],scenario['test']['file']))
#        file_params = self._orderParameters(self._getParameters(F))
#        true_delays = []
#        reco_delays = [[] for _ in evaluators] 
#        xminus_delays = [[] for _ in evaluators]
#        xplus_delays = [[] for _ in evaluators]
#        delayDict = {}
#        chi2_delays = [{} for _ in evaluators] 
#
#        pbar = enlighten.Counter(total=len(scenario['test']['hists']), desc='Progress', unit='hist')
#        for histName,histConfig in scenario['test']['hists'].items():
#            true_delays.append(histConfig['delay'])
#            hist = deepcopy(F.Get(histConfig['dir'] + "/" + histName))
#            delayDict[histConfig['delay']] = hist
#            for i,evaluator in enumerate(evaluators):
#                reco_delay, reco_minus, reco_plus, gCoarse, gFine = evaluator(hist,return_graphs=True,verbose=False)
#                reco_delays[i].append(reco_delay)
#                xminus_delays[i].append(reco_minus)
#                xplus_delays[i].append(reco_plus)
#                chi2_delays[i][histConfig['delay']] = (gCoarse, gFine)
#            pbar.update()
#        F.Close()
# 
#        # Produce TGraphs #
#        n = len(true_delays)
#        x = np.array(true_delays)
#        for i,evaluator in enumerate(evaluators):
#            y      = np.array(reco_delays[i])
#            e_down = np.array(xminus_delays[i])
#            e_up   = np.array(xplus_delays[i])
#            y_down = y - e_down
#            y_up   = y + e_up
#
#            evaluator_params = self._orderParameters(evaluator.params)
#            test_cls = evaluator.__class__.__name__
#            graphs = {}
#            graph_test = None
#            graph_reco = None
#            if test_cls == 'MeanEvaluator':
#                graph_test = getattr(buildEvaluator,test_cls).makeEvaluator(delayDict)['Mean']
#                graph_reco = evaluator.mGraph
##            if test_cls == 'MorphingEvaluator':
##                graph_test = getattr(buildEvaluator,test_cls).makeEvaluator(delayDict)['Shape2D']
##                graph_reco = evaluator.shape
##                graph_test.SetTitle("True;BX ID;delay [ns]")
##                graph_reco.SetTitle("Reco;BX ID;delay [ns]")
#
#            # Text #
#            text = ['Histogram parameters :']
#            for paramName,paramVal in self._translateParameters(file_params).items():
#                text.append(f'  {paramName} = {paramVal}')
#            text.append('')
#            text.append(f'Reco : {evaluator.name}')
#            for paramName,paramVal in self._translateParameters(evaluator_params).items():
#                text.append(f'  {paramName} = {paramVal}')
#
#            pt = ROOT.TPaveText(.81,.5,.99,.9,"NDC")
#            for line in text:
#                pt.AddText(line)    
#            pt.SetTextAlign(11)
#            pt.SetFillStyle(4001)
#            pt.SetBorderSize(0)
#
#            # Plot comparison #
#            pdfName = f"PDF/{test_cls}_Test_{self._formatParameters(file_params)}_Evaluator_{self._formatParameters(evaluator_params)}.pdf"
#            background, g_central, g_band = self._getBands(x,y,y_up,y_down)
#                
#            background.GetXaxis().SetRangeUser(-1,51)
#            background.GetYaxis().SetRangeUser(-1,51)
#            background.GetXaxis().SetTitle("True delay [ns]")
#            background.GetYaxis().SetTitle("Reco delay [ns]")
#
#            C = ROOT.TCanvas("C","C",1000, 600)
#            C.SetLeftMargin(0.1)
#            C.SetRightMargin(0.2)
#            C.SetTopMargin(0.05)
#            C.SetBottomMargin(0.10)
#            C.SetGridx()
#            C.SetGridy()
#            C.SetTickx()
#            C.SetTicky()    
#
#            C.Print(pdfName+'[')
#
#            # Central + bands #
#            background.Draw()
#            g_band.Draw("fe3same")
#            g_central.Draw("lpsame")
#            pt.Draw()
#
#            # Diagonal #
#            line = ROOT.TLine(0.,0.,50.,50.)
#            line.SetLineStyle(2)
#            line.Draw("same")
#
#            C.Print(pdfName,'Title:Reco versus True delay')
#
#            # Plot delta #
#            C.Clear()
#            delta      = x-y
#            delta_down = delta - e_down
#            delta_up   = delta + e_up
#            background, g_central, g_band = self._getBands(x,delta,delta_up,delta_down)
#            background.GetXaxis().SetRangeUser(-1,51)
#            background.GetYaxis().SetRangeUser(-3.,3.)
#            background.GetXaxis().SetTitle("True delay [ns]")
#            background.GetYaxis().SetTitle("True-reco delay [ns]")
#
#            background.Draw()
#            g_band.Draw("fe3same")
#            g_central.Draw("lpsame")
#            pt.Draw()
#
#            C.Print(pdfName,'Title:True-Reco')
#
#
#            # Additional graphs #
#            if graph_test is not None or graph_reco is not None:
#                C.Clear()
#                if isinstance(graph_test,ROOT.TGraph):
#                    maxy = max([g.GetHistogram().GetMaximum() for g in [graph_test,graph_reco]])
#                    miny = max([g.GetHistogram().GetMinimum() for g in [graph_test,graph_reco]])
#                    background.GetYaxis().SetRangeUser(miny*0.9,maxy*1.1)
#                    background.Draw()
#                    graph_test.SetLineWidth(2)
#                    graph_test.SetLineColor(635)
#                    graph_reco.SetLineWidth(2)
#                    graph_reco.SetLineStyle(2)
#                    graph_reco.SetLineColor(603)
#                    leg = ROOT.TLegend(.81,.2,.99,.4)
#                    leg.SetFillStyle(4001)
#                    leg.SetBorderSize(0)
#                    leg.AddEntry(graph_test,'True')
#                    leg.AddEntry(graph_reco,'Reco')
#                    graph_test.Draw("same")
#                    graph_reco.Draw("same")
#                    leg.Draw()
#                if isinstance(graph_test,ROOT.TH2):
#                    C.Divide(2)
#                    c1 = C.cd(1)
#                    graph_test.Draw("colz")
#                    c2 = C.cd(2)
#                    graph_reco.Draw("colz")
#
#                pt.Draw()
#                C.Print(pdfName,'Title:Method plots')
#
#            # Chi2 graphs #
#            for j,(delay,(gCoarse,gFine)) in enumerate(chi2_delays[i].items()):
#                C.Clear()
#                C.Divide(2)
#                c1 = C.cd(1)
#                c1.SetLogy(True)
#                gCoarse.SetTitle(f"Coarse scan (true delay = {delay})")
#                gCoarse.Draw()
#                c2 = C.cd(2)
#                gFine.SetTitle(f"Fine scan (true delay = {delay})")
#                gFine.GetHistogram().GetYaxis().SetRangeUser(0.,2.)
#                gFine.Draw()
#                # Draw 1 sigma line #
#                line = ROOT.TLine(y_down[j],1.,y_up[j],1.)
#                # Print values #
#                pt = ROOT.TPaveText(.5,.8,.9,.9,"NDC")
#                pt.AddText(f"delay = {y[j]:.3f}^{{+{e_up[j]:.3f}}}_{{-{e_down[j]:.3f}}}")
#                pt.SetTextAlign(11)
#                pt.SetFillStyle(4001)
#                pt.SetBorderSize(0)
#                pt.Draw()
#
#
#                line.Draw()
#
#                C.Print(pdfName,f'Title:Chi2 (delay = {delay})')
#            
#            
#            C.Print(pdfName+']')
#            del C

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
    parser.add_argument('--debug', action='store_true', required=False, default=False,
                        help='Debug for submit jobs')
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
        submit(args.yaml,args.submit,args.split,['*.pdf'],['scenarios'],args.debug)
    else:
        f = parseYaml(args.yaml,args.custom)
        RunCalibration(**f)

if __name__ == "__main__":
    main()
