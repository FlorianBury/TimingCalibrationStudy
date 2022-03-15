import os
import yaml
import glob
import ROOT
import inspect
import argparse
import enlighten
import numpy as np
from IPython import embed
from copy import deepcopy

import useEvaluator
import buildEvaluator

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

clsmembers = {clsName:cls for clsName,cls in  inspect.getmembers(useEvaluator, inspect.isclass)}

class RunEvaluator:
    def __init__(self,parameters,scenarios,**kwargs):
        self.parameters = parameters
        for i,scenario in enumerate(scenarios,1):
            print (f'Running scenario {i}/{len(scenarios)} [{i*100./len(scenarios):3.2f}%]')
            self._runScenario(scenario)
        
    def _runScenario(self,scenario):
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

        # Initialize evaluators #
        evaluators = []
        for evaluatorConfig in scenario['evaluators']:
            if evaluatorConfig['method'] not in clsmembers.keys():
                raise RuntimeError(f"Evaluator {evaluatorConfig['method']} class not found")
            evaluator = getattr(useEvaluator,evaluatorConfig['method'])(os.path.join(evaluatorConfig['dir'],evaluatorConfig['file']))
            evaluators.append(evaluator)

        # Loop over histograms #
        F = ROOT.TFile(os.path.join(scenario['test']['dir'],scenario['test']['file']))
        file_params = self._orderParameters(self._getParameters(F))
        true_delays = []
        reco_delays = [[] for _ in evaluators] 
        xminus_delays = [[] for _ in evaluators]
        xplus_delays = [[] for _ in evaluators]
        delayDict = {}
        chi2_delays = [{} for _ in evaluators] 

        pbar = enlighten.Counter(total=len(scenario['test']['hists']), desc='Progress', unit='hist')
        for histName,histConfig in scenario['test']['hists'].items():
            true_delays.append(histConfig['delay'])
            hist = deepcopy(F.Get(histConfig['dir'] + "/" + histName))
            delayDict[histConfig['delay']] = hist
            for i,evaluator in enumerate(evaluators):
                reco_delay, reco_minus, reco_plus, gCoarse, gFine = evaluator(hist,return_graphs=True,verbose=False)
#                if histConfig['delay'] == 13.8:
#                    ROOT.gROOT.SetBatch(True)
#                    embed()
#                if sigma_delay is None:
#                    sigma_delay = 0.
                reco_delays[i].append(reco_delay)
                xminus_delays[i].append(reco_minus)
                xplus_delays[i].append(reco_plus)
                chi2_delays[i][histConfig['delay']] = (gCoarse, gFine)
            pbar.update()
        F.Close()
 
        # Produce TGraphs #
        n = len(true_delays)
        x = np.array(true_delays)
        for i,evaluator in enumerate(evaluators):
            y      = np.array(reco_delays[i])
            e_down = np.array(xminus_delays[i])
            e_up   = np.array(xplus_delays[i])
            y_down = y - e_down
            y_up   = y + e_up

            evaluator_params = self._orderParameters(evaluator.params)
            test_cls = evaluator.__class__.__name__
            graphs = {}
            graph_test = None
            graph_reco = None
            if test_cls == 'MeanEvaluator':
                graph_test = getattr(buildEvaluator,test_cls).makeEvaluator(delayDict)['Mean']
                graph_reco = evaluator.mGraph
#            if test_cls == 'MorphingEvaluator':
#                graph_test = getattr(buildEvaluator,test_cls).makeEvaluator(delayDict)['Shape2D']
#                graph_reco = evaluator.shape
#                graph_test.SetTitle("True;BX ID;delay [ns]")
#                graph_reco.SetTitle("Reco;BX ID;delay [ns]")

            # Text #
            text = ['Histogram parameters :']
            for paramName,paramVal in self._translateParameters(file_params).items():
                text.append(f'  {paramName} = {paramVal}')
            text.append('')
            text.append(f'Reco : {evaluator.name}')
            for paramName,paramVal in self._translateParameters(evaluator_params).items():
                text.append(f'  {paramName} = {paramVal}')

            pt = ROOT.TPaveText(.81,.5,.99,.9,"NDC")
            for line in text:
                pt.AddText(line)    
            pt.SetTextAlign(11)
            pt.SetFillStyle(4001)
            pt.SetBorderSize(0)

            # Plot comparison #
            pdfName = f"PDF/{test_cls}_Test_{self._formatParameters(file_params)}_Evaluator_{self._formatParameters(evaluator_params)}.pdf"
            background, g_central, g_band = self._getBands(x,y,y_up,y_down)
                
            background.GetXaxis().SetRangeUser(-1,51)
            background.GetYaxis().SetRangeUser(-1,51)
            background.GetXaxis().SetTitle("True delay [ns]")
            background.GetYaxis().SetTitle("Reco delay [ns]")

            C = ROOT.TCanvas("C","C",1000, 600)
            C.SetLeftMargin(0.1)
            C.SetRightMargin(0.2)
            C.SetTopMargin(0.05)
            C.SetBottomMargin(0.10)
            C.SetGridx()
            C.SetGridy()
            C.SetTickx()
            C.SetTicky()    

            C.Print(pdfName+'[')

            # Central + bands #
            background.Draw()
            g_band.Draw("fe3same")
            g_central.Draw("lpsame")
            pt.Draw()

            # Diagonal #
            line = ROOT.TLine(0.,0.,50.,50.)
            line.SetLineStyle(2)
            line.Draw("same")

            C.Print(pdfName,'Title:Reco versus True delay')

            # Plot delta #
            C.Clear()
            delta      = x-y
            delta_down = delta - e_down
            delta_up   = delta + e_up
            background, g_central, g_band = self._getBands(x,delta,delta_up,delta_down)
            background.GetXaxis().SetRangeUser(-1,51)
            background.GetYaxis().SetRangeUser(-3.,3.)
            background.GetXaxis().SetTitle("True delay [ns]")
            background.GetYaxis().SetTitle("True-reco delay [ns]")

            background.Draw()
            g_band.Draw("fe3same")
            g_central.Draw("lpsame")
            pt.Draw()

            C.Print(pdfName,'Title:True-Reco')


            # Additional graphs #
            if graph_test is not None or graph_reco is not None:
                C.Clear()
                if isinstance(graph_test,ROOT.TGraph):
                    maxy = max([g.GetHistogram().GetMaximum() for g in [graph_test,graph_reco]])
                    miny = max([g.GetHistogram().GetMinimum() for g in [graph_test,graph_reco]])
                    background.GetYaxis().SetRangeUser(miny*0.9,maxy*1.1)
                    background.Draw()
                    graph_test.SetLineWidth(2)
                    graph_test.SetLineColor(635)
                    graph_reco.SetLineWidth(2)
                    graph_reco.SetLineStyle(2)
                    graph_reco.SetLineColor(603)
                    leg = ROOT.TLegend(.81,.2,.99,.4)
                    leg.SetFillStyle(4001)
                    leg.SetBorderSize(0)
                    leg.AddEntry(graph_test,'True')
                    leg.AddEntry(graph_reco,'Reco')
                    graph_test.Draw("same")
                    graph_reco.Draw("same")
                    leg.Draw()
                if isinstance(graph_test,ROOT.TH2):
                    C.Divide(2)
                    c1 = C.cd(1)
                    graph_test.Draw("colz")
                    c2 = C.cd(2)
                    graph_reco.Draw("colz")

                pt.Draw()
                C.Print(pdfName,'Title:Method plots')

            # Chi2 graphs #
            for j,(delay,(gCoarse,gFine)) in enumerate(chi2_delays[i].items()):
                C.Clear()
                C.Divide(2)
                c1 = C.cd(1)
                c1.SetLogy(True)
                gCoarse.SetTitle(f"Coarse scan (true delay = {delay})")
                gCoarse.Draw()
                c2 = C.cd(2)
                gFine.SetTitle(f"Fine scan (true delay = {delay})")
                gFine.GetHistogram().GetYaxis().SetRangeUser(0.,2.)
                gFine.Draw()
                # Draw 1 sigma line #
                line = ROOT.TLine(y_down[j],1.,y_up[j],1.)
                # Print values #
                pt = ROOT.TPaveText(.5,.8,.9,.9,"NDC")
                pt.AddText(f"delay = {y[j]:.3f}^{{+{e_up[j]:.3f}}}_{{-{e_down[j]:.3f}}}")
                pt.SetTextAlign(11)
                pt.SetFillStyle(4001)
                pt.SetBorderSize(0)
                pt.Draw()


                line.Draw()

                C.Print(pdfName,f'Title:Chi2 (delay = {delay})')
            
            
            C.Print(pdfName+']')
            del C
            
    def _getBands(self,x,central,up,down):
        # Create band #
        band = np.concatenate((up,down[::-1]))
        x_all = np.concatenate((x,x[::-1]))

        # Create graphs #
        g_central = ROOT.TGraph(x.shape[0], x, central)
        g_band = ROOT.TGraph(x_all.shape[0], x_all, band)

        g_band.SetFillColorAlpha(ROOT.kGreen+2, 0.5)
        g_band.SetLineStyle(1)
        g_central.SetLineWidth(2)
        g_central.SetLineColor(9)


        # Plot #
#        minx = x.min()-1
#        maxx = x.max()+1
        background = ROOT.TH1F("b","", x.shape[0]*2, -1.,51.)
#        background.SetTitle(" ")
#        background.GetXaxis().SetTitle("True delay [ns]")
#        background.GetYaxis().SetRangeUser(minx,maxx)
#        background.GetYaxis().SetTitle("Reco delay [ns]")
#        background.SetStats(0)

        return background, g_central, g_band


    def _getParameters(self,F):
        params = {}
        for key in F.GetListOfKeys():
            if key.GetName() not in self.parameters:
                continue
            if 'TParameter' in key.GetClassName():
                params[key.GetName()] = F.Get(key.GetName()).GetVal() 
            if 'TNamed' in key.GetClassName():
                params[key.GetName()] = F.Get(key.GetName()).GetTitle()
        return params

    def _orderParameters(self,params):
        return {inFileName:params[inFileName] for inFileName in self.parameters.keys()}

    def _translateParameters(self,params):
        return {paramName:params[inFileName] for inFileName,paramName in self.parameters.items()}

    @staticmethod
    def _formatParameters(params):
        return "_".join([f"{paramName}_{params[paramName]}" for paramName,paramVal in params.items()]).replace('.','p')

  
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Produce datacards')
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
        submit(args.yaml,args.submit,args.split,['*.pdf'],['scenarios'],args.debug)
    else:
        with open(args.yaml,'r') as handle:
            f = yaml.load(handle,Loader=yaml.FullLoader)
        RunEvaluator(**f)
