import glob
import os
import json
import sys
import csv
import math as m
import numpy as np
import ROOT as r

r.gROOT.SetBatch(True)

from optparse import OptionParser

def parse_arguments():
    usage = """ 
    usage: %prog [options] Calculates expected pseudosignificance
    """

    parser = OptionParser(usage=usage)

    parser.add_option(  "-i", "--input",
                        help="""Input .root file""",
                        dest = "input",
                        default = "outputFromSeedEff.root"
                    )

    (options, args) = parser.parse_args()
    return options, args

def createDir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def main(options, paths):

    r.gStyle.SetOptStat(0)
    nSimHits = r.TH1F("nSimHits",";Number of BPIX Layers with Hits",6,0,6) 
    nSimHitsL1 = r.TH1F("nSimHitsL1",";Number of BPIX Hits in Layer1 ",8,0,8) 
    nSimHitsL2 = r.TH1F("nSimHitsL2",";Number of BPIX Hits in Layer2",8,0,8) 
    nSimHitsL3 = r.TH1F("nSimHitsL3",";Number of BPIX Hits in Layer3",8,0,8) 
    nSimHitsL4 = r.TH1F("nSimHitsL4",";Number of BPIX Hits in Layer4",8,0,8)  
    nSimHitsFPIX = r.TH1F("nSimHitsFPIX",";Number of FPIX Layers with Hits",6,0,6) 
    nSimHitsL1FPIX = r.TH1F("nSimHitsL1FPIX",";Number of FPIX Hits in Layer1 ",8,0,8) 
    nSimHitsL2FPIX = r.TH1F("nSimHitsL2FPIX",";Number of FPIX Hits in Layer2",8,0,8)  
    nSimHitsL3FPIX = r.TH1F("nSimHitsL3FPIX",";Number of FPIX Hits in Layer3",8,0,8)  
    nSimHitsL4FPIX = r.TH1F("nSimHitsL4FPIX",";Number of FPIX Hits in Layer4",8,0,8)
    nSimHitsTotal= r.TH1F("nSimHitsTotal",";Number of Total Layers with Hits",9,0,9) 

    simElepT = r.TH1F("simPt",";p_{T} [GeV]",20,0.,100.) 
    simEleEta = r.TH1F("simEta",";#eta",20,-3.,3.) 
    simElePhi = r.TH1F("simPhi",";#phi",21,-3.15,3.15) 

    simZpeak = r.TH1F("simMass",";mass [GeV]",20,50.,150.)

    nRecoHits = r.TH1F("nRecoHits",";Number of BPIX Layers with Hits",6,0,6) 
    nRecoHitsL1 = r.TH1F("nRecoHitsL1",";Number of BPIX Hits in Layer1 ",8,0,8) 
    nRecoHitsL2 = r.TH1F("nRecoHitsL2",";Number of BPIX Hits in Layer2",8,0,8)
    nRecoHitsL3 = r.TH1F("nRecoHitsL3",";Number of BPIX Hits in Layer3",8,0,8) 
    nRecoHitsL4 = r.TH1F("nRecoHitsL4",";Number of BPIX Hits in Layer4",8,0,8) 
    nRecoHitsFPIX= r.TH1F("nRecoHitsFPIX",";Number of FPIX Layers with Hits",6,0,6)
    nRecoHitsL1FPIX = r.TH1F("nRecoHitsL1FPIX",";Number of FPIX Hits in Layer1 ",8,0,8)
    nRecoHitsL2FPIX = r.TH1F("nRecoHitsL2FPIX",";Number of FPIX Hits in Layer2",8,0,8)
    nRecoHitsL3FPIX = r.TH1F("nRecoHitsL3FPIX",";Number of FPIX Hits in Layer3",8,0,8) 
    nRecoHitsL4FPIX = r.TH1F("nRecoHitsL4FPIX",";Number of FPIX Hits in Layer4",8,0,8) 

    nRecoHitsTotal= r.TH1F("nRecoHitsTotal",";Number of Total Layers with Hits",9,0,9) 

    nSimVsRecoHitsBPIX = r.TH2F("recovssimBPIX",";Number of BPIX Layers with Sim Hits;Number of BPIX Layers with Reco Hits",5,0.,5.,5,0.,5.) 
    nSimVsRecoHitsFPIX = r.TH2F("recovssimFPIX",";Number of FPIX Layers with Sim Hits;Number of FPIX Layers with Reco Hits",5,0.,5.,5,0.,5.) 
    nSimvsRecoDiffBPIX= r.TH1F("nSimvsRecoDiffBPIX",";Number of BPIX Layers with Sim Hits - Number of BPIX Layers with Reco Hits",12,-6,6) 
    nSimvsRecoDiffFPIX= r.TH1F("nSimvsRecoDiffFPIX",";Number of FPIX Layers with Sim Hits - Number of FPIX Layers with Reco Hits",12,-6,6) 


    recoElepT = r.TH1F("recoPt",";p_{T} [GeV]",20,0.,100.) 
    recoEleEta = r.TH1F("recoEta",";#eta",20,-3.,3.) 
    recoElePhi = r.TH1F("recoPhi",";#phi",21,-3.15,3.15) 

    recoZpeak = r.TH1F("recoMass",";mass [GeV]",20,50.,150.)


    f = r.TFile.Open(options.input)
    for event in f.Get("egammaReconstruction/tree"):
        for electron in range(0,len(event.nSimHitLayersBPIX)):
            #if(event.isMatched[electron]):
            nSimHits.Fill(event.nSimHitLayersBPIX[electron])  
            nSimHitsL1.Fill(event.nSimHitsLayer1_BPIX[electron])  
            nSimHitsL2.Fill(event.nSimHitsLayer2_BPIX[electron])  
            nSimHitsL3.Fill(event.nSimHitsLayer3_BPIX[electron])  
            nSimHitsL4.Fill(event.nSimHitsLayer4_BPIX[electron])  
            nSimHitsFPIX.Fill(event.nSimHitLayersFPIX[electron])  
            nSimHitsL1FPIX.Fill(event.nSimHitsLayer1_FPIX[electron])  
            nSimHitsL2FPIX.Fill(event.nSimHitsLayer2_FPIX[electron])  
            nSimHitsL3FPIX.Fill(event.nSimHitsLayer3_FPIX[electron])  
            nSimHitsL4FPIX.Fill(event.nSimHitsLayer4_FPIX[electron])  
            simElepT.Fill(event.simEle_pt[electron]) 
            simEleEta.Fill(event.simEle_eta[electron]) 
            simElePhi.Fill(event.simEle_phi[electron]) 

            totalLayers = 0
            if(event.nSimHitsLayer1_BPIX[electron]>0):
                totalLayers = totalLayers+1
            if(event.nSimHitsLayer2_BPIX[electron]>0):
                totalLayers = totalLayers+1
            if(event.nSimHitsLayer3_BPIX[electron]>0):
                totalLayers = totalLayers+1
            if(event.nSimHitsLayer4_BPIX[electron]>0):
                totalLayers = totalLayers+1
            if(event.nSimHitsLayer1_FPIX[electron]>0):
                totalLayers = totalLayers+1
            if(event.nSimHitsLayer2_FPIX[electron]>0):
                totalLayers = totalLayers+1
            if(event.nSimHitsLayer3_FPIX[electron]>0):
                totalLayers = totalLayers+1
            if(event.nSimHitsLayer4_FPIX[electron]>0):
                totalLayers = totalLayers+1
            nSimHitsTotal.Fill(totalLayers)

        if(len(event.simEle_pt)==2):
            v1 = r.TLorentzVector()
            v2 = r.TLorentzVector()
            v1.SetPtEtaPhiE(event.simEle_pt[0],event.simEle_eta[0],event.simEle_phi[0],event.simEle_E[0])
            v2.SetPtEtaPhiE(event.simEle_pt[1],event.simEle_eta[1],event.simEle_phi[1],event.simEle_E[1])
            simZpeak.Fill((v1+v2).M())

        for electron in range(0,len(event.nRecoHitLayersBPIX)):
            nRecoHits.Fill(event.nRecoHitLayersBPIX[electron])  
            nRecoHitsL1.Fill(event.nRecoHitsLayer1_BPIX[electron])  
            nRecoHitsL2.Fill(event.nRecoHitsLayer2_BPIX[electron])  
            nRecoHitsL3.Fill(event.nRecoHitsLayer3_BPIX[electron])  
            nRecoHitsL4.Fill(event.nRecoHitsLayer4_BPIX[electron])  
            nRecoHitsFPIX.Fill(event.nRecoHitLayersFPIX[electron])  
            nRecoHitsL1FPIX.Fill(event.nRecoHitsLayer1_FPIX[electron])  
            nRecoHitsL2FPIX.Fill(event.nRecoHitsLayer2_FPIX[electron])  
            nRecoHitsL3FPIX.Fill(event.nRecoHitsLayer3_FPIX[electron])  
            nRecoHitsL4FPIX.Fill(event.nRecoHitsLayer4_FPIX[electron])  
            index=0
            if(event.isMatched[electron]):
                recoElepT.Fill(event.recoEle_pt[index]) 
                recoEleEta.Fill(event.recoEle_eta[index]) 
                recoElePhi.Fill(event.recoEle_phi[index]) 
                index = index+1
            if(len(event.nSimHitLayersBPIX)==len(event.nRecoHitLayersBPIX)):
                nSimVsRecoHitsBPIX.Fill(event.nSimHitLayersBPIX[electron],event.nRecoHitLayersBPIX[electron])
                nSimVsRecoHitsFPIX.Fill(event.nSimHitLayersFPIX[electron],event.nRecoHitLayersFPIX[electron])
                nSimvsRecoDiffBPIX.Fill(event.nSimHitLayersBPIX[electron]-event.nRecoHitLayersBPIX[electron])
                nSimvsRecoDiffFPIX.Fill(event.nSimHitLayersFPIX[electron]-event.nRecoHitLayersFPIX[electron])

            totalLayers = 0
            if(event.nRecoHitsLayer1_BPIX[electron]>0):
                totalLayers = totalLayers+1
            if(event.nRecoHitsLayer2_BPIX[electron]>0):
                totalLayers = totalLayers+1
            if(event.nRecoHitsLayer3_BPIX[electron]>0):
                totalLayers = totalLayers+1
            if(event.nRecoHitsLayer4_BPIX[electron]>0):
                totalLayers = totalLayers+1
            if(event.nRecoHitsLayer1_FPIX[electron]>0):
                totalLayers = totalLayers+1
            if(event.nRecoHitsLayer2_FPIX[electron]>0):
                totalLayers = totalLayers+1
            if(event.nRecoHitsLayer3_FPIX[electron]>0):
                totalLayers = totalLayers+1
            if(event.nRecoHitsLayer4_FPIX[electron]>0):
                totalLayers = totalLayers+1
            nRecoHitsTotal.Fill(totalLayers)


        if(len(event.recoEle_pt)==2):
            v1 = r.TLorentzVector()
            v2 = r.TLorentzVector()
            v1.SetPtEtaPhiE(event.recoEle_pt[0],event.recoEle_eta[0],event.recoEle_phi[0],event.recoEle_E[0])
            v2.SetPtEtaPhiE(event.recoEle_pt[1],event.recoEle_eta[1],event.recoEle_phi[1],event.recoEle_E[1])
            recoZpeak.Fill((v1+v2).M())

    paveCMS = r.TPaveText(0.07,0.93,0.92,0.96,"NDC");
    paveCMS.AddText("#bf{CMS Simulation} (#it{Work In Progress)}                       Run-3 (14 TeV)")
    paveCMS.SetFillColor(0)
    paveCMS.SetBorderSize(0)
    paveCMS.SetTextSize(0.04)
    paveCMS.SetTextFont(42)

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimHits.SetLineColor(r.kRed)
    nSimHits.SetLineWidth(2)
    nRecoHits.SetLineColor(r.kBlue)
    nRecoHits.SetLineWidth(2)
    l0 = r.TLegend(.7, .72, .89, .85)
    l0.AddEntry(nSimHits,"Simulated","l")
    l0.AddEntry(nRecoHits,"Reconstructed","l")
    nSimHits.Draw("HIST")
    nRecoHits.Draw("HIST same")
    paveCMS.Draw("same")
    l0.Draw("same")
    c.SaveAs('nHitsBPIX.png')


    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimHitsTotal.SetLineColor(r.kRed)
    nSimHitsTotal.SetLineWidth(2)
    nRecoHitsTotal.SetLineColor(r.kBlue)
    nRecoHitsTotal.SetLineWidth(2)
    nSimHitsTotal.Draw("HIST")
    nRecoHitsTotal.Draw("HIST same")
    paveCMS.Draw("same")
    l0.Draw("same")
    c.SaveAs('nHitsTotal.png')


    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimHitsL1.SetLineColor(r.kRed)
    nSimHitsL1.SetLineWidth(2)
    nRecoHitsL1.SetLineColor(r.kBlue)
    nRecoHitsL1.SetLineWidth(2)
    nRecoHitsL1.Draw("HIST")
    nSimHitsL1.Draw("HIST same")
    paveCMS.Draw("same")
    l0.Draw("same")
    c.SaveAs('nHitsL1BPIX.png')

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimHitsL2.SetLineColor(r.kRed)
    nSimHitsL2.SetLineWidth(2)
    nRecoHitsL2.SetLineColor(r.kBlue)
    nRecoHitsL2.SetLineWidth(2)
    nRecoHitsL2.Draw("HIST")
    nSimHitsL2.Draw("HIST same")
    paveCMS.Draw("same") 
    l0.Draw("same")
    c.SaveAs('nHitsL2BPIX.png')

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimHitsL3.SetLineColor(r.kRed)
    nSimHitsL3.SetLineWidth(2)
    nRecoHitsL3.SetLineColor(r.kBlue)
    nRecoHitsL3.SetLineWidth(2)
    nRecoHitsL3.Draw("HIST")
    nSimHitsL3.Draw("HIST same")
    paveCMS.Draw("same")
    l0.Draw("same") 
    c.SaveAs('nHitsL3BPIX.png')    

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimHitsL4.SetLineColor(r.kRed)
    nSimHitsL4.SetLineWidth(2)
    nRecoHitsL4.SetLineColor(r.kBlue)
    nRecoHitsL4.SetLineWidth(2)
    nRecoHitsL4.Draw("HIST")
    nSimHitsL4.Draw("HIST same")
    paveCMS.Draw("same")
    l0.Draw("same")
    c.SaveAs('nHitsL4BPIX.png')  

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimHitsFPIX.SetLineColor(r.kRed)
    nSimHitsFPIX.SetLineWidth(2)
    nRecoHitsFPIX.SetLineColor(r.kBlue)
    nRecoHitsFPIX.SetLineWidth(2)
    l0 = r.TLegend(.7, .72, .89, .85)
    l0.AddEntry(nSimHitsFPIX,"Simulated","l")
    l0.AddEntry(nRecoHitsFPIX,"Reconstructed","l")
    nSimHitsFPIX.Draw("HIST")
    nRecoHitsFPIX.Draw("HIST same")
    paveCMS.Draw("same")
    l0.Draw("same")
    c.SaveAs('nHitsFPIX.png')

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimHitsL1FPIX.SetLineColor(r.kRed)
    nSimHitsL1FPIX.SetLineWidth(2)
    nRecoHitsL1FPIX.SetLineColor(r.kBlue)
    nRecoHitsL1FPIX.SetLineWidth(2)
    nRecoHitsL1FPIX.Draw("HIST")
    nSimHitsL1FPIX.Draw("HIST same")
    paveCMS.Draw("same")
    l0.Draw("same")
    c.SaveAs('nHitsL1FPIX.png')

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimHitsL2FPIX.SetLineColor(r.kRed)
    nSimHitsL2FPIX.SetLineWidth(2)
    nRecoHitsL2FPIX.SetLineColor(r.kBlue)
    nRecoHitsL2FPIX.SetLineWidth(2)
    nRecoHitsL2FPIX.Draw("HIST")
    nSimHitsL2FPIX.Draw("HIST same")
    paveCMS.Draw("same") 
    l0.Draw("same")
    c.SaveAs('nHitsL2FPIX.png')

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimHitsL3FPIX.SetLineColor(r.kRed)
    nSimHitsL3FPIX.SetLineWidth(2)
    nRecoHitsL3FPIX.SetLineColor(r.kBlue)
    nRecoHitsL3FPIX.SetLineWidth(2)
    nRecoHitsL3FPIX.Draw("HIST")
    nSimHitsL3FPIX.Draw("HIST same")
    paveCMS.Draw("same")
    l0.Draw("same") 
    c.SaveAs('nHitsL3FPIX.png')    

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimHitsL4FPIX.SetLineColor(r.kRed)
    nSimHitsL4FPIX.SetLineWidth(2)
    nRecoHitsL4FPIX.SetLineColor(r.kBlue)
    nRecoHitsL4FPIX.SetLineWidth(2)
    nRecoHitsL4FPIX.Draw("HIST")
    nSimHitsL4FPIX.Draw("HIST same")
    paveCMS.Draw("same")
    l0.Draw("same")
    c.SaveAs('nHitsL4FPIX.png')  

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    simElepT.SetLineColor(r.kRed)
    simElepT.SetLineWidth(2)
    recoElepT.SetLineColor(r.kBlue)
    recoElepT.SetLineWidth(2)
    simElepT.Draw("ep")
    recoElepT.Draw("ep same")
    paveCMS.Draw("same")
    l0.Draw("same")
    c.SaveAs('elePt.png')  

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    simEleEta.SetLineColor(r.kRed)
    simEleEta.SetLineWidth(2)
    recoEleEta.SetLineColor(r.kBlue)
    recoEleEta.SetLineWidth(2)
    simEleEta.Draw("ep")
    recoEleEta.Draw("ep same")
    paveCMS.Draw("same")
    l0.Draw("same")
    c.SaveAs('eleEta.png')  

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    simElePhi.SetLineColor(r.kRed)
    simElePhi.SetLineWidth(2)
    recoElePhi.SetLineColor(r.kBlue)
    recoElePhi.SetLineWidth(2)
    simElePhi.SetMaximum(150)
    simElePhi.SetMinimum(70)
    simElePhi.Draw("ep")
    recoElePhi.Draw("ep same")
    paveCMS.Draw("same")
    l0.Draw("same")
    c.SaveAs('elePhi.png')  

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    simZpeak.SetLineColor(r.kRed)
    simZpeak.SetLineWidth(2)
    recoZpeak.SetLineColor(r.kBlue)
    recoZpeak.SetLineWidth(2)
    simZpeak.Draw("ep")
    recoZpeak.Draw("ep same")
    paveCMS.Draw("same")
    l0.Draw("same")
    c.SaveAs('zpeak.png')  

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimVsRecoHitsBPIX.Draw("COLZ")
    paveCMS.Draw("same")
    c.SaveAs('SimVsRecoBPIX.png')  

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimVsRecoHitsFPIX.Draw("COLZ")
    paveCMS.Draw("same")
    c.SaveAs('SimVsRecoFPIX.png')  

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimvsRecoDiffBPIX.Draw("HIST")
    paveCMS.Draw("same")
    c.SaveAs('SimDiffRecoBPIX.png')  

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimvsRecoDiffFPIX.Draw("HIST")
    paveCMS.Draw("same")
    c.SaveAs('SimDiffRecoFPIX.png')  


if __name__ == '__main__':
    options, paths = parse_arguments()
    main(options = options, paths = paths)