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

def plotHitsPerModule(hists,name):

    paveCMS = r.TPaveText(0.07,0.93,0.92,0.96,"NDC");
    paveCMS.AddText("#bf{CMS Simulation} (#it{Work In Progress)}                       Run-3 (14 TeV)")
    paveCMS.SetFillColor(0)
    paveCMS.SetBorderSize(0)
    paveCMS.SetTextSize(0.04)
    paveCMS.SetTextFont(42)

    colors = [1,2,4,6,7,8,9,46,30]
    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()

    for m in range(0,8):
        hists[m].SetLineColor(colors[m])
        hists[m].SetLineWidth(2)

    l1 = r.TLegend(.7, .55, .89, .85)
    l1.AddEntry(hists[0],"Module 1","l")   
    l1.AddEntry(hists[1],"Module 2","l")   
    l1.AddEntry(hists[2],"Module 3","l")   
    l1.AddEntry(hists[3],"Module 4","l")   
    l1.AddEntry(hists[4],"Module 5","l")   
    l1.AddEntry(hists[5],"Module 6","l")   
    l1.AddEntry(hists[6],"Module 7","l")   
    l1.AddEntry(hists[7],"Module 8","l")   

    hists[7].Draw("HIST")

    for m in range(0,7):        
        hists[m].Draw("HIST same")
        
    paveCMS.Draw("same")
    l1.Draw("same")

    c.SaveAs(name)

def doubletComb(BPX1,BPX2,BPX3,BPX4,FPXp1,FPXp2,FPXp3,FPXn1,FPXn2,FPXn3):

    doublets = [0]*20

    if(BPX1>0 and BPX2>0):
        doublets[0]=1
    if(BPX1>0 and FPXp1>0):
        doublets[1]=1
    if(BPX1>0 and FPXn1>0):
        doublets[2]=1
    if(BPX2>0 and BPX3>0):
        doublets[3]=1
    if(BPX2>0 and FPXp1>0):
        doublets[4]=1
    if(BPX2>0 and FPXn1>0):
        doublets[5]=1
    if(FPXp1>0 and FPXp2>0):
        doublets[6]=1
    if(FPXn1>0 and FPXn2>0):
        doublets[7]=1
    if(BPX3>0 and BPX4>0):
        doublets[8]=1
    if(BPX3>0 and FPXp1>0):
        doublets[9]=1
    if(BPX3>0 and FPXn1>0):
        doublets[10]=1
    if(FPXp2>0 and FPXp3>0):
        doublets[11]=1
    if(FPXn2>0 and FPXn3>0):
        doublets[12]=1
    if(BPX1>0 and BPX3>0):
        doublets[13]=1
    if(BPX2>0 and BPX4>0):
        doublets[14]=1
    if(BPX1>0 and FPXp2>0):
        doublets[15]=1
    if(BPX1>0 and FPXn2>0):
        doublets[16]=1
    if(FPXp1>0 and FPXp3>0):
        doublets[17]=1
    if(FPXn1>0 and FPXn3>0):
        doublets[18]=1

    # Check if there are any other type of combinations
    if(sum(doublets)==0 and ((min(BPX1,1)+min(BPX2,1)+min(BPX3,1)+min(BPX4,1)+min(FPXp1,1)+min(FPXp2,1)+min(FPXp3,1)+min(FPXn1,1)+min(FPXn2,1)+min(FPXn3,1))>1)):  
        doublets[19]=1

    return doublets

# Triplet types 
def buildTriplets(pixelDoubletTypes):
    numberOfTriplets = 0
    tripletList = [[3,4,5,14],[6,7],[7,18],[8,9,10],[6,17],[7,18],[11],[12],[],[6,17],[7,18],[],[],[8,9,10],[],[11],[12],[],[],[]]
    for doublet1 in range(0,len(pixelDoubletTypes)):
        for doublet2 in range(doublet1,len(pixelDoubletTypes)):
            if pixelDoubletTypes[doublet2] in tripletList[pixelDoubletTypes[doublet1]]:
                #print('Lala',pixelDoubletTypes[doublet1],pixelDoubletTypes[doublet2])
                #print(tripletList[pixelDoubletTypes[doublet1]])
                numberOfTriplets = numberOfTriplets + 1

    return numberOfTriplets


def doubletQual(hits_x,hits_y,hits_z,hits_l):

    minz = [-20., 0., -30., -22., 10., -30., -70., -70., -22., 15., -30, -70., -70., -20., -22., 0, -30., -70., -70.,0]
    maxz = [20., 30., 0., 22., 30., -10., 70., 70., 22., 30., -15., 70., 70., 20., 22., 30., 0., 70., 70.,0]
    maxr = [20., 9., 9., 20., 7., 7., 5., 5., 20., 6., 6., 5., 5., 20., 20., 9., 9., 9., 9.,0]

    comb = []
    collection = [0]*20

    for ihit in range(0,len(hits_l)-1):
        for jhit in range(ihit+1,len(hits_l)):
            hitIdLayer = [0]*10
            # Skip checking hits in the same layer
            if(hits_l[ihit]==hits_l[jhit]):
                continue
            hitIdLayer[hits_l[ihit]-1] = 1
            hitIdLayer[hits_l[jhit]-1] = 1
            doublet = doubletComb(hitIdLayer[0],hitIdLayer[1],hitIdLayer[2],hitIdLayer[3],hitIdLayer[4],hitIdLayer[6],hitIdLayer[8],hitIdLayer[5],hitIdLayer[7],hitIdLayer[9])
            combination = doublet.index(max(doublet))
            # Duplicate combinations
            if(collection[combination]==1):
                continue
            else:
                collection[combination]=1

            # Quality cuts
            idIn = jhit
            idOut = ihit
            if(hits_l[ihit]<hits_l[jhit]):
                idIn = ihit
                idOut = jhit
            z = hits_z[idIn]
            if((minz[combination]>z) or (maxz[combination]<z)):
                continue
            dr = m.sqrt(hits_z[idOut]*hits_z[idOut] + hits_x[idOut]*hits_x[idOut] + hits_y[idOut]*hits_y[idOut]) - m.sqrt(hits_z[idIn]*hits_z[idIn] + hits_x[idIn]*hits_x[idIn] + hits_y[idIn]*hits_y[idIn])
            if(dr> maxr[combination]):
                continue
            #print("combination, z , Dr : ",combination, z , dr)
            comb.append(combination)
    
    # Remove "Other combinations" encoded with 19 for tracking particles that have at least one doublet in the patatrack list
    hasCombination = False
    hasAtLeastOne = False
    for combination in range(0,len(comb)):
        if(comb[combination] < 19 and comb[combination] > 0):
            hasAtLeastOne = True
        if(comb[combination] == 19):
            hasCombination = True
    if(hasCombination and hasAtLeastOne):
        comb.remove(19)

    return comb


def main(options, paths):

    noGsfEleMatching = False
    tree = "egammaReconstructionEE"

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

    nRecoHitsTotal= r.TH1F("nRecoHitsTotal",";Number of Total Layers with Hits",9,0,9) 

    nSimVsRecoHitsBPIX = r.TH2F("recovssimBPIX",";Number of BPIX Layers with Sim Hits;Number of BPIX Layers with Reco Hits",5,0.,5.,5,0.,5.) 
    nSimVsRecoHitsFPIX = r.TH2F("recovssimFPIX",";Number of FPIX Layers with Sim Hits;Number of FPIX Layers with Reco Hits",5,0.,5.,5,0.,5.) 
    nSimvsRecoDiffBPIX= r.TH1F("nSimvsRecoDiffBPIX",";Number of BPIX Layers with Sim Hits - Number of BPIX Layers with Reco Hits",12,-6,6) 
    nSimvsRecoDiffFPIX= r.TH1F("nSimvsRecoDiffFPIX",";Number of FPIX Layers with Sim Hits - Number of FPIX Layers with Reco Hits",12,-6,6) 

    nSimVsRecoHitsTotal = r.TH2F("recovssimTotal",";Total number of Layers with Sim Hits;Total number of Layers with Reco Hits",8,0.,8.,8,0.,8.) 
    nSimvsRecoDiffTotal= r.TH1F("nSimvsRecoDiffTotal",";Total number of Layers with Sim Hits - Total number of Layers with Reco Hits",16,-8,8) 

    numberOfPixelDoubletsSim = r.TH1F("numberOfPixelDoubletsSim",";Number of Pixel doublets per TP",10,0,10) 
    numberOfPixelDoubletsReco = r.TH1F("numberOfPixelDoubletsReco",";Number of Pixel doublets per TP",10,0,10) 

    pixelDoubletTypeSim = r.TH1F("pixelDoubletTypeSim",";Type of Pixel doublet ",20,0,20) 
    pixelDoubletTypeReco  = r.TH1F("pixelDoubletTypeReco",";Type of Pixel doublet ",20,0,20) 

    test = r.TH1F("test",";Type of Pixel doublet ",20,0,20) 


    recoElepT = r.TH1F("recoPt",";p_{T} [GeV]",20,0.,100.) 
    recoEleEta = r.TH1F("recoEta",";#eta",20,-3.,3.) 
    recoElePhi = r.TH1F("recoPhi",";#phi",21,-3.15,3.15) 

    recoZpeak = r.TH1F("recoMass",";mass [GeV]",20,50.,150.)

    sim_moduleHistsBPIX_L1=[]
    sim_moduleHistsBPIX_L2=[]
    sim_moduleHistsBPIX_L3=[]
    sim_moduleHistsBPIX_L4=[]
    sim_moduleHistsFPIX_L1=[]
    sim_moduleHistsFPIX_L2=[]
    sim_moduleHistsFPIX_L3=[]

    for module in range(0,8):  
        sim_moduleHistsBPIX_L1.append(r.TH1F("nSimHitsPerModuleL1BPIX_"+str(module),";Number of hits per module in BPIX Layer1 ",8,0,8))
        sim_moduleHistsBPIX_L2.append(r.TH1F("nSimHitsPerModuleL2BPIX_"+str(module),";Number of hits per module in BPIX Layer2 ",8,0,8))
        sim_moduleHistsBPIX_L3.append(r.TH1F("nSimHitsPerModuleL3BPIX_"+str(module),";Number of hits per module in BPIX Layer3 ",8,0,8))
        sim_moduleHistsBPIX_L4.append(r.TH1F("nSimHitsPerModuleL4BPIX_"+str(module),";Number of hits per module in BPIX Layer4 ",8,0,8))
        sim_moduleHistsFPIX_L1.append(r.TH1F("nSimHitsPerModuleL1FPIX_"+str(module),";Number of hits per module in FPIX Layer1 ",8,0,8))
        sim_moduleHistsFPIX_L2.append(r.TH1F("nSimHitsPerModuleL2FPIX_"+str(module),";Number of hits per module in FPIX Layer2 ",8,0,8))
        sim_moduleHistsFPIX_L3.append(r.TH1F("nSimHitsPerModuleL3FPIX_"+str(module),";Number of hits per module in FPIX Layer3 ",8,0,8))

    reco_moduleHistsBPIX_L1=[]
    reco_moduleHistsBPIX_L2=[]
    reco_moduleHistsBPIX_L3=[]
    reco_moduleHistsBPIX_L4=[]
    reco_moduleHistsFPIX_L1=[]
    reco_moduleHistsFPIX_L2=[]
    reco_moduleHistsFPIX_L3=[]

    for module in range(0,8):  
        reco_moduleHistsBPIX_L1.append(r.TH1F("nRecoHitsPerModuleL1BPIX_"+str(module),";Number of hits per module in BPIX Layer1 ",8,0,8))
        reco_moduleHistsBPIX_L2.append(r.TH1F("nRecoHitsPerModuleL2BPIX_"+str(module),";Number of hits per module in BPIX Layer2 ",8,0,8))
        reco_moduleHistsBPIX_L3.append(r.TH1F("nRecoHitsPerModuleL3BPIX_"+str(module),";Number of hits per module in BPIX Layer3 ",8,0,8))
        reco_moduleHistsBPIX_L4.append(r.TH1F("nRecoHitsPerModuleL4BPIX_"+str(module),";Number of hits per module in BPIX Layer4 ",8,0,8))
        reco_moduleHistsFPIX_L1.append(r.TH1F("nRecoHitsPerModuleL1FPIX_"+str(module),";Number of hits per module in FPIX Layer1 ",8,0,8))
        reco_moduleHistsFPIX_L2.append(r.TH1F("nRecoHitsPerModuleL2FPIX_"+str(module),";Number of hits per module in FPIX Layer2 ",8,0,8))
        reco_moduleHistsFPIX_L3.append(r.TH1F("nRecoHitsPerModuleL3FPIX_"+str(module),";Number of hits per module in FPIX Layer3 ",8,0,8))


    f = r.TFile.Open(options.input)
    for event in f.Get(tree+"/tree"):
        for electron in range(0,len(event.nSimHitLayersBPIX)):
            if(event.isMatched[electron] or noGsfEleMatching):

                nSimHits.Fill(event.nSimHitLayersBPIX[electron])  
                nSimHitsL1.Fill(event.nSimHitsLayer1_BPIX[electron])  
                nSimHitsL2.Fill(event.nSimHitsLayer2_BPIX[electron])  
                nSimHitsL3.Fill(event.nSimHitsLayer3_BPIX[electron])  
                nSimHitsL4.Fill(event.nSimHitsLayer4_BPIX[electron])  
                nSimHitsFPIX.Fill(event.nSimHitLayersFPIX[electron])  
                nSimHitsL1FPIX.Fill(event.nSimHitsLayer1_FPIX_pos[electron] + event.nSimHitsLayer1_FPIX_neg[electron])  
                nSimHitsL2FPIX.Fill(event.nSimHitsLayer2_FPIX_pos[electron] + event.nSimHitsLayer2_FPIX_neg[electron])  
                nSimHitsL3FPIX.Fill(event.nSimHitsLayer3_FPIX_pos[electron] + event.nSimHitsLayer3_FPIX_neg[electron])  

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
                if(event.nSimHitsLayer1_FPIX_pos[electron]>0):
                    totalLayers = totalLayers+1
                if(event.nSimHitsLayer1_FPIX_neg[electron]>0):
                    totalLayers = totalLayers+1    
                if(event.nSimHitsLayer2_FPIX_pos[electron]>0):
                    totalLayers = totalLayers+1    
                if(event.nSimHitsLayer2_FPIX_neg[electron]>0):
                    totalLayers = totalLayers+1
                if(event.nSimHitsLayer3_FPIX_pos[electron]>0):
                    totalLayers = totalLayers+1
                if(event.nSimHitsLayer3_FPIX_neg[electron]>0):
                    totalLayers = totalLayers+1
                nSimHitsTotal.Fill(totalLayers)

                doubletsSim = doubletComb(event.nSimHitsLayer1_BPIX[electron],event.nSimHitsLayer2_BPIX[electron],event.nSimHitsLayer3_BPIX[electron],event.nSimHitsLayer4_BPIX[electron],
                                    event.nSimHitsLayer1_FPIX_pos[electron],event.nSimHitsLayer2_FPIX_pos[electron],event.nSimHitsLayer3_FPIX_pos[electron],
                                    event.nSimHitsLayer1_FPIX_neg[electron],event.nSimHitsLayer2_FPIX_neg[electron],event.nSimHitsLayer3_FPIX_neg[electron])
            
                for index in range(0,len(doubletsSim)):
                    pixelDoubletTypeSim.Fill(index,doubletsSim[index])
    
                numberOfPixelDoubletsSim.Fill(sum(doubletsSim))

        ## Number of hits per module for each Pixel layer 
        for part in range(0,len(event.nSimHitsLayer1_BPIX)):
            if(event.isMatched[part] or noGsfEleMatching):
                mod = [0] * 8
                for hit in range(0,len(event.nSimHitsPerModuleLayer1_BPIX[part])):
                    for m in range(1,9):
                        if(event.nSimHitsPerModuleLayer1_BPIX[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for m in range(0,8):
                    sim_moduleHistsBPIX_L1[m].Fill(mod[m])
        ##
                mod = [0] * 8
                for hit in range(0,len(event.nSimHitsPerModuleLayer2_BPIX[part])):
                    for m in range(1,9):
                        if(event.nSimHitsPerModuleLayer2_BPIX[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for m in range(0,8):
                    sim_moduleHistsBPIX_L2[m].Fill(mod[m])
        ##
                mod = [0] * 8
                for hit in range(0,len(event.nSimHitsPerModuleLayer3_BPIX[part])):
                    for m in range(1,9):
                        if(event.nSimHitsPerModuleLayer3_BPIX[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for m in range(0,8):
                    sim_moduleHistsBPIX_L3[m].Fill(mod[m])
        ##
                mod = [0] * 8
                for hit in range(0,len(event.nSimHitsPerModuleLayer4_BPIX[part])):
                    for m in range(1,9):
                        if(event.nSimHitsPerModuleLayer4_BPIX[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for m in range(0,8):
                    sim_moduleHistsBPIX_L4[m].Fill(mod[m])                
        ##
                mod = [0] * 8
                for hit in range(0,len(event.nSimHitsPerModuleLayer1_FPIX_pos[part])):
                    for m in range(1,9):
                        if(event.nSimHitsPerModuleLayer1_FPIX_pos[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for hit in range(0,len(event.nSimHitsPerModuleLayer1_FPIX_neg[part])):
                    for m in range(1,9):
                        if(event.nSimHitsPerModuleLayer1_FPIX_neg[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for m in range(0,8):
                    sim_moduleHistsFPIX_L1[m].Fill(mod[m])
        ##
                mod = [0] * 8
                for hit in range(0,len(event.nSimHitsPerModuleLayer2_FPIX_pos[part])):
                    for m in range(1,9):
                        if(event.nSimHitsPerModuleLayer2_FPIX_pos[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for hit in range(0,len(event.nSimHitsPerModuleLayer2_FPIX_neg[part])):
                    for m in range(1,9):
                        if(event.nSimHitsPerModuleLayer2_FPIX_neg[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for m in range(0,8):
                    sim_moduleHistsFPIX_L2[m].Fill(mod[m])
        ##
                mod = [0] * 8
                for hit in range(0,len(event.nSimHitsPerModuleLayer3_FPIX_pos[part])):
                    for m in range(1,9):
                        if(event.nSimHitsPerModuleLayer3_FPIX_pos[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for hit in range(0,len(event.nSimHitsPerModuleLayer3_FPIX_neg[part])):
                    for m in range(1,9):
                        if(event.nSimHitsPerModuleLayer3_FPIX_neg[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for m in range(0,8):
                    sim_moduleHistsFPIX_L3[m].Fill(mod[m])              

        if(len(event.simEle_pt)==2):
            v1 = r.TLorentzVector()
            v2 = r.TLorentzVector()
            v1.SetPtEtaPhiE(event.simEle_pt[0],event.simEle_eta[0],event.simEle_phi[0],event.simEle_E[0])
            v2.SetPtEtaPhiE(event.simEle_pt[1],event.simEle_eta[1],event.simEle_phi[1],event.simEle_E[1])
            simZpeak.Fill((v1+v2).M())

        for electron in range(0,len(event.nRecoHitLayersBPIX)):
            if(event.isMatched[electron] or noGsfEleMatching):

                # Build doublets and triplets from hits and apply quality cuts 
                doubletsQuality = doubletQual(event.recHit_x[electron],event.recHit_y[electron],event.recHit_z[electron],event.recHit_l[electron])
                for i in range(0,len(doubletsQuality)):
                    test.Fill(doubletsQuality[i],1)

                # Build triplets
                triplets = buildTriplets(doubletsQuality)
                print('Number of Triplets : ',triplets)

                nRecoHits.Fill(event.nRecoHitLayersBPIX[electron])  
                nRecoHitsL1.Fill(event.nRecoHitsLayer1_BPIX[electron])  
                nRecoHitsL2.Fill(event.nRecoHitsLayer2_BPIX[electron])  
                nRecoHitsL3.Fill(event.nRecoHitsLayer3_BPIX[electron])  
                nRecoHitsL4.Fill(event.nRecoHitsLayer4_BPIX[electron])  
                nRecoHitsFPIX.Fill(event.nRecoHitLayersFPIX[electron])  
                nRecoHitsL1FPIX.Fill(event.nRecoHitsLayer1_FPIX_pos[electron] + event.nRecoHitsLayer1_FPIX_neg[electron])  
                nRecoHitsL2FPIX.Fill(event.nRecoHitsLayer2_FPIX_pos[electron] + event.nRecoHitsLayer2_FPIX_neg[electron])  
                nRecoHitsL3FPIX.Fill(event.nRecoHitsLayer3_FPIX_pos[electron] + event.nRecoHitsLayer3_FPIX_neg[electron])  

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

                totalLayersReco = 0
                if(event.nRecoHitsLayer1_BPIX[electron]>0):
                    totalLayersReco = totalLayersReco+1
                if(event.nRecoHitsLayer2_BPIX[electron]>0):
                    totalLayersReco = totalLayersReco+1
                if(event.nRecoHitsLayer3_BPIX[electron]>0):
                    totalLayersReco = totalLayersReco+1
                if(event.nRecoHitsLayer4_BPIX[electron]>0):
                    totalLayersReco = totalLayersReco+1
                if(event.nRecoHitsLayer1_FPIX_pos[electron]>0):
                    totalLayersReco = totalLayersReco+1
                if(event.nRecoHitsLayer1_FPIX_neg[electron]>0):
                    totalLayersReco = totalLayersReco+1                
                if(event.nRecoHitsLayer2_FPIX_pos[electron]>0):
                    totalLayersReco = totalLayersReco+1
                if(event.nRecoHitsLayer2_FPIX_neg[electron]>0):
                    totalLayersReco = totalLayersReco+1                
                if(event.nRecoHitsLayer3_FPIX_pos[electron]>0):
                    totalLayersReco = totalLayersReco+1
                if(event.nRecoHitsLayer3_FPIX_neg[electron]>0):
                    totalLayersReco = totalLayersReco+1

                nRecoHitsTotal.Fill(totalLayersReco)

                doubletsReco = doubletComb(event.nRecoHitsLayer1_BPIX[electron],event.nRecoHitsLayer2_BPIX[electron],event.nRecoHitsLayer3_BPIX[electron],event.nRecoHitsLayer4_BPIX[electron],
                                    event.nRecoHitsLayer1_FPIX_pos[electron],event.nRecoHitsLayer2_FPIX_pos[electron],event.nRecoHitsLayer3_FPIX_pos[electron],
                                    event.nRecoHitsLayer1_FPIX_neg[electron],event.nRecoHitsLayer2_FPIX_neg[electron],event.nRecoHitsLayer3_FPIX_neg[electron])

                for index in range(0,len(doubletsReco)):
                    pixelDoubletTypeReco.Fill(index,doubletsReco[index])

                numberOfPixelDoubletsReco.Fill(sum(doubletsReco))

                if(len(event.nSimHitLayersBPIX)==len(event.nRecoHitLayersBPIX)):
                    nSimVsRecoHitsTotal.Fill(totalLayers,totalLayersReco)
                    nSimvsRecoDiffTotal.Fill(totalLayers-totalLayersReco)

        ## Number of hits per module for each Pixel layer 
        for part in range(0,len(event.nRecoHitsLayer1_BPIX)):
            if(event.isMatched[part] or noGsfEleMatching):
                mod = [0] * 8
                for hit in range(0,len(event.nRecoHitsPerModuleLayer1_BPIX[part])):
                    for m in range(1,9):
                        if(event.nRecoHitsPerModuleLayer1_BPIX[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for m in range(0,8):
                    reco_moduleHistsBPIX_L1[m].Fill(mod[m])
        ##
                mod = [0] * 8
                for hit in range(0,len(event.nRecoHitsPerModuleLayer2_BPIX[part])):
                    for m in range(1,9):
                        if(event.nRecoHitsPerModuleLayer2_BPIX[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for m in range(0,8):
                    reco_moduleHistsBPIX_L2[m].Fill(mod[m])
        ##
                mod = [0] * 8
                for hit in range(0,len(event.nRecoHitsPerModuleLayer3_BPIX[part])):
                    for m in range(1,9):
                        if(event.nRecoHitsPerModuleLayer3_BPIX[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for m in range(0,8):
                    reco_moduleHistsBPIX_L3[m].Fill(mod[m])
        ##
                mod = [0] * 8
                for hit in range(0,len(event.nRecoHitsPerModuleLayer4_BPIX[part])):
                    for m in range(1,9):
                        if(event.nRecoHitsPerModuleLayer4_BPIX[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for m in range(0,8):
                    reco_moduleHistsBPIX_L4[m].Fill(mod[m])                
        ##
                mod = [0] * 8
                for hit in range(0,len(event.nRecoHitsPerModuleLayer1_FPIX_pos[part])):
                    for m in range(1,9):
                        if(event.nRecoHitsPerModuleLayer1_FPIX_pos[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for hit in range(0,len(event.nRecoHitsPerModuleLayer1_FPIX_neg[part])):
                    for m in range(1,9):
                        if(event.nRecoHitsPerModuleLayer1_FPIX_neg[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for m in range(0,8):
                    reco_moduleHistsFPIX_L1[m].Fill(mod[m])
        ##
                mod = [0] * 8
                for hit in range(0,len(event.nRecoHitsPerModuleLayer2_FPIX_pos[part])):
                    for m in range(1,9):
                        if(event.nRecoHitsPerModuleLayer2_FPIX_pos[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for hit in range(0,len(event.nRecoHitsPerModuleLayer2_FPIX_neg[part])):
                    for m in range(1,9):
                        if(event.nRecoHitsPerModuleLayer2_FPIX_neg[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for m in range(0,8):
                    reco_moduleHistsFPIX_L2[m].Fill(mod[m])
        ##
                mod = [0] * 8
                for hit in range(0,len(event.nRecoHitsPerModuleLayer3_FPIX_pos[part])):
                    for m in range(1,9):
                        if(event.nRecoHitsPerModuleLayer3_FPIX_pos[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for hit in range(0,len(event.nRecoHitsPerModuleLayer3_FPIX_neg[part])):
                    for m in range(1,9):
                        if(event.nRecoHitsPerModuleLayer3_FPIX_neg[part][hit]==m):
                            mod[m-1] = mod[m-1] + 1
                for m in range(0,8):
                    reco_moduleHistsFPIX_L3[m].Fill(mod[m])   

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

    dir = 'plots'

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
    c.SaveAs(dir+'/nHitsBPIX.png')

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
    c.SaveAs(dir+'/nHitsTotal.png')

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
    c.SaveAs(dir+'/nHitsL1BPIX.png')

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
    c.SaveAs(dir+'/nHitsL2BPIX.png')

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
    c.SaveAs(dir+'/nHitsL3BPIX.png')    

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
    c.SaveAs(dir+'/nHitsL4BPIX.png')  

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
    c.SaveAs(dir+'/nHitsFPIX.png')

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
    c.SaveAs(dir+'/nHitsL1FPIX.png')

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
    c.SaveAs(dir+'/nHitsL2FPIX.png')

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
    c.SaveAs(dir+'/nHitsL3FPIX.png')    

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
    c.SaveAs(dir+'/elePt.png')  

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
    c.SaveAs(dir+'/eleEta.png')  

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
    c.SaveAs(dir+'/elePhi.png')  

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
    c.SaveAs(dir+'/zpeak.png')  

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimVsRecoHitsBPIX.Draw("COLZ")
    paveCMS.Draw("same")
    c.SaveAs(dir+'/SimVsRecoBPIX.png')  

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimVsRecoHitsFPIX.Draw("COLZ")
    paveCMS.Draw("same")
    c.SaveAs(dir+'/SimVsRecoFPIX.png')  

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimVsRecoHitsTotal.Draw("COLZ TEXT")
    paveCMS.Draw("same")
    c.SaveAs(dir+'/SimVsRecoTotal.png')  

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimvsRecoDiffBPIX.Draw("HIST")
    paveCMS.Draw("same")
    c.SaveAs(dir+'/SimDiffRecoBPIX.png')  

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimvsRecoDiffFPIX.Draw("HIST")
    paveCMS.Draw("same")
    c.SaveAs(dir+'/SimDiffRecoFPIX.png')  

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    nSimvsRecoDiffTotal.Draw("HIST")
    paveCMS.Draw("same")
    c.SaveAs(dir+'/SimDiffRecoTotal.png')  

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    pixelDoubletTypeSim.SetLineColor(r.kRed)
    pixelDoubletTypeSim.SetLineWidth(2)
    pixelDoubletTypeReco.SetLineColor(r.kBlue)
    pixelDoubletTypeReco.SetLineWidth(2)
    pixelDoubletTypeSim.Draw("HIST")
    pixelDoubletTypeReco.Draw("HIST same")
    paveCMS.Draw("same")
    l0.Draw("same")
    c.SaveAs(dir+'/pixelDoubletTypes.png')  

    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    numberOfPixelDoubletsSim.SetLineColor(r.kRed)
    numberOfPixelDoubletsSim.SetLineWidth(2)
    numberOfPixelDoubletsReco.SetLineColor(r.kBlue)
    numberOfPixelDoubletsReco.SetLineWidth(2)
    numberOfPixelDoubletsSim.Draw("HIST")
    numberOfPixelDoubletsReco.Draw("HIST same")
    paveCMS.Draw("same")
    l0.Draw("same")
    c.SaveAs(dir+'/NumberOfDoublets.png')  

    c = plotHitsPerModule(sim_moduleHistsBPIX_L1,dir+"/nSimHitsPerModule_L1_BPIX.png")
    c = plotHitsPerModule(sim_moduleHistsBPIX_L2,dir+"/nSimHitsPerModule_L2_BPIX.png")
    c = plotHitsPerModule(sim_moduleHistsBPIX_L3,dir+"/nSimHitsPerModule_L3_BPIX.png")
    c = plotHitsPerModule(sim_moduleHistsBPIX_L4,dir+"/nSimHitsPerModule_L4_BPIX.png")
    c = plotHitsPerModule(sim_moduleHistsFPIX_L1,dir+"/nSimHitsPerModule_L1_FPIX.png")
    c = plotHitsPerModule(sim_moduleHistsFPIX_L2,dir+"/nSimHitsPerModule_L2_FPIX.png")
    c = plotHitsPerModule(sim_moduleHistsFPIX_L3,dir+"/nSimHitsPerModule_L3_FPIX.png")

    c = plotHitsPerModule(reco_moduleHistsBPIX_L1,dir+"/nRecoHitsPerModule_L1_BPIX.png")
    c = plotHitsPerModule(reco_moduleHistsBPIX_L2,dir+"/nRecoHitsPerModule_L2_BPIX.png")
    c = plotHitsPerModule(reco_moduleHistsBPIX_L3,dir+"/nRecoHitsPerModule_L3_BPIX.png")
    c = plotHitsPerModule(reco_moduleHistsBPIX_L4,dir+"/nRecoHitsPerModule_L4_BPIX.png")
    c = plotHitsPerModule(reco_moduleHistsFPIX_L1,dir+"/nRecoHitsPerModule_L1_FPIX.png")
    c = plotHitsPerModule(reco_moduleHistsFPIX_L2,dir+"/nRecoHitsPerModule_L2_FPIX.png")
    c = plotHitsPerModule(reco_moduleHistsFPIX_L3,dir+"/nRecoHitsPerModule_L3_FPIX.png")


    c = r.TCanvas("c", "canvas", 900, 700)
    c.cd()
    test.SetLineColor(r.kRed)
    test.Draw("HIST")
    paveCMS.Draw("same")
    l0.Draw("same")
    c.SaveAs(dir+'/test.png')  


if __name__ == '__main__':
    options, paths = parse_arguments()
    main(options = options, paths = paths)