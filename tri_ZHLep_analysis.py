#!/usr/bin/env python

# Author: Jonathon Langford
#         Imperial College London
#         CMS, Hgg IC group

# Description: To perform trilinear analysis on ZHLeptonic tagged events
#              Use reco information to tag events with diphoton, dilepton consistent with ZH kinematics
#              		> The lepton selection is lacking: 2nd bullet point in PAS for ttHLeptonic Category, discuss with Nick how to implement this
#			  This incluse: Loose requirements on electron and tight on Muon
#              Output pT(H) and m_ZH distribution: reco and gen-level
#              		> Compatible with extract C1 values using event re-weighting
#              m_gg distribution in each bin of kinematic distribution
#              To be run on generated MC samples: signal + background (expect background to be very small)

#		UPDATE: Now use gen-level Higgs (instead of diphoton pair as Delphes card does not extract genPhotons)
###############################################################################
#	PRELIMINARIES

#Import Libraries
import sys
import ROOT
import math
from array import array

#Check for correct input
if len(sys.argv) != 3:
  print " Usage: Example1.py input_file input_file_rwgt output_file"
  sys.exit(1)

#Load delphes ROOT libraries
ROOT.gSystem.Load("libDelphes")

#Including packages to read TTree
try:
  ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
  ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
except:
  pass

#Take input file from command line
inputFile = sys.argv[1]
inputFile_rwgt = sys.argv[2]

# Create chain of root trees
chain = ROOT.TChain("Delphes")
chain_rwgt = ROOT.TChain("Delphes")
chain.Add(inputFile)
chain_rwgt.Add(inputFile_rwgt)

# Create objects of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
treeReader_rwgt = ROOT.ExRootTreeReader(chain_rwgt)
numberOfEntries = treeReader.GetEntries()
#numberOfEntries = 10000

# Get pointers to branches used in this analysis
branchEvent = treeReader.UseBranch("Event")
branchEvent_rwgt = treeReader_rwgt.UseBranch("Event")
branchGenParticle = treeReader.UseBranch("Particle")
branchPhoton =  treeReader.UseBranch("Photon")
branchElectron = treeReader.UseBranch("Electron")
branchMuon = treeReader.UseBranch("Muon")


##############################################################################
#	CONFIGURE OUTPUT
#Open .root file to write histograms to
f = ROOT.TFile.Open( "output_ZHLep_1e5.root" ,"RECREATE")

# Book histograms

#pT_H
#Gen-level
hist_pTH_gen_LO = ROOT.TH1F("h_pTH_gen_LO", "LO Diphoton p_{T} (gen level)", 30, 0, 300 )
hist_pTH_gen_O3 = ROOT.TH1F("h_pTH_gen_O3", "O3 Diphoton p_{T} (gen level) (rwgt)", 30, 0, 300 )
hist_pTH_gen_LO_lb = ROOT.TH1F("h_pTH_gen_LO_lb", "LO Diphoton p_{T} (gen-level)", 6, 0, 300 )
hist_pTH_gen_O3_lb = ROOT.TH1F("h_pTH_gen_O3_lb", "O3 Diphoton p_{T} (gen-level) (rwgt)", 6, 0, 300 )
#Reco-level
hist_pTH_reco_LO = ROOT.TH1F("h_pTH_reco_LO", "LO Diphoton p_{T} (reco level)", 30, 0, 300 )
hist_pTH_reco_O3 = ROOT.TH1F("h_pTH_reco_O3", "O3 Diphoton p_{T} (reco level) (rwgt)", 30, 0, 300 )
hist_pTH_reco_LO_lb = ROOT.TH1F("h_pTH_reco_LO_lb", "LO Diphoton p_{T} (reco-level)", 6, 0, 300 )
hist_pTH_reco_O3_lb = ROOT.TH1F("h_pTH_reco_O3_lb", "O3 Diphoton p_{T} (reco-level) (rwgt)", 6, 0, 300 )

#m_llgg
#Gen-level
hist_mZH_gen_LO = ROOT.TH1F("h_mZH_gen_LO", "LO m_{ZH} (gen level)", 30, 200, 500 )
hist_mZH_gen_O3 = ROOT.TH1F("h_mZH_gen_O3", "O3 m_{ZH} (gen level) (rwgt)", 30, 200, 500 )
hist_mZH_gen_LO_lb = ROOT.TH1F("h_mZH_gen_LO_lb", "LO m_{ZH} (gen-level)", 6, 200, 500 )
hist_mZH_gen_O3_lb = ROOT.TH1F("h_mZH_gen_O3_lb", "O3 m_{ZH} (gen-level) (rwgt)", 6, 200, 500 )
#Reco-level
hist_mZH_reco_LO = ROOT.TH1F("h_mZH_reco_LO", "LO m_{ZH} (reco level)", 30, 200, 500 )
hist_mZH_reco_O3 = ROOT.TH1F("h_mZH_reco_O3", "O3 m_{ZH} (reco level) (rwgt)", 30, 200, 500 )
hist_mZH_reco_LO_lb = ROOT.TH1F("h_mZH_reco_LO_lb", "LO m_{ZH} (reco-level)", 6, 200, 500 )
hist_mZH_reco_O3_lb = ROOT.TH1F("h_mZH_reco_O3_lb", "O3 m_{ZH} (reco-level) (rwgt)", 6, 200, 500 )

#m_gg
hist_mgg_reco = ROOT.TH1F("h_mgg_reco", "Diphoton invariant mass spectrum (reco level)", 30, 110, 140 )
hist_mgg_gen = ROOT.TH1F("h_mgg_gen", "Diphoton invariant mass spectrum (gen level)", 30, 110, 140 )
#m_ll
hist_mll_reco = ROOT.TH1F("h_mll_reco", "Dilepton invariant mass spectrum (reco level)", 50, 65, 115 )
hist_mll_gen = ROOT.TH1F("h_mll_gen", "Dilepton invariant mass spectrum (gen level)", 50, 65, 115 )

#m_gg: In bins of pTH
hist_mgg_reco_1 = ROOT.TH1F("h_mgg_reco_1", "m_{#gamma#gamma}^{reco}  for  p_{T}^{reco}(#gamma#gamma) #in [0,50] GeV", 30, 110, 140 )
hist_mgg_reco_2 = ROOT.TH1F("h_mgg_reco_2", "m_{#gamma#gamma}^{reco}  for  p_{T}^{reco}(#gamma#gamma) #in [50,100] GeV", 30, 110, 140 )
hist_mgg_reco_3 = ROOT.TH1F("h_mgg_reco_3", "m_{#gamma#gamma}^{reco}  for  p_{T}^{reco}(#gamma#gamma) #in [100,150] GeV", 30, 110, 140 )
hist_mgg_reco_4 = ROOT.TH1F("h_mgg_reco_4", "m_{#gamma#gamma}^{reco}  for  p_{T}^{reco}(#gamma#gamma) #in [150,200] GeV", 30, 110, 140 )
hist_mgg_reco_5 = ROOT.TH1F("h_mgg_reco_5", "m_{#gamma#gamma}^{reco}  for  p_{T}^{reco}(#gamma#gamma) #in [200,250] GeV", 30, 110, 140 )
hist_mgg_reco_6 = ROOT.TH1F("h_mgg_reco_6", "m_{#gamma#gamma}^{reco}  for  p_{T}^{reco}(#gamma#gamma) #in [250,350] GeV", 30, 110, 140 )

#Response matrix plots
hist_pTH_responseMatrix = ROOT.TH2F("h_pTH_responseMatrix","; p_{T}^{reco}(#gamma#gamma)   GeV; p_{T}^{gen}(#gamma#gamma)   GeV", 6, 0, 300, 6, 0, 300 )
hist_mZH_responseMatrix = ROOT.TH2F("h_mZH_responseMatrix","; m_{ZH}^{reco}   GeV; m_{ZH}^{gen}   GeV", 6, 200, 500, 6, 200, 500 )
###############################################################################
#	FUNCTIONS FOR KINEMATICS

def deltaR( eta1, phi1, eta2, phi2 ):
  return math.sqrt( (eta1-eta2)*(eta1-eta2) + (phi1-phi2)*(phi1-phi2) )

def pT_vector_calc( part1, part2 ):
  Px1 = part1.PT*math.cos( part1.Phi )
  Px2 = part2.PT*math.cos( part2.Phi )
  Py1 = part1.PT*math.sin( part1.Phi )
  Py2 = part2.PT*math.sin( part2.Phi )
  return math.sqrt( (Px1+Px2)*(Px1+Px2) + (Py1+Py2)*(Py1+Py2) )

###############################################################################
#	FUNCTIONS FOR EVENT SELECTION

def SelectPhoton( _photon, photonPtThreshold, photonEtaThresholds, phoIsoChRelThreshold ):
  photon_pass = True
  if( _photon.PT < photonPtThreshold ): photon_pass = False
  #Eta: inc outside transition region between barrel and endcap 
  if( ( abs( _photon.Eta ) > photonEtaThresholds[2] ) | ( ( abs( _photon.Eta ) > photonEtaThresholds[0] ) & ( abs( _photon.Eta ) < photonEtaThresholds[1] ) ) ): photon_pass = False
  #Isolation: currently only using Ich, need Ipho and Itrk
  if( (_photon.SumPtCharged/_photon.PT) > phoIsoChRelThreshold ): photon_pass = False
  #SHOWER SHAPE VARIABLES: R9, sigma_etaeta, need to access from Delphes in some way

  return photon_pass


def SelectDiPhoton( _leadPhoton, _subleadPhoton, leadPhoPTOverMassThreshold, subleadPhoPTOverMassThreshold ):
  diphoton_pass = True
  q_gg = _leadPhoton.P4()+_subleadPhoton.P4()
  m_gg = math.sqrt( q_gg*q_gg )
  if( _leadPhoton.PT/m_gg < leadPhoPTOverMassThreshold ) | ( _subleadPhoton.PT/m_gg < subleadPhoPTOverMassThreshold ): diphoton_pass = False
  #DIPHOTON MVA EQUIVALENT
  return diphoton_pass


def SelectMuon( _muon, _dipho, muonPtThreshold, muonEtaThreshold, muPFIsoSumRelThreshold, deltaRMuonPhoThreshold ):
  muon_pass = True
  if( _muon.PT < muonPtThreshold ): muon_pass = False
  if( abs( _muon.Eta ) > muonEtaThreshold ): muon_pass = False
  #Vertex: missing, require vertex info in CMS card, copy isTightMuon() (see implementation on git)
  #Isolation: using sumPt variable: assuming same as hard sum in flashgg::LeptonSelection.cc
  if( (_muon.SumPt/_muon.PT) > muPFIsoSumRelThreshold ): muon_pass = False
  
  #if muon passed then calc dR between leadPho and subleadPho
  if muon_pass:
    dR_Muon_LeadPho = deltaR( _dipho[0][0].Eta, _dipho[0][0].Phi, _muon.Eta, _muon.Phi )
    dR_Muon_SubleadPho = deltaR( _dipho[0][1].Eta, _dipho[0][1].Phi, _muon.Eta, _muon.Phi )
    if( dR_Muon_LeadPho < deltaRMuonPhoThreshold ) | ( dR_Muon_SubleadPho < deltaRMuonPhoThreshold ): muon_pass = False

  return muon_pass


def SelectElectron( _electron, _dipho, electronPtThreshold , electronEtaThresholds, electronPhoMassThreshold , deltaRElectronPhoThreshold ):
  electron_pass = True
  if( _electron.PT < electronPtThreshold ): electron_pass = False
  #Eta: inc outside transition region between barrel and endcap 
  if( ( abs( _electron.Eta ) > electronEtaThresholds[2] ) | ( ( abs( _electron.Eta ) > electronEtaThresholds[0] ) & ( abs( _electron.Eta ) < electronEtaThresholds[1] ) ) ): electron_pass = False 
  #Vertex: missing, require vertex info
  #ID: flashgg::passLooseID()

  #mass of electron+photon not close to Z mass: fasely recon electrons
  if electron_pass:
    m_eLeadPho = math.sqrt( abs((_dipho[0][0].P4()+_electron.P4())*(_dipho[0][0].P4()+_electron.P4())) )
    m_eSubleadPho = math.sqrt( abs((_dipho[0][1].P4()+_electron.P4())*(_dipho[0][1].P4()+_electron.P4())) )
    if( abs( m_eLeadPho-91.2 ) < 5. ) | ( abs( m_eSubleadPho-91.2 ) < 5. ): electron_pass = False

  #if electron passed then calc dR between leadPho and subleadPho
  if electron_pass:
    dR_Electron_LeadPho = deltaR( _dipho[0][0].Eta, _dipho[0][0].Phi, _electron.Eta, _electron.Phi )
    dR_Electron_SubleadPho = deltaR( _dipho[0][1].Eta, _dipho[0][1].Phi, _electron.Eta, _electron.Phi )
    if( dR_Electron_LeadPho < deltaRElectronPhoThreshold ) | ( dR_Electron_SubleadPho < deltaRElectronPhoThreshold ): electron_pass = False

  return electron_pass

###############################################################################
#	COUNTERS FOR DEBUGGING
N_dipho = 0
N_diMuon = 0
N_diElectron = 0
N_selection = 0
N_Zll = 0

###############################################################################
#	EVENTS LOOP

# Loop over all events
for entry in range(0, numberOfEntries ):
  
  if entry % 10000 == 0: print "Processing event: (", entry, "/", numberOfEntries, ")"
  #############################################################################
  #Define boolean for event passing selection
  event_pass = False
  # Load selected branches with data from specified event
  treeReader.ReadEntry(entry)
  treeReader_rwgt.ReadEntry(entry)

  #############################################################################
  #Event branch: get event weight
  _event = branchEvent.At(0)
  LO_weight = _event.Weight

  _event_rwgt = branchEvent_rwgt.At(0)
  O3_weight = _event_rwgt.Weight

  #############################################################################
  #ZHLeptonicTag:
  #	> For now: Using gen-level photons. This needs to be changed to reco-level, how can I mimic diphotonMVA?
  #     > 2 same-flavour leptons at reco-level. Passing same selection requirements as in flashgg::ZHLeptonicTagDumper
  # 		> Use selectMuons and selectElectrons function defined above
  
  #list to hold photon and lepton candidates
  photons = []
  diphotons = []
  photon_pair = []
  muons = []
  electrons = []

  #booleans describing event passing different stages of selection
  photon_selection = False
  isDiMuon = False
  isDiElectron = False
  isZ = False

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Photon Selection
  if branchPhoton.GetEntries() > 0:

    #Loop over photons in event and apply photon selection
    for i in range( branchPhoton.GetEntries() ):
      photon = branchPhoton.At(i)
      #Apply selection on single photonsL Pt threshold currently relaxed to 20GeV
      if( SelectPhoton( photon, photonPtThreshold=20., photonEtaThresholds=[1.4442,1.566,2.5], phoIsoChRelThreshold=0.3 ) ): photons.append( photon )
  
    #if atleast 2 photons in event
    if len( photons ) >= 2:
    
      #sort photons according to pT (descending)
      photons.sort( key=lambda g: g.PT, reverse=True )

      #Loop over photon pairs in event and apply diphoton selection
      for leadPho_idx in range( len( photons ) ):
        for subleadPho_idx in range( len( photons ) ):
          #Only once for each pair
          if subleadPho_idx > leadPho_idx:
            if( SelectDiPhoton( photons[leadPho_idx], photons[subleadPho_idx], leadPhoPTOverMassThreshold=0.375, subleadPhoPTOverMassThreshold=0.25) ): diphotons.append( [photons[leadPho_idx],photons[subleadPho_idx]] )
      
      #If atleast one diphoton pair passing selection then set photon_selection to true
      #if >1 diphoton passing selection: choose pair with highest sum pT
      if( len( diphotons ) > 1 ):
        pT_max = -999.
        dipho_idx_opt = -999
        for dipho_idx in range( len(diphotons) ):
          pT_H = pT_vector_calc( diphotons[dipho_idx][0], diphotons[dipho_idx][1] )
          if pT_H > pT_max:
            pT_max = pT_H
            dipho_idx_opt = dipho_idx
        photon_pair.append( diphotons[ dipho_idx ] )
        photon_selection = True
      #else if = 1 then append photon_pair list
      elif( len( diphotons ) == 1 ):
        photon_pair.append( diphotons[0] )
        photon_selection = True

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Lepton Selection: to be performed if diphoton pair found
  if( photon_selection ):
    if branchMuon.GetEntries() > 0:
      #loop over Muons in event and extract those which satisfy criteria
      for i in range( branchMuon.GetEntries() ):
        muon = branchMuon.At(i)
      
        #Muon Selection
        if( SelectMuon( muon, photon_pair, muonPtThreshold=20., muonEtaThreshold=2.4, muPFIsoSumRelThreshold=0.25, deltaRMuonPhoThreshold=0.5 ) ):
          muons.append( muon )

    if branchElectron.GetEntries() > 0:
      #loop over Electrons in event and extract those which satisfy criteria
      for i in range( branchElectron.GetEntries() ):
        electron = branchElectron.At(i)

        #Electron selection
        if( SelectElectron( electron, photon_pair, electronPtThreshold=20., electronEtaThresholds=[1.4442,1.566,2.5], electronPhoMassThreshold=5., deltaRElectronPhoThreshold=1. ) ):
          electrons.append( electron )
  
  #check for size of vectors
  if len(muons) >= 2: isDiMuon = True
  if len(electrons) >= 2: isDiElectron = True

  #Z mass window: invariant mass of lepton pair
  if( isDiMuon ):
    q_ll = muons[0].P4()+muons[1].P4()
    m_ll = math.sqrt( q_ll*q_ll )
    if( m_ll > 70 ) & ( m_ll < 110 ): isZ = True
  if( isDiElectron ):
    q_ll = electrons[0].P4()+electrons[1].P4()
    m_ll = math.sqrt( q_ll*q_ll )
    if( m_ll > 70 ) & ( m_ll < 110 ): isZ = True

##############################################################################
  #Couters for debugging
  if( photon_selection ):
    N_dipho += 1
    if( isDiMuon ): N_diMuon += 1
    if( isDiElectron ): N_diElectron += 1
    if( isDiMuon | isDiElectron ): 
      N_selection += 1
      if( isZ ): N_Zll += 1
 
##############################################################################
#	EVENTS PASSING SELECTION
  if( photon_selection & ( isDiMuon | isDiElectron ) & isZ ):

    #Define final photons and leptons
    leadPhoton = photon_pair[0][0]
    subleadPhoton = photon_pair[0][1]
    if( isDiMuon ): 
      lep1 = muons[0] 
      lep2 = muons[1]
    else: 
      lep1 = electrons[0] 
      lep2 = electrons[1]   

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # GEN PARTICLE EXTRACTION
    dR_genPartleadPho_min = 9999.
    dR_genPartsubleadPho_min = 9999.
    dR_genPartLep1_min = 9999.
    dR_genPartLep2_min = 9999.

    #initialise
    genleadPho_idx = -999
    gensubleadPho_idx = -999
    genLep1_idx = -999
    genLep2_idx = -999

    if branchGenParticle.GetEntries() > 0:
    #loop over GenParticles in event
      for i in range( branchGenParticle.GetEntries() ):
        genPart = branchGenParticle.At(i)

        #Photon pair extraction: use PID and Status
        if( genPart.PID == 22 ):
          dR = deltaR( genPart.Eta, genPart.Phi, leadPhoton.Eta, leadPhoton.Phi ) # LeadPhoton
          if( dR < dR_genPartleadPho_min ):
            dR_genPartleadPho_min = dR
            genleadPho_idx = i
          dR = deltaR( genPart.Eta, genPart.Phi, subleadPhoton.Eta, subleadPhoton.Phi ) # SubleadPhoton
          if( dR < dR_genPartsubleadPho_min ):
            dR_genPartsubleadPho_min = dR
            gensubleadPho_idx = i
        
        #Lepton pair extraction
        if( isDiMuon ): #Muons
          if( abs( genPart.PID ) == 13 ):
            dR = deltaR( genPart.Eta, genPart.Phi, lep1.Eta, lep1.Phi ) # lepton 1
            if( dR < dR_genPartLep1_min ):
              dR_genPartLep1_min = dR
              genLep1_idx = i
            dR = deltaR( genPart.Eta, genPart.Phi, lep2.Eta, lep2.Phi ) # lepton 2
            if( dR < dR_genPartLep2_min ):
              dR_genPartLep2_min = dR
              genLep2_idx = i
        else: #Electrons
          if( abs( genPart.PID ) == 11 ):
            dR = deltaR( genPart.Eta, genPart.Phi, lep1.Eta, lep1.Phi ) # lepton 1
            if( dR < dR_genPartLep1_min ):
              dR_genPartLep1_min = dR
              genLep1_idx = i
            dR = deltaR( genPart.Eta, genPart.Phi, lep2.Eta, lep2.Phi ) # lepton 2
            if( dR < dR_genPartLep2_min ):
              dR_genPartLep2_min = dR
              genLep2_idx = i 

      leadPhoton_gen = branchGenParticle.At( genleadPho_idx )
      subleadPhoton_gen = branchGenParticle.At( gensubleadPho_idx )
      lep1_gen = branchGenParticle.At( genLep1_idx )
      lep2_gen = branchGenParticle.At( genLep2_idx )
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    

    #Calculate kinematics
    pT_H_reco = pT_vector_calc( leadPhoton, subleadPhoton )
    pT_H_gen = pT_vector_calc( leadPhoton_gen, subleadPhoton_gen )

    q_ggll_reco = leadPhoton.P4()+subleadPhoton.P4()+lep1.P4()+lep2.P4() 
    q_ggll_gen = leadPhoton_gen.P4()+subleadPhoton_gen.P4()+lep1_gen.P4()+lep2_gen.P4()
    m_ZH_reco = math.sqrt( q_ggll_reco*q_ggll_reco )
    m_ZH_gen = math.sqrt( q_ggll_gen*q_ggll_gen ) 

    m_gg_reco = math.sqrt( (leadPhoton.P4()+subleadPhoton.P4())*(leadPhoton.P4()+subleadPhoton.P4()) )
    m_gg_gen = math.sqrt( (leadPhoton_gen.P4()+subleadPhoton_gen.P4())*(leadPhoton_gen.P4()+subleadPhoton_gen.P4()) )
    m_ll_reco = math.sqrt( (lep1.P4()+lep2.P4())*(lep1.P4()+lep2.P4()) )
    m_ll_gen = math.sqrt( (lep1_gen.P4()+lep2_gen.P4())*(lep1_gen.P4()+lep2_gen.P4()) )

    #Fill histograms
    hist_pTH_reco_LO.Fill( pT_H_reco, LO_weight )
    hist_pTH_reco_O3.Fill( pT_H_reco, O3_weight )
    hist_pTH_reco_LO_lb.Fill( pT_H_reco, LO_weight )
    hist_pTH_reco_O3_lb.Fill( pT_H_reco, O3_weight )
    hist_pTH_gen_LO.Fill( pT_H_gen, LO_weight )
    hist_pTH_gen_O3.Fill( pT_H_gen, O3_weight )
    hist_pTH_gen_LO_lb.Fill( pT_H_gen, LO_weight )
    hist_pTH_gen_O3_lb.Fill( pT_H_gen, O3_weight )

    hist_mZH_reco_LO.Fill( m_ZH_reco, LO_weight )
    hist_mZH_reco_O3.Fill( m_ZH_reco, O3_weight )
    hist_mZH_reco_LO_lb.Fill( m_ZH_reco, LO_weight )
    hist_mZH_reco_O3_lb.Fill( m_ZH_reco, O3_weight )
    hist_mZH_gen_LO.Fill( m_ZH_gen, LO_weight )
    hist_mZH_gen_O3.Fill( m_ZH_gen, O3_weight )
    hist_mZH_gen_LO_lb.Fill( m_ZH_gen, LO_weight )
    hist_mZH_gen_O3_lb.Fill( m_ZH_gen, O3_weight )

    hist_mgg_reco.Fill( m_gg_reco, LO_weight )
    hist_mgg_gen.Fill( m_gg_gen, LO_weight )
    hist_mll_reco.Fill( m_ll_reco, LO_weight )
    hist_mll_gen.Fill( m_ll_gen, LO_weight )

    if( pT_H_reco >= 0. ) & ( pT_H_reco < 50. ): hist_mgg_reco_1.Fill( m_gg_reco, LO_weight )
    elif( pT_H_reco >= 50. ) & ( pT_H_reco < 100. ): hist_mgg_reco_2.Fill( m_gg_reco, LO_weight )
    elif( pT_H_reco >= 100. ) & ( pT_H_reco < 150. ): hist_mgg_reco_3.Fill( m_gg_reco, LO_weight )
    elif( pT_H_reco >= 150. ) & ( pT_H_reco < 200. ): hist_mgg_reco_4.Fill( m_gg_reco, LO_weight )
    elif( pT_H_reco >= 200. ) & ( pT_H_reco < 250. ): hist_mgg_reco_5.Fill( m_gg_reco, LO_weight )
    elif( pT_H_reco >= 250. ) & ( pT_H_reco < 300. ): hist_mgg_reco_6.Fill( m_gg_reco, LO_weight )

    #Fill response matrix
    hist_pTH_responseMatrix.Fill( pT_H_reco, pT_H_gen, LO_weight )
    hist_mZH_responseMatrix.Fill( m_ZH_reco, m_ZH_gen, LO_weight )
      
#
###############################################################################
#	FINAL OUTPUT CONFIG

#normalise matrix by column: i.e. see what percentage gen level falls in each reco bin
#loop over columns
for i in range(1,hist_pTH_responseMatrix.GetNbinsX()+1):
  column_sum = 0
  #loop over rows and sum up
  for j in range(1,hist_pTH_responseMatrix.GetNbinsY()+1):
    column_sum += hist_pTH_responseMatrix.GetBinContent( i, j )

  #loop over rows again and scale value
  for j in range(1, hist_pTH_responseMatrix.GetNbinsY()+1):
    hist_pTH_responseMatrix.SetBinContent( i, j, (hist_pTH_responseMatrix.GetBinContent(i,j)*100)/column_sum )

#loop over columns
for i in range(1,hist_mZH_responseMatrix.GetNbinsX()+1):
  column_sum = 0
  #loop over rows and sum up
  for j in range(1,hist_mZH_responseMatrix.GetNbinsY()+1):
    column_sum += hist_mZH_responseMatrix.GetBinContent( i, j )

  #loop over rows again and scale value
  for j in range(1, hist_mZH_responseMatrix.GetNbinsY()+1):
    hist_mZH_responseMatrix.SetBinContent( i, j, (hist_mZH_responseMatrix.GetBinContent(i,j)*100)/column_sum )

#	WRITE HISTOGRAMS TO FILE AND CLOSE
f.Write()
f.Close()

#Print out info
print "########################################################## "
print "                      CUT COUNTERS                         "
print "Total events:", numberOfEntries
print "	-> Diphoton events:", N_dipho
print "	-> DiMuon events:", N_diMuon
print "	-> DiElec events:", N_diElectron
print "	-> DiSele events:", N_selection
print "	-> DiSele in Z mass window events:", N_Zll
print "########################################################## "


raw_input("Press Enter to continue...")
