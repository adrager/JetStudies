#! /usr/bin/python

import os

config="""
#
# Configuration File for the Calibration Program
# Hamburg, 2007/08/15
#

#---------------------------------------------------------------------------------
#   Mode
#---------------------------------------------------------------------------------
# Running mode: Calibration (0), Jet smearing (1)
Mode = 0


#---------------------------------------------------------------------------------
#Fit method: LVMini=1, use existing calibration=2, write start values=3
Fit method = 2
#Parametrization
Parametrization Class = L2L3JetParametrization

Number of IO Threads = 5

#Error Parametrization
tower error parametrization = const
jet error parametrization   = jet et
#jet error parametrization   = toy

start values = 1.0
jet start values        = 1. 0. 0. 0. 0. 0.
global jet start values = 0.998 4.997 3.084 2.048

#---------------------------------------------------------------------------------
#   Geometry / Binning for fitting
#---------------------------------------------------------------------------------
maximum eta twr used  = 82
granularity in eta    = 2   # allowed values are: 1,3,5,11,21,41 (*2)
granularity in phi    = 1   # allowed values are: 1,2,3,6,9,18,36,72
symmetry in eta       = false

jet granularity in eta = 82 #  allowed values are: 1,3,5,11,21,41 (*2)
jet granularity in phi = 1 #   1 : default

track granularity in eta = 1
track granularity in phi = 1

#---------------------------------------------------------------------------------
#   Analyses (GammaJet, TrackTower, ... )
#---------------------------------------------------------------------------------
Et cut on jet              = 0.0
Et cut on gamma            = 0.0
Et cut on Z                = 20.0
Et cut on tower            = 0.0
Et cut on cluster          = 0.0
Et cut on track            = 0.0
Et cut on n+1 Jet          = 0.0
Eta max cut on jet         = 6.0
Relative Rest Jet Cut      = 0.2      #NonLeadingJetsEt / PhotonEt
Min had fraction           = -0.05    #Default: 0.07
Max had fraction           = 9999    #Default: 0.95
Et genJet min              = 0.0
Et genJet max              = 7000.0
DeltaR cut on jet matching = 0.25
Max cut on relative n+1 Jet Et = 0.2 ## to get more dijet-like events

#---------------------------------------------------------------------------------
#   Input / Output
#---------------------------------------------------------------------------------
Number of IO Threads = 5

#jet correction source = JetMETCor
#jet correction name   = Spring10_AK5TRK
Default Tree Name      = CalibTree

# List of input files:
Gamma-Jet tree         = GammaJetTree
Z-Jet tree             = ZJetTree
Track-Tower tree       = TrackTowerTree
Track-Cluster tree     = TrackClusterTree
Di-Jet tree            = DiJetTree
Di-Jet Control1 tree   = DiJetTree
Di-Jet Control2 tree   = DiJetTree
Tri-Jet tree           = TriJetTree
Top tree               = TopTree

#Di-Jet input file = input/dijetlistspring10Track
#Di-Jet input file =  input/toy_dijet_const_5-500_uniform.root

#Top input file = /scratch/current/cms/user/stadie/Top_Madgraph.root
#Z-Jet input file = /scratch/current/cms/user/stadie/ZJet_Track_230_300_rereco.root; /scratch/current/cms/user/stadie/ZJet_Track_300_INF_rereco.root; 

# -1: use all events, 1000: run over 1000 events
use Gamma-Jet events       = 0
use Track-Tower events     = 0
use Track-Cluster events   = 0
use Di-Jet events          = -1
use Di-Jet Control1 events = -1
use Di-Jet Control2 events = -1
use Tri-Jet events         = 0
use Z-Jet events           = 0
use Top events             = 0

Gamma-Jet data class     = 1
Z-Jet data class     = 3
Di-Jet data class    = 11
Top data class       = 1


#---------------------------------------------------------------------------------
#   Event processor setup
#---------------------------------------------------------------------------------

PUTruthReweighting = false

PU TruthWeighting Reweight all eventvectors (for MC validation) = true
PU TruthWeighting = adraeger/kalibri/PUDistributions/Cert_2012_190456-208686_ReReco
UseNaf2.0Directories = true
PU TruthWeighting MC distribution = adraeger/kalibri/PUDistributions/TrueDistributions/Summer12S10CMSSW53
Di-Jet trigger names = HLT_PFJet400 ## will reweight to PUTruth distribution of unprescaled PF400-trigger as determined by Kristin
Di-Jet trigger thresholds = 5 ## set threshold low so all JetTruthEvents should pass



#Di-Jet prescale = 1000
#-----------------------------------------------------------------
#  Control plots
#-----------------------------------------------------------------
#  General parameters
create plots                     = true
#plots output directory           = L2L3Plots
#plots format                      = pdf

# JetTruthEvent plots
create JetTruthEvent plots    =  true

JetTruthEvent plots names =  MCTruthResponseVsNPUTruth; MCTruthResponseVsNPUTruthEta0; MCTruthResponseVsNPUTruthEta1; MCTruthResponseVsNPUTruthEta2; MCTruthResponseVsNPUTruthEta3; MCTruthResponseVsNPUTruthEta4; MCTruthResponseVsGenJetPt; MCTruthResponseVsEta; MCTruthResolPUEta0; MCTruthResolPUEta1; MCTruthResolPUEta2; MCTruthResolPUEta3; MCTruthResolPUEta4; MCTruthResolPU; MCTruthResolPUAll; MCTruthResolPUEta00;MCTruthResolPUEta01;MCTruthResolPUEta02;MCTruthResolPUEta03;MCTruthResolPUEta04;MCTruthResolPUEta05;MCTruthResolPUEta06;MCTruthResolPUEta07;MCTruthResolPUEta08;MCTruthResolPUEta09;MCTruthResolPUEta10;MCTruthResolPUEta11;MCTruthResolPUEta12;MCTruthResolPUEta13;MCTruthResolPUEta14;MCTruthResolPUEta15;MCTruthResolPUEta16;MCTruthResolPUEta17;MCTruthResolPUEta18;MCTruthResolPUEta19;MCTruthResolPUEta20;MCTruthResolPUEta21;MCTruthResolPUEta22;MCTruthResolPUEta23;MCTruthResolPUEta24;MCTruthResolPUEta25;MCTruthResolPUEta26;MCTruthResolPUEta27;MCTruthResolPUEta28
MCTruthResponseVsNPUTruth x variable        =  NPUTruth
MCTruthResponseVsNPUTruth x edges           =  50 0 50
MCTruthResponseVsNPUTruth y variable        =  GenJetResponse
MCTruthResponseVsNPUTruth y edges           =  51 0 2.0 0.9 1.1 0.0 0.5
MCTruthResponseVsNPUTruth bin variable      =  Eta
MCTruthResponseVsNPUTruth bin edges         =  -5.0 5.0
MCTruthResponseVsNPUTruth correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResponseVsNPUTruth profile types     =  GaussFitMean; GaussFitWidth
MCTruthResponseVsNPUTruth legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResponseVsNPUTruth input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResponseVsNPUTruthEta0 x variable        =  NPUTruth
MCTruthResponseVsNPUTruthEta0 x edges           =  50 0 50
MCTruthResponseVsNPUTruthEta0 y variable        =  GenJetResponse
MCTruthResponseVsNPUTruthEta0 y edges           =  51 0 2.0 0.9 1.1 0.0 0.5
MCTruthResponseVsNPUTruthEta0 bin variable      =  AbsEta
MCTruthResponseVsNPUTruthEta0 bin edges         =  0 .5
MCTruthResponseVsNPUTruthEta0 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResponseVsNPUTruthEta0 profile types     =  GaussFitMean; GaussFitWidth
MCTruthResponseVsNPUTruthEta0 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResponseVsNPUTruthEta0 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResponseVsNPUTruthEta1 x variable        =  NPUTruth
MCTruthResponseVsNPUTruthEta1 x edges           =  50 0 50
MCTruthResponseVsNPUTruthEta1 y variable        =  GenJetResponse
MCTruthResponseVsNPUTruthEta1 y edges           =  51 0 2.0 0.9 1.1 0.0 0.5
MCTruthResponseVsNPUTruthEta1 bin variable      =  AbsEta
MCTruthResponseVsNPUTruthEta1 bin edges         =  .5 1.1
MCTruthResponseVsNPUTruthEta1 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResponseVsNPUTruthEta1 profile types     =  GaussFitMean; GaussFitWidth
MCTruthResponseVsNPUTruthEta1 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResponseVsNPUTruthEta1 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResponseVsNPUTruthEta2 x variable        =  NPUTruth
MCTruthResponseVsNPUTruthEta2 x edges           =  50 0 50
MCTruthResponseVsNPUTruthEta2 y variable        =  GenJetResponse
MCTruthResponseVsNPUTruthEta2 y edges           =  51 0 2.0 0.9 1.1 0.0 0.5
MCTruthResponseVsNPUTruthEta2 bin variable      =  AbsEta
MCTruthResponseVsNPUTruthEta2 bin edges         =  1.1 1.7
MCTruthResponseVsNPUTruthEta2 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResponseVsNPUTruthEta2 profile types     =  GaussFitMean; GaussFitWidth
MCTruthResponseVsNPUTruthEta2 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResponseVsNPUTruthEta2 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResponseVsNPUTruthEta3 x variable        =  NPUTruth
MCTruthResponseVsNPUTruthEta3 x edges           =  50 0 50
MCTruthResponseVsNPUTruthEta3 y variable        =  GenJetResponse
MCTruthResponseVsNPUTruthEta3 y edges           =  51 0 2.0 0.9 1.1 0.0 0.5
MCTruthResponseVsNPUTruthEta3 bin variable      =  AbsEta
MCTruthResponseVsNPUTruthEta3 bin edges         =  1.7 2.3
MCTruthResponseVsNPUTruthEta3 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResponseVsNPUTruthEta3 profile types     =  GaussFitMean; GaussFitWidth
MCTruthResponseVsNPUTruthEta3 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResponseVsNPUTruthEta3 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResponseVsNPUTruthEta4 x variable        =  NPUTruth
MCTruthResponseVsNPUTruthEta4 x edges           =  50 0 50
MCTruthResponseVsNPUTruthEta4 y variable        =  GenJetResponse
MCTruthResponseVsNPUTruthEta4 y edges           =  51 0 2.0 0.9 1.1 0.0 0.5
MCTruthResponseVsNPUTruthEta4 bin variable      =  AbsEta
MCTruthResponseVsNPUTruthEta4 bin edges         =  2.3 5.0
MCTruthResponseVsNPUTruthEta4 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResponseVsNPUTruthEta4 profile types     =  GaussFitMean; GaussFitWidth
MCTruthResponseVsNPUTruthEta4 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResponseVsNPUTruthEta4 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResponseVsGenJetPt x variable        =  GenJetPt;  log
MCTruthResponseVsGenJetPt x edges           =  50 10 3000
MCTruthResponseVsGenJetPt y variable        =  GenJetResponse
MCTruthResponseVsGenJetPt y edges           =  51 0 2.0 0.9 1.1 0.0 0.5
MCTruthResponseVsGenJetPt bin variable      =  Eta
MCTruthResponseVsGenJetPt bin edges         =  -5.0 -3.0 -1.3 1.3 3.0 5.0
MCTruthResponseVsGenJetPt correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResponseVsGenJetPt profile types     =  GaussFitMean; GaussFitWidth
MCTruthResponseVsGenJetPt legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResponseVsGenJetPt input samples     =  0:Z2star; #1:Z2; 2:Herwig

MCTruthResponseVsEta x variable         =  Eta
MCTruthResponseVsEta x edges            =  40 -5 5
MCTruthResponseVsEta y variable         =  GenJetResponse
MCTruthResponseVsEta y edges            =  51 0 2.0 0.9 1.1
MCTruthResponseVsEta bin variable       =  GenJetPt
MCTruthResponseVsEta bin edges          =  20 50 100 500 2000
MCTruthResponseVsEta correction types   =  Uncorrected; L2L3
#; L2L3L4 
MCTruthResponseVsEta profile types      =  GaussFitMean
MCTruthResponseVsEta legend label       =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResponseVsEta input samples     =  0:Z2star; #1:Z2; 2:Herwig

MCTruthResponsePU x variable         =  Eta
MCTruthResponsePU x edges            =  40 -5 5
MCTruthResponsePU y variable         =  GenJetResponse
MCTruthResponsePU y edges            =  51 0 2.0 0.9 1.1
MCTruthResponsePU bin variable       =  NPU
MCTruthResponsePU bin edges          =  0 1 10 20 30
MCTruthResponsePU cut variable       =  GenJetPt
MCTruthResponsePU cut edges          =  50 100
MCTruthResponsePU correction types   =  Uncorrected; L2L3
#; L2L3L4 
MCTruthResponsePU profile types      =  GaussFitMean
MCTruthResponsePU legend label       =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResponsePU input samples     =  0:Z2star; #1:Z2; 2:Herwig




MCTruthResolPUEta00 x variable        =  GenJetPt;  log
MCTruthResolPUEta00 x edges           =  25 10 3000
MCTruthResolPUEta00 y variable        =  GenJetResponse
MCTruthResolPUEta00 y edges           =  51 0 2.0 0 0.5
MCTruthResolPUEta00 bin variable      =  NPU
MCTruthResolPUEta00 bin edges         =  0 10 20 30 50
MCTruthResolPUEta00 cut variable      =  AbsEta
MCTruthResolPUEta00 cut edges         =  0 0.1
MCTruthResolPUEta00 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta00 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta00 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta00 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta01 x variable        =  GenJetPt;  log
MCTruthResolPUEta01 x edges           =  25 10 3001
MCTruthResolPUEta01 y variable        =  GenJetResponse
MCTruthResolPUEta01 y edges           =  51 0 2.0 0 0.5
MCTruthResolPUEta01 bin variable      =  NPU
MCTruthResolPUEta01 bin edges         =  0 10 20 30 50
MCTruthResolPUEta01 cut variable   =  AbsEta
MCTruthResolPUEta01 cut edges      =  0.1 0.2
MCTruthResolPUEta01 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta01 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta01 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta01 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta02 x variable        =  GenJetPt;  log
MCTruthResolPUEta02 x edges           =  25 10 3002
MCTruthResolPUEta02 y variable        =  GenJetResponse
MCTruthResolPUEta02 y edges           =  51 0 2.0 0 0.5
MCTruthResolPUEta02 bin variable      =  NPU
MCTruthResolPUEta02 bin edges         =  0 10 20 30 50
MCTruthResolPUEta02 cut variable   =  AbsEta
MCTruthResolPUEta02 cut edges      =  0.2 0.3
MCTruthResolPUEta02 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta02 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta02 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta02 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta03 x variable        =  GenJetPt;  log
MCTruthResolPUEta03 x edges           =  25 10 3003
MCTruthResolPUEta03 y variable        =  GenJetResponse
MCTruthResolPUEta03 y edges           =  51 0 2.0 0 0.5
MCTruthResolPUEta03 bin variable      =  NPU
MCTruthResolPUEta03 bin edges         =  0 10 20 30 50
MCTruthResolPUEta03 cut variable   =  AbsEta
MCTruthResolPUEta03 cut edges      =  0.3 0.4
MCTruthResolPUEta03 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta03 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta03 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta03 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta04 x variable        =  GenJetPt;  log
MCTruthResolPUEta04 x edges           =  25 10 3004
MCTruthResolPUEta04 y variable        =  GenJetResponse
MCTruthResolPUEta04 y edges           =  51 0 2.0 0 0.5
MCTruthResolPUEta04 bin variable      =  NPU
MCTruthResolPUEta04 bin edges         =  0 10 20 30 50
MCTruthResolPUEta04 cut variable   =  AbsEta
MCTruthResolPUEta04 cut edges      =  0.4 0.5
MCTruthResolPUEta04 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta04 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta04 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta04 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta05 x variable        =  GenJetPt;  log
MCTruthResolPUEta05 x edges           =  25 10 3005
MCTruthResolPUEta05 y variable        =  GenJetResponse
MCTruthResolPUEta05 y edges           =  51 0 2.0 0 0.5
MCTruthResolPUEta05 bin variable      =  NPU
MCTruthResolPUEta05 bin edges         =  0 10 20 30 50
MCTruthResolPUEta05 cut variable   =  AbsEta
MCTruthResolPUEta05 cut edges      =  0.5 0.6
MCTruthResolPUEta05 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta05 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta05 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta05 input samples     =  0:Z2star; #1:Z2; 2:Herwig



MCTruthResolPUEta06 x variable        =  GenJetPt;  log
MCTruthResolPUEta06 x edges           =  25 10 3006
MCTruthResolPUEta06 y variable        =  GenJetResponse
MCTruthResolPUEta06 y edges           =  51 0 2.0 0 0.5
MCTruthResolPUEta06 bin variable      =  NPU
MCTruthResolPUEta06 bin edges         =  0 10 20 30 50
MCTruthResolPUEta06 cut variable   =  AbsEta
MCTruthResolPUEta06 cut edges      =  0.6 0.7
MCTruthResolPUEta06 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta06 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta06 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta06 input samples     =  0:Z2star; #1:Z2; 2:Herwig



MCTruthResolPUEta07 x variable        =  GenJetPt;  log
MCTruthResolPUEta07 x edges           =  25 10 3007
MCTruthResolPUEta07 y variable        =  GenJetResponse
MCTruthResolPUEta07 y edges           =  51 0 2.0 0 0.5
MCTruthResolPUEta07 bin variable      =  NPU
MCTruthResolPUEta07 bin edges         =  0 10 20 30 50
MCTruthResolPUEta07 cut variable   =  AbsEta
MCTruthResolPUEta07 cut edges      =  0.7 0.8
MCTruthResolPUEta07 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta07 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta07 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta07 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta08 x variable        =  GenJetPt;  log
MCTruthResolPUEta08 x edges           =  25 10 3008
MCTruthResolPUEta08 y variable        =  GenJetResponse
MCTruthResolPUEta08 y edges           =  51 0 2.0 0 0.5
MCTruthResolPUEta08 bin variable      =  NPU
MCTruthResolPUEta08 bin edges         =  0 10 20 30 50
MCTruthResolPUEta08 cut variable   =  AbsEta
MCTruthResolPUEta08 cut edges      =  0.8 0.9
MCTruthResolPUEta08 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta08 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta08 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta08 input samples     =  0:Z2star; #1:Z2; 2:Herwig



MCTruthResolPUEta09 x variable        =  GenJetPt;  log
MCTruthResolPUEta09 x edges           =  25 10 3009
MCTruthResolPUEta09 y variable        =  GenJetResponse
MCTruthResolPUEta09 y edges           =  51 0 2.0 0 0.5
MCTruthResolPUEta09 bin variable      =  NPU
MCTruthResolPUEta09 bin edges         =  0 10 20 30 50
MCTruthResolPUEta09 cut variable   =  AbsEta
MCTruthResolPUEta09 cut edges      =  0.9 1.0
MCTruthResolPUEta09 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta09 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta09 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta09 input samples     =  0:Z2star; #1:Z2; 2:Herwig





MCTruthResolPUEta10 x variable        =  GenJetPt;  log
MCTruthResolPUEta10 x edges           =  25 10 3000
MCTruthResolPUEta10 y variable        =  GenJetResponse
MCTruthResolPUEta10 y edges           =  51 0 2.0 0 1.5
MCTruthResolPUEta10 bin variable      =  NPU
MCTruthResolPUEta10 bin edges         =  0 10 20 30 50
MCTruthResolPUEta10 cut variable   =  AbsEta
MCTruthResolPUEta10 cut edges      =  1.0 1.1
MCTruthResolPUEta10 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta10 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta10 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta10 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta11 x variable        =  GenJetPt;  log
MCTruthResolPUEta11 x edges           =  25 10 3001
MCTruthResolPUEta11 y variable        =  GenJetResponse
MCTruthResolPUEta11 y edges           =  51 0 2.0 0 1.5
MCTruthResolPUEta11 bin variable      =  NPU
MCTruthResolPUEta11 bin edges         =  0 10 20 30 50
MCTruthResolPUEta11 cut variable   =  AbsEta
MCTruthResolPUEta11 cut edges      =  1.1 1.2
MCTruthResolPUEta11 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta11 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta11 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta11 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta12 x variable        =  GenJetPt;  log
MCTruthResolPUEta12 x edges           =  25 10 3002
MCTruthResolPUEta12 y variable        =  GenJetResponse
MCTruthResolPUEta12 y edges           =  51 0 2.0 0 1.5
MCTruthResolPUEta12 bin variable      =  NPU
MCTruthResolPUEta12 bin edges         =  0 10 20 30 50
MCTruthResolPUEta12 cut variable   =  AbsEta
MCTruthResolPUEta12 cut edges      =  1.2 1.3
MCTruthResolPUEta12 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta12 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta12 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta12 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta13 x variable        =  GenJetPt;  log
MCTruthResolPUEta13 x edges           =  25 10 3003
MCTruthResolPUEta13 y variable        =  GenJetResponse
MCTruthResolPUEta13 y edges           =  51 0 2.0 0 1.5
MCTruthResolPUEta13 bin variable      =  NPU
MCTruthResolPUEta13 bin edges         =  0 10 20 30 50
MCTruthResolPUEta13 cut variable   =  AbsEta
MCTruthResolPUEta13 cut edges      =  1.3 1.4
MCTruthResolPUEta13 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta13 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta13 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta13 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta14 x variable        =  GenJetPt;  log
MCTruthResolPUEta14 x edges           =  25 10 3004
MCTruthResolPUEta14 y variable        =  GenJetResponse
MCTruthResolPUEta14 y edges           =  51 0 2.0 0 1.5
MCTruthResolPUEta14 bin variable      =  NPU
MCTruthResolPUEta14 bin edges         =  0 10 20 30 50
MCTruthResolPUEta14 cut variable   =  AbsEta
MCTruthResolPUEta14 cut edges      =  1.4 1.5
MCTruthResolPUEta14 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta14 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta14 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta14 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta15 x variable        =  GenJetPt;  log
MCTruthResolPUEta15 x edges           =  25 10 3005
MCTruthResolPUEta15 y variable        =  GenJetResponse
MCTruthResolPUEta15 y edges           =  51 0 2.0 0 1.5
MCTruthResolPUEta15 bin variable      =  NPU
MCTruthResolPUEta15 bin edges         =  0 10 20 30 50
MCTruthResolPUEta15 cut variable   =  AbsEta
MCTruthResolPUEta15 cut edges      =  1.5 1.6
MCTruthResolPUEta15 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta15 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta15 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta15 input samples     =  0:Z2star; #1:Z2; 2:Herwig



MCTruthResolPUEta16 x variable        =  GenJetPt;  log
MCTruthResolPUEta16 x edges           =  25 10 3006
MCTruthResolPUEta16 y variable        =  GenJetResponse
MCTruthResolPUEta16 y edges           =  51 0 2.0 0 1.5
MCTruthResolPUEta16 bin variable      =  NPU
MCTruthResolPUEta16 bin edges         =  0 10 20 30 50
MCTruthResolPUEta16 cut variable   =  AbsEta
MCTruthResolPUEta16 cut edges      =  1.6 1.7
MCTruthResolPUEta16 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta16 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta16 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta16 input samples     =  0:Z2star; #1:Z2; 2:Herwig



MCTruthResolPUEta17 x variable        =  GenJetPt;  log
MCTruthResolPUEta17 x edges           =  25 10 3007
MCTruthResolPUEta17 y variable        =  GenJetResponse
MCTruthResolPUEta17 y edges           =  51 0 2.0 0 1.5
MCTruthResolPUEta17 bin variable      =  NPU
MCTruthResolPUEta17 bin edges         =  0 10 20 30 50
MCTruthResolPUEta17 cut variable   =  AbsEta
MCTruthResolPUEta17 cut edges      =  1.7 1.8
MCTruthResolPUEta17 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta17 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta17 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta17 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta18 x variable        =  GenJetPt;  log
MCTruthResolPUEta18 x edges           =  25 10 3008
MCTruthResolPUEta18 y variable        =  GenJetResponse
MCTruthResolPUEta18 y edges           =  51 0 2.0 0 1.5
MCTruthResolPUEta18 bin variable      =  NPU
MCTruthResolPUEta18 bin edges         =  0 10 20 30 50
MCTruthResolPUEta18 cut variable   =  AbsEta
MCTruthResolPUEta18 cut edges      =  1.8 1.9
MCTruthResolPUEta18 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta18 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta18 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta18 input samples     =  0:Z2star; #1:Z2; 2:Herwig



MCTruthResolPUEta19 x variable        =  GenJetPt;  log
MCTruthResolPUEta19 x edges           =  25 10 3009
MCTruthResolPUEta19 y variable        =  GenJetResponse
MCTruthResolPUEta19 y edges           =  51 0 2.0 0 1.5
MCTruthResolPUEta19 bin variable      =  NPU
MCTruthResolPUEta19 bin edges         =  0 10 20 30 50
MCTruthResolPUEta19 cut variable   =  AbsEta
MCTruthResolPUEta19 cut edges      =  1.9 2.0
MCTruthResolPUEta19 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta19 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta19 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta19 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta20 x variable        =  GenJetPt;  log
MCTruthResolPUEta20 x edges           =  25 10 3000
MCTruthResolPUEta20 y variable        =  GenJetResponse
MCTruthResolPUEta20 y edges           =  51 0 2.0 0 2.5
MCTruthResolPUEta20 bin variable      =  NPU
MCTruthResolPUEta20 bin edges         =  0 10 20 30 50
MCTruthResolPUEta20 cut variable   =  AbsEta
MCTruthResolPUEta20 cut edges      =  2.0 2.1
MCTruthResolPUEta20 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta20 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta20 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta20 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta21 x variable        =  GenJetPt;  log
MCTruthResolPUEta21 x edges           =  25 10 3001
MCTruthResolPUEta21 y variable        =  GenJetResponse
MCTruthResolPUEta21 y edges           =  51 0 2.0 0 2.5
MCTruthResolPUEta21 bin variable      =  NPU
MCTruthResolPUEta21 bin edges         =  0 10 20 30 50
MCTruthResolPUEta21 cut variable   =  AbsEta
MCTruthResolPUEta21 cut edges      =  2.1 2.2
MCTruthResolPUEta21 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta21 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta21 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta21 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta22 x variable        =  GenJetPt;  log
MCTruthResolPUEta22 x edges           =  25 10 3002
MCTruthResolPUEta22 y variable        =  GenJetResponse
MCTruthResolPUEta22 y edges           =  51 0 2.0 0 2.5
MCTruthResolPUEta22 bin variable      =  NPU
MCTruthResolPUEta22 bin edges         =  0 10 20 30 50
MCTruthResolPUEta22 cut variable   =  AbsEta
MCTruthResolPUEta22 cut edges      =  2.2 2.3
MCTruthResolPUEta22 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta22 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta22 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta22 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta23 x variable        =  GenJetPt;  log
MCTruthResolPUEta23 x edges           =  25 10 3003
MCTruthResolPUEta23 y variable        =  GenJetResponse
MCTruthResolPUEta23 y edges           =  51 0 2.0 0 2.5
MCTruthResolPUEta23 bin variable      =  NPU
MCTruthResolPUEta23 bin edges         =  0 10 20 30 50
MCTruthResolPUEta23 cut variable   =  AbsEta
MCTruthResolPUEta23 cut edges      =  2.3 2.4
MCTruthResolPUEta23 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta23 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta23 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta23 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta24 x variable        =  GenJetPt;  log
MCTruthResolPUEta24 x edges           =  25 10 3004
MCTruthResolPUEta24 y variable        =  GenJetResponse
MCTruthResolPUEta24 y edges           =  51 0 2.0 0 2.5
MCTruthResolPUEta24 bin variable      =  NPU
MCTruthResolPUEta24 bin edges         =  0 10 20 30 50
MCTruthResolPUEta24 cut variable   =  AbsEta
MCTruthResolPUEta24 cut edges      =  2.4 2.5
MCTruthResolPUEta24 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta24 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta24 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta24 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta25 x variable        =  GenJetPt;  log
MCTruthResolPUEta25 x edges           =  25 10 3005
MCTruthResolPUEta25 y variable        =  GenJetResponse
MCTruthResolPUEta25 y edges           =  51 0 2.0 0 2.5
MCTruthResolPUEta25 bin variable      =  NPU
MCTruthResolPUEta25 bin edges         =  0 10 20 30 50
MCTruthResolPUEta25 cut variable   =  AbsEta
MCTruthResolPUEta25 cut edges      =  2.5 2.8
MCTruthResolPUEta25 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta25 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta25 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta25 input samples     =  0:Z2star; #1:Z2; 2:Herwig



MCTruthResolPUEta26 x variable        =  GenJetPt;  log
MCTruthResolPUEta26 x edges           =  25 10 3006
MCTruthResolPUEta26 y variable        =  GenJetResponse
MCTruthResolPUEta26 y edges           =  51 0 2.0 0 2.5
MCTruthResolPUEta26 bin variable      =  NPU
MCTruthResolPUEta26 bin edges         =  0 10 20 30 50
MCTruthResolPUEta26 cut variable   =  AbsEta
MCTruthResolPUEta26 cut edges      =  2.8 3.0
MCTruthResolPUEta26 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta26 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta26 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta26 input samples     =  0:Z2star; #1:Z2; 2:Herwig


MCTruthResolPUEta27 x variable        =  GenJetPt;  log
MCTruthResolPUEta27 x edges           =  25 10 3005
MCTruthResolPUEta27 y variable        =  GenJetResponse
MCTruthResolPUEta27 y edges           =  51 0 2.0 0 2.5
MCTruthResolPUEta27 bin variable      =  NPU
MCTruthResolPUEta27 bin edges         =  0 10 20 30 50
MCTruthResolPUEta27 cut variable   =  AbsEta
MCTruthResolPUEta27 cut edges      =  3.0 3.2
MCTruthResolPUEta27 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta27 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta27 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta27 input samples     =  0:Z2star; #1:Z2; 2:Herwig



MCTruthResolPUEta28 x variable        =  GenJetPt;  log
MCTruthResolPUEta28 x edges           =  25 10 3006
MCTruthResolPUEta28 y variable        =  GenJetResponse
MCTruthResolPUEta28 y edges           =  51 0 2.0 0 2.5
MCTruthResolPUEta28 bin variable      =  NPU
MCTruthResolPUEta28 bin edges         =  0 10 20 30 50
MCTruthResolPUEta28 cut variable   =  AbsEta
MCTruthResolPUEta28 cut edges      =  3.2 5.0
MCTruthResolPUEta28 correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUEta28 profile types     =  GaussFitWidth; RMS
MCTruthResolPUEta28 legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUEta28 input samples     =  0:Z2star; #1:Z2; 2:Herwig




## pileup inclusive
MCTruthResolPUAll x variable        =  GenJetPt;  log
MCTruthResolPUAll x edges           =  25 10 3000
MCTruthResolPUAll y variable        =  GenJetResponse
MCTruthResolPUAll y edges           =  51 0 2.0 0 0.5
MCTruthResolPUAll bin variable      =  NPU
MCTruthResolPUAll bin edges         =  0 50
MCTruthResolPUAll cut variable   =  AbsEta
MCTruthResolPUAll cut edges      =  0.0 9999
MCTruthResolPUAll correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPUAll profile types     =  GaussFitWidth; RMS
MCTruthResolPUAll legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPUAll input samples     =  0:Z2star; #1:Z2; 2:Herwig

MCTruthResolVsGenJetPt x variable        =  GenJetPt;  log
MCTruthResolVsGenJetPt x edges           =  50 10 3000
MCTruthResolVsGenJetPt y variable        =  GenJetResponse
MCTruthResolVsGenJetPt y edges           =  51 0 1.0 0 0.5
MCTruthResolVsGenJetPt bin variable      =  Eta
MCTruthResolVsGenJetPt bin edges         =  -5.0 -3.0 -1.3 1.3 3.0 5.0
MCTruthResolVsGenJetPt correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolVsGenJetPt profile types     =  GaussFitWidth
MCTruthResolVsGenJetPt legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolVsGenJetPt input samples     =  0:Z2star; #1:Z2; 2:Herwig

MCTruthRespFlavorVsGenJetPt x variable        =  GenJetPt;  log
MCTruthRespFlavorVsGenJetPt x edges           =  50 10 3000
MCTruthRespFlavorVsGenJetPt y variable        =  GenJetResponse
MCTruthRespFlavorVsGenJetPt y edges           =  51 0 2 0.9 1.1 0.9 1.1
MCTruthRespFlavorVsGenJetPt bin variable      =  Flavor
MCTruthRespFlavorVsGenJetPt bin edges         =  -0.5 0.5 1.5 2.5 3.5
MCTruthRespFlavorVsGenJetPt correction types  =  Uncorrected; L2L3
MCTruthRespFlavorVsGenJetPt profile types     =  Mean; GaussFitMean
MCTruthRespFlavorVsGenJetPt legend label      =  L2L3:CMS L1L2L3
MCTruthRespFlavorVsGenJetPt input samples     =  0:Z2star; #1:Z2; 2:Herwig

MCTruthResponseVsMeanWidth x variable         =  meanMoment
MCTruthResponseVsMeanWidth x edges            =  30 0 0.5
MCTruthResponseVsMeanWidth y variable         =  GenJetResponse
MCTruthResponseVsMeanWidth y edges            =  51 0 2 0.9 1.1 0.9 1.1 0.0 0.5
MCTruthResponseVsMeanWidth bin variable       =  GenJetPt
MCTruthResponseVsMeanWidth bin edges          =  10 30 50 80 120 300 600 2000
MCTruthResponseVsMeanWidth correction types   =  Uncorrected; L2L3
MCTruthResponseVsMeanWidth profile types      =  Mean; GaussFitMean; GaussFitWidth
MCTruthResponseVsMeanWidth legend label       =  L2L3:CMS L1L2L3
MCTruthResponseVsMeanWidth input samples     =  0:Z2star; #1:Z2; 2:Herwig

"""

jettypes = ["ak5FastPF"]#, "ak5PFCHS"]
#jettypes =  ["ak5PFCHS"]
#datadir = "/scratch/hh/dust/naf/cms/user/kriheine/CalibNTupel/MC/Z2star_pythia_v3" #Z2star
datadir = "/nfs/dust/cms/user/adraeger/kalibri/nTuples/MC/Z2star_pythia_v3" # NAF2.0
datadirmc = "/nfs/dust/cms/user/adraeger/kalibri/nTuples/MC/Z2_pythia_v3" # NAF2.0
datadir2mc = "/nfs/dust/cms/user/adraeger/kalibri/nTuples/MC/EE3C_herwigpp_v3" # NAF2.0
#datadirmc = "/scratch/hh/dust/naf/cms/user/kriheine/CalibNTupel/MC/Z2_pythia_v3" #Z2
#datadir2mc = "/scratch/hh/dust/naf/cms/user/kriheine/CalibNTupel/MC/EE3C_herwigpp_v3" #Herwig
jecname = "Dec12_Z2star"
datasetname="Z2star_TEST"

correctJets=False
CutAwayPU=False


for jettype in jettypes:
    print "make plots for jettype "+jettype
    if os.path.exists("tempdijetlist"):
        os.remove("tempdijetlist")
        os.remove("mcdijetlist")
##        os.remove("mc2dijetlist")
    if os.path.exists("tempplots"):
        os.system("rm tempplots/*")

    os.system("ls "+datadir+"/*_"+jettype+".root > tempdijetlist");
    os.system("ls "+datadirmc+"/*_"+jettype+".root > mcdijetlist");
    os.system("ls "+datadir2mc+"/*_"+jettype+".root > mc2dijetlist");
#    os.system("ls "+datadir+"*_"+jettype+"*.root > tempdijetlist");
#    os.system("ls "+datadirmc+"*_"+jettype+"*.root > mcdijetlist");
#    os.system("ls "+datadir2mc+"*_"+jettype+"*.root > mc2dijetlist");
    fcfg = open("valid.cfg", "w")
    #change labels
    if CutAwayPU:
        config = config.replace('CMS L1L2L3','CMS L2L3')

    fcfg.write(config)
    fcfg.write("plots output directory = tempplots\n")
    fcfg.write("Di-Jet input file = tempdijetlist\n")
 #   fcfg.write("Di-Jet Control1 input file = mcdijetlist\n")
 #   fcfg.write("Di-Jet Control2 input file = mc2dijetlist\n")
    jetalgo = jettype[0:3]
    
    if correctJets:
        fcfg.write("jet correction source = JetMETCor\n");
        fcfg.write("jet correction name   = "+jecname+"_"+jetalgo.upper()+jettype[3:len(jettype)]+"\n");
        #fcfg.write("jet correction name   = "+jecname+"_"+jetalgo.upper()+jettype[3:len(jettype)]+"NoOffset\n");
    if CutAwayPU:
        #fcfg.write("MAX n PU from MC = 8\n");
        fcfg.write("correct jets L1 = false\n")
        jec = jecname+'NoPU'
    else:
        fcfg.write("correct jets L1 = true\n")
        jec = jecname
    fcfg.close()
    
    kalibricmd = "./junk valid.cfg"
    print "running "+kalibricmd
    os.system(kalibricmd)
    tarball=jec+"plots"+jettype+".tar"
    tarcmd = "cd tempplots; tar cf ../"+tarball+" *Eta[0-9].eps *Pt[0-9].eps *Eta[0-9]_zoom.eps *Pt[0-9]_zoom.eps *Flavor[0-9].eps *Flavor[0-9]_zoom.eps *NPU[0-9].eps *NPU[0-9]_zoom.eps *.root; cd -"
    print "running "+tarcmd
    os.system(tarcmd)

print "please run:"
for jettype in jettypes:
    tarball=os.getcwd()+"/"+jec+"plots"+jettype+".tar"
    
    webcmd = "./scripts/createJECValidationHtmlPage.sh jetmet "+tarball+" \""+jec+"\" \""+jecname+"\" "+datasetname+" "+jettype[0:3] +" "+jettype[3:len(jettype)].lower()
    print webcmd
    #os.system(webcmd)
    

