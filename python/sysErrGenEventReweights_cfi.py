import FWCore.ParameterSet.Config as cms

# produce event (re)weights for estimating PDF uncertainties (maximum is 3)
# (copied from ElectroWeakAnalysis/Utilities/test/PdfSystematicsAnalyzer.py;
#  as soon as possible to be replaced by "officially blessed" configuration parameterset)
pdfWeights = cms.EDProducer("PdfWeightProducer",
    PdfInfoTag = cms.untracked.InputTag("generator"),
    PdfSetNames = cms.untracked.vstring(
        "cteq66.LHgrid",
        #"MSTW2008nlo68cl.LHgrid",
        #"NNPDF20_100.LHgrid"                      
    )
)

# produce event (re)weights to estimate ISR/FSR uncertainties
# (copied from ElectroWeakAnalysis/Utilities/test/SimpleSystematicsAnalyzer.py;
#  as soon as possible to be replaced by "officially blessed" configuration parameterset)
isrWeight = cms.EDProducer("ISRWeightProducer",
    GenTag = cms.untracked.InputTag("genParticles"),
    ISRBinEdges = cms.untracked.vdouble(
         0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,
        10., 11., 12., 13., 14., 15., 16., 17., 18., 19.,
        20., 21., 22., 23., 24., 25., 26., 27., 28., 29.,
        30., 31., 32., 33., 34., 35., 36., 37., 38., 39.,
        40., 41., 42., 43., 44., 45., 46., 47., 48., 49.,
        999999.
    ),
    PtWeights = cms.untracked.vdouble( 
        0.800665, 0.822121, 0.851249, 0.868285, 0.878733,
        0.953853, 0.928108, 0.982021, 1.00659 , 1.00648 ,
        1.03218 , 1.04924 , 1.03621 , 1.08743 , 1.01951 ,
        1.10519 , 0.984263, 1.04853 , 1.06724 , 1.10183 ,
        1.0503  , 1.13162 , 1.03837 , 1.12936 , 0.999173,
        1.01453 , 1.11435 , 1.10545 , 1.07199 , 1.04542 ,
        1.00828 , 1.0822  , 1.09667 , 1.16144 , 1.13906 ,
        1.27974 , 1.14936 , 1.23235 , 1.06667 , 1.06363 ,
        1.14225 , 1.22955 , 1.12674 , 1.03944 , 1.04639 ,
        1.13667 , 1.20493 , 1.09349 , 1.2107  , 1.21073
    )
)

fsrWeight = cms.EDProducer("FSRWeightProducer",
    GenTag = cms.untracked.InputTag("genParticles")
)

produceSysErrGenEventReweights = cms.Sequence(
    ##pdfWeights
    isrWeight
   + fsrWeight
)
