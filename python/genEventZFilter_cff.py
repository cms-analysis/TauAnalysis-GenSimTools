import FWCore.ParameterSet.Config as cms

from TauAnalysis.GenSimTools.genEventZFilter_cfi import * 


genEventZeeFilter = genEventZFilter.clone()
genEventZeeFilter.Z0DaughtersPdgId.append( 11 )

genEventZmumuFilter = genEventZFilter.clone()
genEventZmumuFilter.Z0DaughtersPdgId.append( 13 )

genEventZtautauFilter = genEventZFilter.clone()
genEventZtautauFilter.Z0DaughtersPdgId.append( 15 )

pathZee = cms.Path(
    genEventZeeFilter
    )
pathZmumu = cms.Path(
    genEventZmumuFilter
    )
pathZtautau = cms.Path(
    genEventZtautauFilter
    )
