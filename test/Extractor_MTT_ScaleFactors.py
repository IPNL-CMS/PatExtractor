import FWCore.ParameterSet.Config as cms

def loadMuonScaleFactor(filename):
  import pickle
  f = open(filename)

  dict = pickle.load(f)

  suffix = "2012ABCD"
  etas = [[0, 0.9], [0.9, 1.2], [1.2, 2.1]]

  mainSet = cms.VPSet()

  for eta in etas:
    if eta[0] == 0:
      key = "ptabseta<%.1f_%s" % (eta[1], suffix) 
    else:
      key = "ptabseta%.1f-%.1f_%s" % (eta[0], eta[1], suffix) 

    etaSet = cms.VPSet()

    for pt, sf in dict["Tight"][key].items():
      pt_array = pt.split("_", 1)

      pset = cms.PSet(
          pt = cms.vdouble(float(pt_array[0]), float(pt_array[1])),
          value = cms.double(float(sf["data/mc"]["efficiency_ratio"])),
          error_low = cms.double(float(sf["data/mc"]["err_low"])),
          error_high = cms.double(float(sf["data/mc"]["err_hi"]))
          )

      etaSet.append(pset)

    mainSet.append(
        cms.PSet(
          eta = cms.vdouble(eta[0], eta[1]),
          SF = etaSet
        )
      )

  return mainSet
