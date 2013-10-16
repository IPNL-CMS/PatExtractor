import FWCore.ParameterSet.Config as cms

def extractScaleFactor(dict, key):
    s = []

    for pt, sf in dict[key].items():
        pt_array = pt.split("_", 1)

        pset = {
                "pt": cms.vdouble(float(pt_array[0]), float(pt_array[1])),
                "value": cms.double(float(sf["data/mc"]["efficiency_ratio"])),
                "error_low": cms.double(float(sf["data/mc"]["err_low"])),
                "error_high": cms.double(float(sf["data/mc"]["err_hi"]))
                }

        s.append(pset)

    return s

def loadMuonScaleFactor(filenameIso, filenameEff, effWorkingPoint, isoWorkingPoint):
  import pickle
  import math
  fIso = open(filenameIso)
  fEff = open(filenameEff)

  dictIso = pickle.load(fIso)
  dictEff = pickle.load(fEff)

  etas = [[0, 0.9], [0.9, 1.2], [1.2, 2.1], [2.1, 2.4]]

  mainSet = cms.VPSet()

  for eta in etas:
    if eta[0] == 0:
      key = "ptabseta<%.1f" % (eta[1])
    else:
      key = "ptabseta%.1f-%.1f" % (eta[0], eta[1])

    etaSetEff = extractScaleFactor(dictEff[effWorkingPoint], key)
    etaSetIso = extractScaleFactor(dictIso[isoWorkingPoint], key) #"combRelIsoPF04dBeta<012_Tight"

    assert len(etaSetEff) == len(etaSetIso)

    etaSet = cms.VPSet()

    for i in range(0, len(etaSetEff)):
        pset = cms.PSet(
                pt = etaSetIso[i]["pt"],
                value = cms.double(etaSetIso[i]["value"]._value * etaSetEff[i]["value"]._value),
                error_high = cms.double(math.sqrt( math.pow(etaSetIso[i]["error_high"]._value, 2) + math.pow(etaSetEff[i]["error_high"]._value, 2)  )),
                error_low = cms.double(math.sqrt( math.pow(etaSetIso[i]["error_low"]._value, 2) + math.pow(etaSetEff[i]["error_low"]._value, 2)  )),
                )

        etaSet.append(pset)

    mainSet.append(
        cms.PSet(
          eta = cms.vdouble(eta[0], eta[1]),
          SF = etaSet
        )
      )

  fIso.close()
  fEff.close()

  return mainSet

def loadElectronScaleFactor(filename):
  import json
  with open(filename) as f:
    data = json.load(f)

    etaBins = data["eta"]
    ptBins = data["pt"]

    mainSet = cms.VPSet()

    for i in range(0, len(etaBins) - 1):
      etaSet = cms.VPSet()

      for j in range(0, len(ptBins) - 1):

        pset = cms.PSet(
            pt = cms.vdouble(float(ptBins[j]), float(ptBins[j + 1])),
            value = cms.double(float(data["sf"][i][j][0])),
            error_high = cms.double(float(data["sf"][i][j][1])),
            error_low = cms.double(float(data["sf"][i][j][2])),
            )

        etaSet.append(pset)

      mainSet.append(
          cms.PSet(
            eta = cms.vdouble(float(etaBins[i]), float(etaBins[i + 1])),
            SF = etaSet
            )
          )

    return mainSet

def loadBTagScaleFactors(process):
  #process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDB2013")
  process.load ("RecoBTag.PerformanceDB.BTagPerformanceDB2013")
