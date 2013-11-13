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

def loadElectronScaleFactor(filename, workingPoint):
  import json
  with open(filename) as f:
    data = json.load(f)

    data = data[workingPoint]
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
  process.load ("RecoBTag.PerformanceDB.BTagPerformanceDBWinter13")

def loadLightJetsScaleFactor():
    # TF1 for each eta bin
    return cms.VPSet(
            cms.PSet(
                eta = cms.vdouble(0, 0.8),
                value = cms.string("((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x)))"),
                error_high = cms.string("((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x)))"),
                error_low = cms.string("((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x)))"),
                ),
            cms.PSet(
                eta = cms.vdouble(0.8, 1.6),
                value = cms.string("((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)))"),
                error_high = cms.string("((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x)))"),
                error_low = cms.string("((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x)))"),
                ),
            cms.PSet(
                eta = cms.vdouble(1.6, 2.4),
                value = cms.string("((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)))"),
                error_high = cms.string("((1.15575+(0.000693344*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x)))"),
                error_low = cms.string("((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x)))"),
                ),
            )
