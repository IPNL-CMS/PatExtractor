Install instructions
===========================

Setup local area
----------------

### CMSSW

    export SCRAM_ARCH=slc6_amd64_gcc472
    cmsrel CMSSW_5_3_18
    cd CMSSW_5_3_18
    cmsenv

    cd src/

CMSSW dependencies
------------------

### PAT

    git cms-addpkg PhysicsTools/PatAlgos

### E/gamma tools

    git cms-addpkg EgammaAnalysis/ElectronTools
    cd EgammaAnalysis/ElectronTools/data
    cat download.url | xargs wget
    cd -

### PU Jet ID

    git cms-merge-topic -u IPNL-CMS:53x_pujetid

Install PatExtractor
--------------------

    mkdir Extractors
    cd Extractors
    git clone https://github.com/IPNL-CMS/PatExtractor.git
    cd PatExtracor
    git checkout cmssw_5_3_18
    cd ..

## Build

    scram b -j8
