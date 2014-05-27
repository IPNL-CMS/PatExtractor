## Dependencies

### PAT

    git cms-addpkg PhysicsTools/PatAlgos

### E/gamma tools

    git cms-addpkg EgammaAnalysis/ElectronTools
    cd EgammaAnalysis/ElectronTools/data
    cat download.url | xargs wget
    cd -

### PU Jet ID

    git cms-merge-topic -u IPNL-CMS:53x_pujetid
