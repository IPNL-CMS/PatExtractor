#! /bin/bash

# JEC
mkdir JECup
cp Extractor_MTT_*.py JECup/
cp createAndRunMCCrab.py JECup/
cp crab_MC.cfg.template.ipnl JECup/
cp Electron_scale_factors.json JECup/
cp Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.pkl JECup/
cp createOutputListForMC.py JECup/

sed -i 's/jesSign = 0/jesSign = 1/' JECup/Extractor_MTT_common.py
sed -i 's/dataset_name = dataset\[1\]/dataset_name = dataset\[1\] + "_JECup"/' JECup/createAndRunMCCrab.py

sed -i '54,69s/^/#/' JECup/createAndRunMCCrab.py

mkdir JECdown
cp Extractor_MTT_*.py JECdown/
cp createAndRunMCCrab.py JECdown/
cp crab_MC.cfg.template.ipnl JECdown/
cp Electron_scale_factors.json JECdown/
cp Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.pkl JECdown/
cp createOutputListForMC.py JECdown/

sed -i 's/jesSign = 0/jesSign = -1/' JECdown/Extractor_MTT_common.py
sed -i 's/dataset_name = dataset\[1\]/dataset_name = dataset\[1\] + "_JECdown"/' JECdown/createAndRunMCCrab.py

sed -i '54,69s/^/#/' JECdown/createAndRunMCCrab.py

# JER
mkdir JERup
cp Extractor_MTT_*.py JERup/
cp createAndRunMCCrab.py JERup/
cp crab_MC.cfg.template.ipnl JERup/
cp Electron_scale_factors.json JERup/
cp Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.pkl JERup/
cp createOutputListForMC.py JERup/

sed -i 's/jerSign = 0/jerSign = 1/' JERup/Extractor_MTT_common.py
sed -i 's/dataset_name = dataset\[1\]/dataset_name = dataset\[1\] + "_JERup"/' JERup/createAndRunMCCrab.py

sed -i '54,69s/^/#/' JERup/createAndRunMCCrab.py

mkdir JERdown
cp Extractor_MTT_*.py JERdown/
cp createAndRunMCCrab.py JERdown/
cp crab_MC.cfg.template.ipnl JERdown/
cp Electron_scale_factors.json JERdown/
cp Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.pkl JERdown/
cp createOutputListForMC.py JERdown/

sed -i 's/jerSign = 0/jerSign = -1/' JERdown/Extractor_MTT_common.py
sed -i 's/dataset_name = dataset\[1\]/dataset_name = dataset\[1\] + "_JERdown"/' JERdown/createAndRunMCCrab.py

sed -i '54,69s/^/#/' JERdown/createAndRunMCCrab.py
