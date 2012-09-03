#!/bin/csh

################################################
#
# runPAT_extraction.csh
#
# Script launching the extraction of a PATified dataset
#
# This script works for any PATuple 
#
# Command is :
#
# source runPAT_extraction.csh P1 P2 P3 P4 P5 P6 
#
# Where parameters are defined as:
#
# P1 : do we run on MC or not?  True or False 
# P2 : the dataset type: HT,TTJets_TuneZ2_7TeV-madgraph-tauola,... 
# P3 : the dataset version: Fall11-PU_S6_START42_V14B-v1,Run2011A-05Aug2011-v1,...
# P4 : the CMSSW release: 4_4_1,...
# P5 : the GlobalTag: START44_V6,...
# P6 : launch the job to batch (P6=DO_BATCH) or not
#
# # Small tutorial on how to use the script:
#
# http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.PHYTuto
#
#
# Author: Seb Viret <viret@in2p3.fr>  (10/11/2011)
#
#################################################

# The directory where data is stored
set OUTPUT_ROOT = "$CASTOR_HOME/CMS/Extraction" 

# The directory where PATuples are retrieved
set INPUT_ROOT  = "$CASTOR_HOME/CMS/PAT"        

# Define your CMSSW working area (where CMSSW releases are installed)
set WORKDIR  = "$HOME/scratch0/testarea/"


##########################################################
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! # 
# YOU ARE NOT SUPPOSED TO TOUCH ANYTHING BELOW THIS LINE #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
##########################################################

set DATA_OUTPUT = $OUTPUT_ROOT"/"${2}_${3} 
set DATA_INPUT  = $INPUT_ROOT"/"${2}_${3}  
set RDIR        = "CMSSW_"${4} # Your CMSSW version
set GTAG        = ${5}"::All"  # Your global TAG 
set CMSSW_PROJECT_SRC = $WORKDIR$RDIR"/src"
set STEP              = "Extractors/PatExtractor_2"

#
# STEP 1: running jobs control
#
# There can't be more than $njoblimit running jobs ib batch
# We need to do this in order to avoid CASTOR overload
# we are rejected if making too many rfio requests 

cd $CMSSW_PROJECT_SRC

setenv SCRAM_ARCH slc5_amd64_gcc434

eval `scramv1 runtime -csh`

cd $CMSSW_PROJECT_SRC/$STEP/batch

set n_running     = `bjobs | wc -l` # Number of running jobs
@ njob            = $n_running 
@ njoblimit       = 30              # Max number of running jobs

if ($njob >= $njoblimit) then
    echo "Too many jobs are already running, you have to wait a bit..."
    #exit
endif


#
# STEP 2: look for data and launch job, if necessary
#


@ ndatafileslimit = 20 # If more than 20 datafiles, the global run is sliced apart              
@ ndatafilestreat = 0 

# How many files on tape
set nfiles   = `nsls $DATA_INPUT | wc -l`


# We launch the pat_extraction jobs

if ($nfiles > 0) then
    while ($nfiles > $ndatafilestreat)
	echo 'Sending job dealing with '$ndatafileslimit' files starting from run '$ndatafilestreat
	rm TMP_FILES/pat_extr_${2}_${3}_$ndatafilestreat.sh
	set rootdir = "$DATA_INPUT"	
	echo "#\!/bin/bash" > TMP_FILES/pat_extr_${2}_${3}_$ndatafilestreat.sh
	echo "source $CMSSW_PROJECT_SRC/$STEP/batch/pat_extractor.sh $DATA_INPUT $DATA_OUTPUT $ndatafilestreat $ndatafileslimit $CMSSW_PROJECT_SRC $GTAG ${1}" >> TMP_FILES/pat_extr_${2}_${3}_$ndatafilestreat.sh
	chmod 755 TMP_FILES/pat_extr_${2}_${3}_$ndatafilestreat.sh
	
	if (${6} == "DO_BATCH") then
	    bsub -q 1nd -e /dev/null -o /tmp/${LOGNAME}_out.txt $CMSSW_PROJECT_SRC/$STEP/batch/TMP_FILES/pat_extr_${2}_${3}_$ndatafilestreat.sh
	endif

	@ ndatafilestreat = $ndatafilestreat + $ndatafileslimit
   end
endif
	

