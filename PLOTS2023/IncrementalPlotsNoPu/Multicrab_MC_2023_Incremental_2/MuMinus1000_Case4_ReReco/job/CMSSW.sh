#!/bin/sh

# put STDERR to STDOUT 
exec 2>&1

echo "This script was generated by crab (version 2.10.5_patch1)."
#
# HEAD
#
#
echo "Running $0 with $# positional parameters: $*"

getRandSeed() {
     den=(0 1 2 3 4 5 6 7 8 9 A B C D E F G H I J K L M N O P Q R S T U V W X Y Z a b c d e f g h i j k l m n o p q r s t u v w x y z)
    nd=${#den[*]}
    randj=${den[$RANDOM % $nd]}${den[$RANDOM % $nd]}${den[$RANDOM % $nd]}
    echo $randj
}

dumpStatus() {
    echo ">>> info for dashboard:"
    echo "***** Cat $1 *****"
    cat $1
    echo "***** End Cat jobreport *****"
    chmod a+x $RUNTIME_AREA/report.py

    $RUNTIME_AREA/report.py $(cat $1)
    rm -f $1
    echo "MonitorJobID=`echo $MonitorJobID`" > $1
    echo "MonitorID=`echo $MonitorID`" >> $1
}


### REMOVE THE WORKING_DIR IN OSG SITES ###
remove_working_dir() {
    cd $RUNTIME_AREA
    echo ">>> working dir = $WORKING_DIR"
    echo ">>> current directory (RUNTIME_AREA): $RUNTIME_AREA"
    echo ">>> Remove working directory: $WORKING_DIR"
    /bin/rm -rf $WORKING_DIR
    if [ -d $WORKING_DIR ] ;then
        echo "ERROR ==> OSG $WORKING_DIR could not be deleted on WN `hostname`"
        job_exit_code=10017
    fi
}

### DUMP ORIGINAL ENVIRONMENT BEFORE CMSSW CUSTOMIZATIOn
dumpEnv(){
echo export PATH=$PATH >> CacheEnv.sh
echo export LD_LIBRARY_PATH=$LD_LIBRARY_PATH >> CacheEnv.sh
echo export PYTHONPATH=$PYTHONPATH >> CacheEnv.sh
}

outOfBound()
{
echo "********** CRAB WRAPPER CMSSW.sh TERMINATED BY A KILL"
if  [ ! -f ${RUNTIME_AREA}/WATCHDOG-SAYS-EXCEEDED-RESOURCE ]
then
  echo "********** KILL WAS NOT ISSUED BY CRAB WATCHDOG. SET GENERIC EXIT CODE"
  job_exit_code=50669
else
  echo "********** KILL signal was issued by Crab Watchdog"
  exceededResource=`cat  ${RUNTIME_AREA}/WATCHDOG-SAYS-EXCEEDED-RESOURCE`
  echo "********** because usage of resource ${exceededResource} was excessive"
  case ${exceededResource} in 
    "RSS"        ) job_exit_code=50660 ;;
    "VSIZE"      ) job_exit_code=50661 ;;
    "DISK"       ) job_exit_code=50662 ;;
    "CPU TIME"   ) job_exit_code=50663 ;;
    "WALL TIME"  ) job_exit_code=50664 ;;
    * ) echo "watchdog kill reason not given, set generic exit code"; job_exit_code=50669 ;;
  esac
  echo "************** JOB EXIT CODE set to: ${job_exit_code}"
fi
func_exit
}

exitBySignal()
{
echo "********** ${externalSignal} detected at `date`  -  `date -u`"
outOfBound
}

detectSIGINT()
{
externalSignal='SIGINT'; exitBySignal
}
detectSIGUSR1()
{
externalSignal='SIGUSR1'; exitBySignal
}
detectSIGUSR2()
{
externalSignal='SIGUSR2'; exitBySignal
}
detectSIGXCPU()
{
externalSignal='SIGXCPU'; exitBySignal
}
detectSIGXFSZ()
{
externalSignal='SIGXFSZ'; exitBySignal
}
detectSIGTERM()
{
externalSignal='SIGTERM'; exitBySignal
}



#
# EXECUTE THIS FUNCTION BEFORE EXIT 
#

func_exit() { 
    echo "CrabWrapper EXIT_FUNC entered at `date`   -  `date -u`"
    if [ -f ${RUNTIME_AREA}/WATCHDOG-SAYS-EXCEEDED-RESOURCE ]; then
       echo "*** Watchdog kicked in, make sure we"
       echo "*** do not go on before Watchdog completes the cleanup"
       wait ${WatchdogPID}
    fi
    if [ ! -s $RUNTIME_AREA/fillCrabFjr.py ]; then 
        echo "WARNING: it is not possible to create crab_fjr.xml to final report" 
    else 
        python $RUNTIME_AREA/fillCrabFjr.py $RUNTIME_AREA/crab_fjr_$NJob.xml --errorcode $job_exit_code $executable_exit_status 
    fi
    cd $RUNTIME_AREA  
    echo "Checking output files:  "
    for file in $filesToCheck ; do
        if [ -e $file ]; then
            echo "will tar file $file in  $out_files"
        else
            echo "WARNING: output file $file not found!"
        fi
    done
    if [ $middleware == OSG ]; then
        if [ $WORKING_DIR ]; then
            remove_working_dir
        fi
        symlinks -d .
    fi
    TIME_WRAP_END=`date +%s`
    let "TIME_WRAP = TIME_WRAP_END - TIME_WRAP_INI" 

    let "MIN_JOB_DURATION = 60*10" 
    let "PADDING_DURATION = MIN_JOB_DURATION - TIME_WRAP" 
    if [ $PADDING_DURATION -gt 0 ]; then 
        echo ">>> padding time: Sleeping the wrapper for $PADDING_DURATION seconds"
        sleep $PADDING_DURATION
        TIME_WRAP_END=`date +%s`
        let "TIME_WRAP = TIME_WRAP_END - TIME_WRAP_INI" 
    else 
        echo ">>> padding time: Wrapper lasting more than $MIN_JOB_DURATION seconds. No sleep required."
    fi

    echo "STOPPING WATCHDOG. CrabWatchdog PID is ${WatchdogPID}"
    kill $WatchdogPID
    echo "********** LAST 50 LINES OF WATCHDOG LOG"
    tail -50 Watchdog_$NJob.log
    echo "********** WATCHDOG LOG ENDED"
    echo "ZIPPING WATCHDOG LOG"
    gzip Watchdog_$NJob.log
    if [ ! -s $RUNTIME_AREA/fillCrabFjr.py ]; then 
        echo "WARNING: it is not possible to create crab_fjr.xml to final report" 
    else 
        set -- $CPU_INFOS 
        echo "CrabUserCpuTime=$1" >>  $RUNTIME_AREA/$repo 
        echo "CrabSysCpuTime=$2" >>  $RUNTIME_AREA/$repo 
        echo "CrabCpuPercentage=$3" >>  $RUNTIME_AREA/$repo 
        python $RUNTIME_AREA/fillCrabFjr.py $RUNTIME_AREA/crab_fjr_$NJob.xml --timing $TIME_WRAP $TIME_EXE $TIME_STAGEOUT \"$CPU_INFOS\" 
        echo "CrabWrapperTime=$TIME_WRAP" >> $RUNTIME_AREA/$repo 
        if [ $TIME_STAGEOUT -lt 0 ]; then 
            export TIME_STAGEOUT=NULL 
        fi
        echo "CrabStageoutTime=$TIME_STAGEOUT" >> $RUNTIME_AREA/$repo 
    fi
    echo "Disk space used:"
    echo "du -sh $RUNTIME_AREA"
    du -sh $RUNTIME_AREA 

    final_list=$filesToCheck
#Check for stdout/err in new location as of condor 7.7
    stdo=`ls -l /proc/$$/fd/1|awk '{print $NF}'`
    stde=`ls -l /proc/$$/fd/2|awk '{print $NF}'`
    if ! [ `basename $stdo` ==  CMSSW_${NJob}.stdout ]; then
      echo "Found unusual stdout, rename for OSB"
      cp -pfv $stdo CMSSW_${NJob}.stdout
    fi
    if ! [ `basename $stde` ==  CMSSW_${NJob}.stderr ]; then
      echo "Found unusual stderr, rename for OSB"
      cp -pfv $stde CMSSW_${NJob}.stderr
    fi
        tar zcvf ${out_files}.tgz  ${final_list}
    osb_size=`ls -gG ${out_files}.tgz | awk '{ print $3 }'`
    size=`expr $osb_size`
    echo "Total Output dimension: $size"
        rm ${out_files}.tgz
    limit=50000000 
    echo "WARNING: output files size limit is set to: $limit"
    if [ "$limit" -lt "$size" ]; then
        job_exit_code=70000
        echo "Output Sanbox too big. Produced output is lost "
        python $RUNTIME_AREA/fillCrabFjr.py $RUNTIME_AREA/crab_fjr_$NJob.xml --errorcode $job_exit_code $executable_exit_status 
        final_list="CMSSW_${NJob}.stdout CMSSW_${NJob}.stderr crab_fjr_${NJob}.xml"
    else
        echo "Total Output dimension $size is fine."
    fi
    echo "JOB_EXIT_STATUS = $job_exit_code"
    echo "JobExitCode=$job_exit_code" >> $RUNTIME_AREA/$repo
    dumpStatus $RUNTIME_AREA/$repo
    echo "EXITING at `date`  -  `date -u`"
    if [ -s _condor_stdout ]; then
      cp -pfv _condor_stdout CMSSW_${NJob}.stdout
      cp -pfv _condor_stderr CMSSW_${NJob}.stderr
    fi
    tar zcvf ${out_files}.tgz  ${final_list}
    exit $job_exit_code
}



RUNTIME_AREA=`pwd`
export RUNTIME_AREA

echo "Today is `date` - `date -u`"
echo "Job submitted on host `hostname`"
uname -a
echo ">>> current directory (RUNTIME_AREA): `pwd`"
#echo ">>> current directory content:"
#ls -Al
echo ">>> current user: `id`"
chmod 0700 -R .
echo ">>> directory permission set to 0700"
# proxy file needs special setting
proxyFile=`voms-proxy-info -path`
chmod 0600 ${proxyFile}

umask 077
echo ">>> umask set to: " `umask -S`

chmod 0704 /proc/$$/fd/1
chmod 0704 /proc/$$/fd/2
echo ">>> stdout stderr permission set to 0704"

echo ">>> voms proxy information:"
voms-proxy-info -all

repo=jobreport.txt
echo "WNHostName=`hostname`" | tee -a $RUNTIME_AREA/$repo


#Written by cms_cmssw::wsUntarSoftware
echo ">>> tar --no-same-permissions -xf $RUNTIME_AREA/default.tgz :" 
tar --no-same-permissions -xf $RUNTIME_AREA/default.tgz
untar_status=$? 
if [ $untar_status -ne 0 ]; then 
   echo "ERROR ==> Untarring .tgz file failed"
   job_exit_code=$untar_status
   func_exit
else 
   echo "Successful untar" 
fi 

echo ">>> Include $RUNTIME_AREA in PYTHONPATH:"
if [ -z "$PYTHONPATH" ]; then
   export PYTHONPATH=$RUNTIME_AREA/
else
   export PYTHONPATH=$RUNTIME_AREA/:${PYTHONPATH}
echo "PYTHONPATH=$PYTHONPATH"
fi


#
# SETUP ENVIRONMENT
#

export TIME_WRAP_INI=`date +%s` 
export TIME_STAGEOUT=-2 

# remoteglidein specific stuff
# strip arguments
echo "strip arguments"
args=("$@")
nargs=$#
shift $nargs
# job number (first parameter for job wrapper)
NJob=${args[0]}; export NJob
NResub=${args[1]}; export NResub
NRand=`getRandSeed`; export NRand
OutUniqueID=_$NRand
OutUniqueID=_$NResub$OutUniqueID
OutUniqueID=$NJob$OutUniqueID; export OutUniqueID
CRAB_UNIQUE_JOB_ID=calabria_MuMinus1000_Case4_ReReco_u805tn_${OutUniqueID}; export CRAB_UNIQUE_JOB_ID
echo env var CRAB_UNIQUE_JOB_ID set to: ${CRAB_UNIQUE_JOB_ID}
out_files=out_files_${NJob}; export out_files
echo $out_files
echo ">>> list of expected files on output sandbox"
echo "output files: crab_fjr_$NJob.xml CMSSW_$NJob.stdout CMSSW_$NJob.stderr Watchdog_$NJob.log.gz"
filesToCheck="crab_fjr_$NJob.xml CMSSW_$NJob.stdout CMSSW_$NJob.stderr Watchdog_$NJob.log.gz"
export filesToCheck
if [ $Glidein_MonitorID ]; then 
   MonitorJobID=${NJob}_${Glidein_MonitorID}
   SyncGridJobId=${Glidein_MonitorID}
else 
   MonitorJobID=${NJob}_https://cmssusy.ba.infn.it/c1b504cd81b2e4507c0b7c332a83ccd0a53ffe0e/${NJob}
   SyncGridJobId=https://cmssusy.ba.infn.it/c1b504cd81b2e4507c0b7c332a83ccd0a53ffe0e/${NJob}
fi
MonitorID=calabria_MuMinus1000_Case4_ReReco_u805tn
echo "MonitorJobID=$MonitorJobID" >> $RUNTIME_AREA/$repo 
echo "SyncGridJobId=$SyncGridJobId" >> $RUNTIME_AREA/$repo 
echo "MonitorID=$MonitorID" >> $RUNTIME_AREA/$repo
echo ">>> GridFlavour discovery: " 
if [ $OSG_GRID ]; then 
    middleware=OSG 
    echo "source OSG GRID setup script" 
    source $OSG_GRID/setup.sh 
elif [ $NORDUGRID_CE ]; then 
    middleware=ARC 
elif [ $VO_CMS_SW_DIR ]; then 
    middleware=LCG 
else 
    echo "ERROR ==> GridFlavour not identified" 
    job_exit_code=10030 
    func_exit 
fi 
echo "GridFlavour=$middleware" 
echo "GridFlavour=$middleware" >> $RUNTIME_AREA/$repo 

echo ">>> SyncSite discovery: " 
if [ $GLIDEIN_CMSSite ]; then 
    SyncSite=$GLIDEIN_CMSSite 
    echo "SyncSite=$SyncSite" 
    echo "SyncSite=$SyncSite" >> $RUNTIME_AREA/$repo ;
else
    echo "not reporting SyncSite"
fi

echo ">>> SyncCE discovery: " 
if [ $OSG_JOB_CONTACT ]; then 
    echo "getting SyncCE from OSG_JOB_CONTACT" 
    SyncCE="$OSG_JOB_CONTACT"; 
elif [ $NORDUGRID_CE ]; then 
    echo "getting SyncCE from NORDUGRID_CE" 
    SyncCE="${NORDUGRID_CE}:2811/nordugrid-GE-${QUEUE:-queue}"
 elif [ $CE_ID ]; then 
    echo "getting SyncCE from CE_ID" 
    SyncCE="${CE_ID}" 
elif [ "$GLIDEIN_Gatekeeper" ]; then 
    echo "getting SyncCE from GLIDEIN_Gaekeeper" 
    GKtmp="`echo $GLIDEIN_Gatekeeper | sed -e s,http's\?'://,,`"
    SyncCE="`echo $GKtmp | cut -d: -f1`" 
else 
    echo "getting SyncCE glite-brokerinfo" 
    SyncCE="`glite-brokerinfo getCE`" 
fi 
echo "SyncCE=$SyncCE" 
echo "SyncCE=$SyncCE" >> $RUNTIME_AREA/$repo ;
dumpStatus $RUNTIME_AREA/$repo 

export VO=cms

dumpEnv


#Written by cms_cmssw::wsSetupEnvironment
echo ">>> setup environment"
echo "set SCRAM ARCH to slc5_amd64_gcc472"
export SCRAM_ARCH=slc5_amd64_gcc472
echo "SCRAM_ARCH = $SCRAM_ARCH"
if [ $middleware == LCG ] || [ $middleware == CAF ] || [ $middleware == LSF ]; then 

#Written by cms_cmssw::wsSetupCMSLCGEnvironment_
    echo ">>> setup CMS LCG environment:"
    echo "set SCRAM ARCH and BUILD_ARCH to slc5_amd64_gcc472 ###"
    export SCRAM_ARCH=slc5_amd64_gcc472
    export BUILD_ARCH=slc5_amd64_gcc472
    if [ ! $VO_CMS_SW_DIR ] ;then
        echo "ERROR ==> CMS software dir not found on WN `hostname`"
        job_exit_code=10031
        func_exit
    else
        echo "Sourcing environment... "
        if [ ! -s $VO_CMS_SW_DIR/cmsset_default.sh ] ;then
            echo "ERROR ==> cmsset_default.sh file not found into dir $VO_CMS_SW_DIR"
            job_exit_code=10020
            func_exit
        fi
        echo "sourcing $VO_CMS_SW_DIR/cmsset_default.sh"
        source $VO_CMS_SW_DIR/cmsset_default.sh
        result=$?
        if [ $result -ne 0 ]; then
            echo "ERROR ==> problem sourcing $VO_CMS_SW_DIR/cmsset_default.sh"
            job_exit_code=10032
            func_exit
        fi
    fi
    
    echo "==> setup cms environment ok"
elif [ $middleware == OSG ]; then
    WORKING_DIR=`/bin/mktemp  -d $OSG_WN_TMP/cms_XXXXXXXXXXXX`
    if [ ! $? == 0 ] ;then
        echo "ERROR ==> OSG $WORKING_DIR could not be created on WN `hostname`"
        job_exit_code=10016
        func_exit
    fi
    echo ">>> Created working directory: $WORKING_DIR"

    echo "Change to working directory: $WORKING_DIR"
    cd $WORKING_DIR
    echo ">>> current directory (WORKING_DIR): $WORKING_DIR"

#Written by cms_cmssw::wsSetupCMSOSGEnvironment_
    echo ">>> setup CMS OSG environment:"
    echo "set SCRAM ARCH to slc5_amd64_gcc472"
    export SCRAM_ARCH=slc5_amd64_gcc472
    echo "SCRAM_ARCH = $SCRAM_ARCH"
    echo "OSG_APP is $OSG_APP"
    if [ -f $OSG_APP/cmssoft/cms/cmsset_default.sh ] ;then
        cmsSetupFile=$OSG_APP/cmssoft/cms/cmsset_default.sh
    elif [ -f $CVMFS/cms.cern.ch/cmsset_default.sh ] ; then 
        cmsSetupFile=$CVMFS/cms.cern.ch/cmsset_default.sh
    elif [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ] ; then 
        cmsSetupFile=/cvmfs/cms.cern.ch/cmsset_default.sh
    else
        echo "CVMSF = $CVMFS"
        echo "/cvmfs/ is"
        echo "ls /"
        ls /
        echo "ls /cvmfs"
        ls /cvmfs
        echo "ls /cvmfs/cms.cern.ch"
        ls /cvmfs/cms.cern.ch
        ls /cvmfs/cms.cern.ch/cmsset*
        ls /cvmfs/cms.cern.ch/cmsset_default.sh
        echo "ERROR ==> cmsset_default.sh file not found"
        job_exit_code=10020
        func_exit
    fi

    echo "sourcing $cmsSetupFile ..."
    source $cmsSetupFile
    result=$?
    if [ $result -ne 0 ]; then
       echo "ERROR ==> problem sourcing $cmsSetupFile"
       job_exit_code=10032
       func_exit
    else
      echo "==> setup cms environment ok"
      echo "SCRAM_ARCH = $SCRAM_ARCH"
    fi
elif [ $middleware == SGE ]; then

#Written by cms_cmssw::wsSetupCMSLCGEnvironment_
    echo ">>> setup CMS LCG environment:"
    echo "set SCRAM ARCH and BUILD_ARCH to slc5_amd64_gcc472 ###"
    export SCRAM_ARCH=slc5_amd64_gcc472
    export BUILD_ARCH=slc5_amd64_gcc472
    if [ ! $VO_CMS_SW_DIR ] ;then
        echo "ERROR ==> CMS software dir not found on WN `hostname`"
        job_exit_code=10031
        func_exit
    else
        echo "Sourcing environment... "
        if [ ! -s $VO_CMS_SW_DIR/cmsset_default.sh ] ;then
            echo "ERROR ==> cmsset_default.sh file not found into dir $VO_CMS_SW_DIR"
            job_exit_code=10020
            func_exit
        fi
        echo "sourcing $VO_CMS_SW_DIR/cmsset_default.sh"
        source $VO_CMS_SW_DIR/cmsset_default.sh
        result=$?
        if [ $result -ne 0 ]; then
            echo "ERROR ==> problem sourcing $VO_CMS_SW_DIR/cmsset_default.sh"
            job_exit_code=10032
            func_exit
        fi
    fi
    
    echo "==> setup cms environment ok"
elif [ $middleware == ARC ]; then

#Written by cms_cmssw::wsSetupCMSLCGEnvironment_
    echo ">>> setup CMS LCG environment:"
    echo "set SCRAM ARCH and BUILD_ARCH to slc5_amd64_gcc472 ###"
    export SCRAM_ARCH=slc5_amd64_gcc472
    export BUILD_ARCH=slc5_amd64_gcc472
    if [ ! $VO_CMS_SW_DIR ] ;then
        echo "ERROR ==> CMS software dir not found on WN `hostname`"
        job_exit_code=10031
        func_exit
    else
        echo "Sourcing environment... "
        if [ ! -s $VO_CMS_SW_DIR/cmsset_default.sh ] ;then
            echo "ERROR ==> cmsset_default.sh file not found into dir $VO_CMS_SW_DIR"
            job_exit_code=10020
            func_exit
        fi
        echo "sourcing $VO_CMS_SW_DIR/cmsset_default.sh"
        source $VO_CMS_SW_DIR/cmsset_default.sh
        result=$?
        if [ $result -ne 0 ]; then
            echo "ERROR ==> problem sourcing $VO_CMS_SW_DIR/cmsset_default.sh"
            job_exit_code=10032
            func_exit
        fi
    fi
    
    echo "==> setup cms environment ok"
elif [ $middleware == PBS ] || [ $middleware == PBSV2 ]; then

#Written by cms_cmssw::wsSetupCMSLCGEnvironment_
    echo ">>> setup CMS LCG environment:"
    echo "set SCRAM ARCH and BUILD_ARCH to slc5_amd64_gcc472 ###"
    export SCRAM_ARCH=slc5_amd64_gcc472
    export BUILD_ARCH=slc5_amd64_gcc472
    if [ ! $VO_CMS_SW_DIR ] ;then
        echo "ERROR ==> CMS software dir not found on WN `hostname`"
        job_exit_code=10031
        func_exit
    else
        echo "Sourcing environment... "
        if [ ! -s $VO_CMS_SW_DIR/cmsset_default.sh ] ;then
            echo "ERROR ==> cmsset_default.sh file not found into dir $VO_CMS_SW_DIR"
            job_exit_code=10020
            func_exit
        fi
        echo "sourcing $VO_CMS_SW_DIR/cmsset_default.sh"
        source $VO_CMS_SW_DIR/cmsset_default.sh
        result=$?
        if [ $result -ne 0 ]; then
            echo "ERROR ==> problem sourcing $VO_CMS_SW_DIR/cmsset_default.sh"
            job_exit_code=10032
            func_exit
        fi
    fi
    
    echo "==> setup cms environment ok"
fi


echo ">>> specific cmssw setup environment:"
echo "CMSSW_VERSION =  CMSSW_6_2_0_SLHC12"
scram project CMSSW CMSSW_6_2_0_SLHC12
status=$?
if [ $status != 0 ] ; then
    echo "ERROR ==> CMSSW CMSSW_6_2_0_SLHC12 not found on `hostname`" 
    job_exit_code=10034
    func_exit
fi 
cd CMSSW_6_2_0_SLHC12
SOFTWARE_DIR=`pwd`; export SOFTWARE_DIR
echo ">>> current directory (SOFTWARE_DIR): $SOFTWARE_DIR" 
eval `scram runtime -sh | grep -v SCRAMRT_LSB_JOBNAME`
if [ $? != 0 ] ; then
    echo "ERROR ==> Problem with the command: "
    echo "eval \`scram runtime -sh | grep -v SCRAMRT_LSB_JOBNAME \` at `hostname`"
    job_exit_code=10034
    func_exit
fi 

## number of arguments (first argument always jobnumber, the second is the resubmission number)

if [ $nargs -lt 2 ]
then
    echo 'ERROR ==> Too few arguments' +$nargs+ 
    job_exit_code=50113
    func_exit
fi


DatasetPath=/calabria_SingleMuPt1000_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case4/calabria-calabria_SingleMuPt1000_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case4-fa5ffd20b461790f8466430606d05ce8/USER
PrimaryDataset=calabria_SingleMuPt1000_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case4
DataTier=calabria-calabria_SingleMuPt1000_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case4-fa5ffd20b461790f8466430606d05ce8
ApplicationFamily=cmsRun

cp  $RUNTIME_AREA/CMSSW.py .
cp  $RUNTIME_AREA/CMSSW.py.pkl .
PreserveSeeds=; export PreserveSeeds
IncrementSeeds=; export IncrementSeeds
echo "PreserveSeeds: <$PreserveSeeds>"
echo "IncrementSeeds:<$IncrementSeeds>"
mv -f CMSSW.py pset.py

dumpStatus $RUNTIME_AREA/$repo

#
# END OF SETUP ENVIRONMENT
#

#
# fork watchdog
#
${RUNTIME_AREA}/crabWatchdog.sh &
export WatchdogPID=$!
echo "crabWatchdog started as process $WatchdogPID"
trap detectSIGINT  SIGINT
trap detectSIGUSR1 SIGUSR1
trap detectSIGUSR2 SIGUSR2
trap detectSIGXCPU SIGXCPU
trap detectSIGXFSZ SIGXFSZ
trap detectSIGTERM SIGTERM

#
# PREPARE AND RUN EXECUTABLE
#


#Written by cms_cmssw::wsBuildExe
echo ">>> moving CMSSW software directories in `pwd`" 
rm -rf lib/ module/ 
mv $RUNTIME_AREA/lib/ . 
mv $RUNTIME_AREA/module/ . 
rm -rf src/ 
mv $RUNTIME_AREA/src/ . 
echo ">>> Include $RUNTIME_AREA in PYTHONPATH:"
if [ -z "$PYTHONPATH" ]; then
   export PYTHONPATH=$RUNTIME_AREA/
else
   export PYTHONPATH=$RUNTIME_AREA/:${PYTHONPATH}
echo "PYTHONPATH=$PYTHONPATH"
fi


edmConfigHash pset.py 
PSETHASH=`edmConfigHash pset.py` 
echo "PSETHASH = $PSETHASH" 
if [ -z "$PSETHASH" ]; then 
   export PSETHASH=null
fi 

executable=cmsRun


#
# END OF PREPARE AND RUN EXECUTABLE
#

#
# COPY INPUT
#


#
# Rewrite cfg or cfgpy file
#

# Rewrite cfg for this job
echo  $RUNTIME_AREA/writeCfg.py  pset.py pset.py
python $RUNTIME_AREA/writeCfg.py  pset.py pset.py

        result=$?
        if [ $result -ne 0 ]; then
            echo "ERROR ==> problem re-writing config file"
            job_exit_code=10040
            func_exit
        fi

          
cat $RUNTIME_AREA/inputsReport.txt  

echo ">>> Executable $executable"
which $executable
res=$?
if [ $res -ne 0 ];then
  echo "ERROR ==> executable not found on WN `hostname`"
  job_exit_code=10035
  func_exit
else
  echo "ok executable found"
fi

echo "ExeStart=$executable" >>  $RUNTIME_AREA/$repo
dumpStatus $RUNTIME_AREA/$repo

echo ">>> $executable started at `date -u`"
start_exe_time=`date +%s`
CPU_INFOS=0 
/usr/bin/time -f "%U %S %P" -o cpu_timing.txt $executable  -j $RUNTIME_AREA/crab_fjr_$NJob.xml -p pset.py > executable.out 2>&1 
executable_exit_status=$?
CPU_INFOS=`tail -n 1 cpu_timing.txt`
stop_exe_time=`date +%s`
echo ">>> $executable ended at `date -u`"

chmod 0704 $RUNTIME_AREA/crab_fjr_*.xml
echo ">>> crab_fjr_*.xml  permission set to 0704"


#### dashboard add timestamp!
echo "ExeEnd=$executable" >> $RUNTIME_AREA/$repo
dumpStatus $RUNTIME_AREA/$repo

let "TIME_EXE = stop_exe_time - start_exe_time"
echo "TIME_EXE = $TIME_EXE sec"
echo "ExeTime=$TIME_EXE" >> $RUNTIME_AREA/$repo

#
# limit executable stdout size to 2K lines
#
exeOutLines=`wc -l executable.out | awk '{print $1'}`
echo ">>> $executable wrote $exeOutLines lines of stdout+stderr"
echo ">>> START OF printout of $exeOutLines lines from stdout+stderr"
if [ $exeOutLines -gt 3000 ]
then
  echo ">>> print out only first 1000 and last 2000 lines:"
  head -1000 executable.out; echo ""; echo ">>>[...BIG SNIP...]";echo "";tail -2000 executable.out
else
  cat executable.out
fi
echo ">>> END OF printout of $exeOutLines lines from stdout+stderr"

grep -q 'Fatal Exception' executable.out
fatal=$?
if [ ${fatal} == "0" ]
then
    echo ">>> ERROR: Fatal Exception from CMSSW:"
    awk '/Begin .* Exception/ { printing =1 } /End .* Exception/ { print $0; printing = 0 } printing { print $0 } ' executable.out 
fi

egrep -i 'segmentation.*fault|segmentation.*violation|segfault' executable.out
segFault=$?
if [ ${segFault} == "0" ]
then
  echo "ERROR ==> user executable segfaulted"
  job_exit_code=50800
  func_exit
fi

#
# if Watchdog killed executable, make sure we
# do not go on before Watchdog completes the cleanup
#
if [ -f ${RUNTIME_AREA}/WATCHDOG-SAYS-EXCEEDED-RESOURCE ]; then wait ${WatchdogPID}; fi


#Written by cms_cmssw::wsParseFJR
echo ">>> Parse FrameworkJobReport crab_fjr.xml"
if [ -s $RUNTIME_AREA/crab_fjr_$NJob.xml ]; then
    if [ -s $RUNTIME_AREA/parseCrabFjr.py ]; then
        cmd_out=`python $RUNTIME_AREA/parseCrabFjr.py --input $RUNTIME_AREA/crab_fjr_$NJob.xml --dashboard $MonitorID,$MonitorJobID `
        cmd_out_1=`python $RUNTIME_AREA/parseCrabFjr.py --input $RUNTIME_AREA/crab_fjr_$NJob.xml --popularity $MonitorID,$MonitorJobID,$RUNTIME_AREA/inputsReport.txt `
        echo "Result of parsing the FrameworkJobReport crab_fjr.xml: $cmd_out_1"
        executable_exit_status=`python $RUNTIME_AREA/parseCrabFjr.py --input $RUNTIME_AREA/crab_fjr_$NJob.xml --exitcode`
        if [ $executable_exit_status -eq 50115 ];then
            echo ">>> crab_fjr.xml contents: "
            cat $RUNTIME_AREA/crab_fjr_$NJob.xml
            echo "Wrong FrameworkJobReport --> does not contain useful info. ExitStatus: $executable_exit_status"
        elif [ $executable_exit_status -eq -999 ];then
            echo "ExitStatus from FrameworkJobReport not available. not available. Using exit code of executable from command line."
        else
            echo "Extracted ExitStatus from FrameworkJobReport parsing output: $executable_exit_status"
        fi
    else
        echo "CRAB python script to parse CRAB FrameworkJobReport crab_fjr.xml is not available, using exit code of executable from command line."
    fi
    if [ $executable_exit_status -eq 0 ];then
        echo ">>> Executable succeded  $executable_exit_status"
    fi
else
    echo "CRAB FrameworkJobReport crab_fjr.xml is not available, using exit code of executable from command line."
fi

if [ $executable_exit_status -ne 0 ];then
    echo ">>> Executable failed  $executable_exit_status"
    echo "ExeExitCode=$executable_exit_status" | tee -a $RUNTIME_AREA/$repo
    echo "EXECUTABLE_EXIT_STATUS = $executable_exit_status"
    job_exit_code=$executable_exit_status
    func_exit
fi

echo "ExeExitCode=$executable_exit_status" | tee -a $RUNTIME_AREA/$repo
echo "EXECUTABLE_EXIT_STATUS = $executable_exit_status"
job_exit_code=$executable_exit_status

#
# PROCESS THE PRODUCED RESULTS
#



#Written by cms_cmssw::wsRenameOutput
echo ">>> current directory $PWD" 
echo ">>>  (SOFTWARE_DIR): $SOFTWARE_DIR" 
echo ">>> (WORKING_DIR): $WORKING_DIR" 
echo ">>> current directory content:"
ls -Al


# check output file
if [ -e ./validationEDM.root ] ; then
    mv validationEDM.root validationEDM_$OutUniqueID.root
    ln -s `pwd`/validationEDM_$OutUniqueID.root $RUNTIME_AREA/validationEDM.root
else
    job_exit_code=60302
    echo "WARNING: Output file validationEDM.root not found"
fi
file_list="$SOFTWARE_DIR/validationEDM_$OutUniqueID.root"

echo ">>> current directory $PWD" 
echo ">>> (SOFTWARE_DIR): $SOFTWARE_DIR" 
echo ">>> (WORKING_DIR): $WORKING_DIR" 
echo ">>> current directory content:"
ls -Al

cd $RUNTIME_AREA
echo ">>> current directory (RUNTIME_AREA):  $RUNTIME_AREA"


#
# COPY OUTPUT FILE TO SE
#

export SE=storm-se-01.ba.infn.it
echo "SE = $SE"
export SE_PATH=/srm/managerv2?SFN=//cms/store/user/calabria/calabria_SingleMuPt1000_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case4/calabria_MuMinus1000_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_VALNOPU2_Case4/${PSETHASH}/
echo "SE_PATH = $SE_PATH"
export LFNBaseName=/store/user/calabria/calabria_SingleMuPt1000_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case4/calabria_MuMinus1000_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_VALNOPU2_Case4/${PSETHASH}/
echo "LFNBaseName = $LFNBaseName"
export USER=calabria
echo "USER = $USER"
export endpoint=srm://storm-se-01.ba.infn.it:8444/srm/managerv2?SFN=//cms/store/user/calabria/calabria_SingleMuPt1000_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_Case4/calabria_MuMinus1000_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_VALNOPU2_Case4/${PSETHASH}/
echo "endpoint = $endpoint"
echo ">>> Copy output files from WN = `hostname` to $SE_PATH :"
export TIME_STAGEOUT_INI=`date +%s` 
copy_exit_status=0
echo "python cmscp.py  --destination $endpoint --inputFileList $file_list --middleware $middleware --se_name $SE --for_lfn $LFNBaseName    "
python cmscp.py  --destination $endpoint --inputFileList $file_list --middleware $middleware --se_name $SE --for_lfn $LFNBaseName    
if [ -f $RUNTIME_AREA/resultCopyFile ] ;then
    cat $RUNTIME_AREA/resultCopyFile
    pwd
else
    echo "ERROR ==> $RUNTIME_AREA/resultCopyFile file not found. Problem during the stageout"
    echo "RUNTIME_AREA content: " 
    ls $RUNTIME_AREA 
    job_exit_code=60318
    func_exit 
fi
if [ -f ${RUNTIME_AREA}/cmscpReport.sh ] ;then
    echo "-------- cat ${RUNTIME_AREA}/cmscpReport.sh "
    cat ${RUNTIME_AREA}/cmscpReport.sh
    echo "-------- end of ${RUNTIME_AREA}/cmscpReport.sh "
    source ${RUNTIME_AREA}/cmscpReport.sh
    source_result=$? 
    if [ $source_result -ne 0 ]; then
        echo "problem with the source of cmscpReport.sh file"
        StageOutExitStatus=60307
    fi
else
    echo "cmscpReport.sh file not found"
    StageOutExitStatus=60307
fi
if [ $StageOutExitStatus -ne 0 ]; then
    echo "Problem copying file to $SE $SE_PATH"
    copy_exit_status=$StageOutExitStatus 
if [ -f .SEinteraction.log ] ;then
    echo "########## contents of SE interaction"
    cat .SEinteraction.log
    echo "#####################################"
else
    echo ".SEinteraction.log file not found"
fi
    job_exit_code=$StageOutExitStatus
fi
export TIME_STAGEOUT_END=`date +%s` 
let "TIME_STAGEOUT = TIME_STAGEOUT_END - TIME_STAGEOUT_INI" 

echo ">>> current dir: `pwd`"
echo ">>> current dir content:"
ls -Al


#Written by cms_cmssw::wsModifyReport
echo ">>> Modify Job Report:" 
chmod a+x $RUNTIME_AREA/ProdCommon/FwkJobRep/ModifyJobReport.py
echo "CMSSW_VERSION = $CMSSW_VERSION"

ProcessedDataset=calabria_MuMinus1000_GEN-SIM-DIGI-RECO_CMSSW_6_2_0_SLHC12_2023Scenario_VALNOPU2_Case4
echo "ProcessedDataset = $ProcessedDataset"
echo "$RUNTIME_AREA/ProdCommon/FwkJobRep/ModifyJobReport.py fjr $RUNTIME_AREA/crab_fjr_$NJob.xml json $RUNTIME_AREA/resultCopyFile n_job $OutUniqueID PrimaryDataset $PrimaryDataset  ApplicationFamily $ApplicationFamily ApplicationName $executable cmssw_version $CMSSW_VERSION psethash $PSETHASH UserProcessedDataset $USER-$ProcessedDataset-$PSETHASH"
$RUNTIME_AREA/ProdCommon/FwkJobRep/ModifyJobReport.py fjr $RUNTIME_AREA/crab_fjr_$NJob.xml json $RUNTIME_AREA/resultCopyFile n_job $OutUniqueID PrimaryDataset $PrimaryDataset  ApplicationFamily $ApplicationFamily ApplicationName $executable cmssw_version $CMSSW_VERSION psethash $PSETHASH UserProcessedDataset $USER-$ProcessedDataset-$PSETHASH
modifyReport_result=$?
if [ $modifyReport_result -ne 0 ]; then
    modifyReport_result=70500
    job_exit_code=$modifyReport_result
    echo "ModifyReportResult=$modifyReport_result" | tee -a $RUNTIME_AREA/$repo
    echo "WARNING: Problem with ModifyJobReport"
else
    mv NewFrameworkJobReport.xml $RUNTIME_AREA/crab_fjr_$NJob.xml
fi

chmod 0704 $RUNTIME_AREA/crab_fjr_*.xml
echo ">>> crab_fjr_*.xml  permission set to 0704"

#
# END OF PROCESS THE PRODUCED RESULTS
#


func_exit
