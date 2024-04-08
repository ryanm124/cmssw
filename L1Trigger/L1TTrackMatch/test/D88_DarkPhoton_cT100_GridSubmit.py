import os
import sys
from itertools import islice
submit = 'universe = vanilla\n' ## writing .sub file
submit += 'arguments = "$(argument)"\n'
submit += 'output = submit01.out\n'
submit += 'error = submit01.err\n'
submit += 'log = submit01.log\n'
submit += '+JobFlavour = "tomorrow"\n' 
submit += 'queue\n'
submitName = 'D88_DarkPhoton_mH125_mGammaD20_cT100/submit01.sub'
sub1 = open(submitName,'w')
sub1.write(submit+'\n')
sub1.close() ##finish writing .sub file
nfile = 2 ##number of sample file processed per each condor job
filename = 'DarkPhoton_mH125_mGammaD20_cT_100.txt' ## location of the txt file used to help locate samples
with open(filename,'r') as f:
    counter1 = 1 ## the nth job
    while True:
        lines = list(islice(f, nfile))
        if not lines:
            break
        counter2 = 1 ## the nth file in one job
        for line in lines:
            if counter2 == 1:
               input='root://cms-xrd-global.cern.ch//'+line.rstrip() ## Remove  'root:..'+  if you don't have grid certificate. 
            else:
               input =input+',root://cms-xrd-global.cern.ch//'+line.rstrip() ## Remove  ',root:..'+  if you don't have grid certificate.
            counter2+=1
        create = '#!/bin/bash\n' ##writng .sh file
        create += 'export CMSSW_PROJECT_SRC=/afs/cern.ch/user/r/rmccarth/private/dispVert/CMSSW_14_0_0_pre2/src\n' 
        create += 'cd $CMSSW_PROJECT_SRC\n' ## go to the src directly 
        create += 'eval `scramv1 runtime -sh`\n' ## this is effectively "cmsenv"
        create += 'export X509_USER_PROXY=/afs/cern.ch/user/r/rmccarth/x509up_u120948\n' ## exporting the proxy
        create += 'cd /afs/cern.ch/user/r/rmccarth/private/dispVert/CMSSW_14_0_0_pre2/src/L1Trigger/L1TTrackMatch/test\n'## go the directly where you have the producer and python config file
        create += 'cmsRun GridL1TrackObjectNtupleMaker_cfg.py inputFiles='+input+' outputFile=/afs/cern.ch/work/r/rmccarth/private/D88_DarkPhoton_mH125_mGammaD20_cT100_NNScores/'+str(160*nfile)+'events_D88_HybridDisp'+str(counter1)+'.root  maxEvents='+str(160*nfile)+'\n'
        createName = 'D88_DarkPhoton_mH125_mGammaD20_cT100/submit'+str(counter1)+'.sh'
        sub2 = open(createName,'w')
        sub2.write(create+'\n')
        sub2.close() ## finish writing .sh file
        counter1+=1
        os.system('chmod 755 '+createName) ## make .sh file executable
        os.system('condor_submit '+ submitName+' executable='+createName) ## submit the job using condor_submit command
