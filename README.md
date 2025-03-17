# NanoAnaCRAB
NanoAna TTreeReader based code to skim or analyze tree using CRAB facility

*Author: Arnab Laha, Prachurjya Hazarika*

This repo consists of TTreeReader based codes to process NanoAOD dataset using CRAB facility.

*Script for analysis*

1) C++ script : VLLAna.C
2) header file: VLLAna.h 
3) parameter settings: ana_crab.C
4) execution macro: runana_crab.C

So, the flow at the execution level is runana_crab.C --> ana_crab.C--> VLLAna.C (VLLAna.h). Check that VLLAna.C is compiled correctly to reduce your stress and coffee count!

*Scripts for CRAB settings:*

1) crab_config.py
2) crab_script.sh (this file will be executed in each worker node)
3) shell_instructions.sh (where we do some parameter settings, feed to runana_crab.C, and execute the runana to launch the job)
4) FrameworkJobReport.xml: CRAB will know the jobs are finished, start collecting outputs from the current session

## Local test
Test locally your code before submitting to CRAB.

$```chmod 777 shell_instructions.sh```

$```./shell_instructions.sh```

This should take the input file from the PSet.py file, and executing the code (VLLAna). Don't forget to create a grid proxy!

## CRAB jobs
change the dataset and full DAS string as needed. Then run,

$ ```crab submit -c crab_config.py```

It is always a good practice to dry run (simulate what will happen in worker nodes after job submission) to detect any CRAB-related dependency.

$ ```crab submit -c crab_config.py --dryrun```


