#!/bin/bash

echo "Listing files in job directory..."
ls -lh


# Access input and output file parameters from PSet.py
input_file=$(python3 -c "import PSet; print(PSet.process.source.fileNames[0])")  # Input file

output_file="skimFile.root"

#Access custom parameters from PSet.py
data_flag="0"
year="2018"
era="NA"
dataset="mu"
samplename="DY"

#Build the ROOT command
root_command="runana_crab.C(\"$input_file\",\"$output_file\",\"$data_flag\",\"$year\",\"$dataset\",\"$era\",\"$samplename\")"

#Echo the command that will be executed
echo "Running ROOT command: $root_command"

# Now pass all these parameters to the ROOT macro
root -q -b -l "$root_command"
