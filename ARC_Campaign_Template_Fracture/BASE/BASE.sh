#!/bin/bash

#SBATCH -t 24:00:00                 # Job time limit
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=64                # Number of tasks (total CPUs across all nodes)
#SBATCH --cpus-per-task=1           # CPUs per task
#SBATCH -p normal_q                 # Queue (partition) to submit to
#SBATCH -A ray_time                   # Account name
#SBATCH --job-name=BASE   # Job name
#SBATCH --output=current_working_directory/%x_%j.out
#SBATCH --error=current_working_directory/%x_%j.err


module reset
module load OpenMPI/4.1.5-GCC-12.3.0


# !!!!! Script is run in the parent SCRATCH FOLDER


# --------------------------------------------------------------------------
# Create log files and functions

# Define the path to the simulation campaign's logging file
OG_SIM_DIR="current_working_directory"
OG_PARENT_DIR=$(dirname "$OG_SIM_DIR")
CAMPAIGN_LOG_DIR="$OG_PARENT_DIR/Automated_Campaign_Output"

SUBMISSION_SCRIPT_LOG="$CAMPAIGN_LOG_DIR/submission_script_log.txt"
ERROR_LOG="$CAMPAIGN_LOG_DIR/error_log.txt"

touch $ERROR_LOG


# Function to format the entry of logging data

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $@" >> "$SUBMISSION_SCRIPT_LOG"
}

log_Error() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $@" >> "$ERROR_LOG"
}

# --------------------------------------------------------------------------


# the location of this ls-dyna job submission script
WORKING_DIR="current_working_directory"

OUTPUT_DIR="${WORKING_DIR/home/projects/energy_tpl}"

# Define the scratch space location
SCRATCH_DIR="$(pwd)/BASE_SlurmID_$SLURM_JOB_ID"

# Create the scratch space directory if it doesn't exist
mkdir -p $SCRATCH_DIR



# Copy files to be run to scratch space
cp "$WORKING_DIR/../projectile_geometry.k" "$SCRATCH_DIR" || { log_Error "projectile_geometry.k file copy to scratch directory failed"; exit 1; }
cp "$WORKING_DIR"/* "$SCRATCH_DIR" || { log_Error "Copy to scratch directory failed"; exit 1; }


# Change to scratch space directory
cd $SCRATCH_DIR || { log_Error "Could not change to scratch directory: $SCRATCH_DIR"; exit 1; }

#Create file containing the path to corresponding simulation folder in the home directory (i.e., the location of this parent bash script)
echo $OUTPUT_DIR > OutputSimFolder.txt
echo $WORKING_DIR > WorkingDirFolder.txt

# Define environment variables for the license server
export LSTC_LICENSE=network
export LSTC_LICENSE_SERVER=lstc.software.vt.edu

# Run simulation
mpirun --mca pml ob1 --mca btl tcp,vader,self /home/harryjang/00_LS_DYNA_Exe/ls-dyna_mpp_d_R15_0_2_x64_centos79_ifort190_avx2_openmpi405 i=BASE.k d=nodump NCPU=64 memory=2000m memory2=1000m > BASE.txt

# Check for errors in LS-DYNA execution
if [ $? -ne 0 ]; then
  log_Error "TINKERCLIFFS: Run error in scratch directory: $SCRATCH_DIR."
  exit 1
fi

# Create the output directory in the projects folder
mkdir -p $OUTPUT_DIR


# Move files we'd like to keep to the projects directory
AnalyzeFolder="./Analyze"

{
    cp *binout* "$OUTPUT_DIR"
    cp *.txt "$OUTPUT_DIR"
    cp *d3plot* "$OUTPUT_DIR"
    cp -r "$AnalyzeFolder" "$OUTPUT_DIR"
} 

if [[ $? -ne 0 ]]; then
    log_Error "Error: One or more copy operations to the output directory failed: $SCRATCH_DIR"
    exit 1
fi


SimFolder=$(pwd)
cd ..

# Remove the scratch directory folder
sleep 30
rm -r "$SimFolder"


# Check for successful removal of simulation file
if [ $? -eq 0 ]; then
    # Reduce 1 from the simulation left to run count

    # Define the file path
    FILE="SimsLeft.txt"

    # Check if the file exists
    if [[ -f "$FILE" ]]; then
        # Read the number from the file
        num=$(cat "$FILE")

        # Check if the number is an integer
        if [[ "$num" =~ ^-?[0-9]+$ ]]; then
            # Subtract 1
            ((num--))

            # Update the file with the new value
            echo "$num" > "$FILE"
        else
            log_Error "Error: $FILE does not contain a valid integer: $SCRATCH_DIR"
            exit 1
        fi
    else
        log_Error "Error: $FILE not found: $SCRATCH_DIR"
        exit 1
    fi
else

    log_Error "Error: Could not successfully remove completed simulation folder: $SCRATCH_DIR"
    sleep 120
    rm -r "$SimFolder"
    if [ $? -eq 0 ]; then
        log_Error "Error: Failed AGAIN to remove completed simulation folder: $SCRATCH_DIR"
    fi
    
fi


exit 0
