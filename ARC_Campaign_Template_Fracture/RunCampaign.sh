#! /bin/bash
#
#SBATCH -t 120:00:00
#SBATCH --nodes=1                   # Number of nodes
#SBATCH -p normal_q                 # Queue (partition) to submit to
#SBATCH -A ray_time                   # Account name
#SBATCH --exclusive

# --------------------------------
# load modules
# --------------------------------
module load python3


# !!!!! Be careful in setting ntasks and cpus-per-task.
# ntasks: at least the number of matlab scripts that can run at once
# cpus-per-task: the number of CPUs each matlab parallel worker pool asks for + 2 (1 matlab program runs the batch script and 1 serves as the base point from which the parallel workers are launched)


# Variable to track the total number of jobs that have been submitted to SLURM
total_jobs=0

MAX_JOBS=23

# Status of the cluster - 0 means nominal anything other than 0 indicates an issue
cluster_status=0

# Top-level directory from which to start the recursive search
TOP_LEVEL_DIR=$(pwd)

#Total number of simulations that will be run in the campaign
total_sims=$(find "$TOP_LEVEL_DIR" -type f -name 'MAIN*.sh' ! -name "*PostProcess.sh" | wc -l)
sims_left_to_run=$total_sims


# --------------------------------------------------------
# Create directories to store campaign status information
# --------------------------------------------------------


# Create folder to store data on campaign running / submission state.
CAMPAIGN_OUTPUT_DIR="$TOP_LEVEL_DIR/Automated_Campaign_Output"
mkdir -p "$CAMPAIGN_OUTPUT_DIR"

# File to keep track of all submitted jobs
SUBMITTED_JOBS_HISTORY="$CAMPAIGN_OUTPUT_DIR/submitted_jobs_history.txt"
touch "$SUBMITTED_JOBS_HISTORY"

# File to keep track of which jobs are actively running or queued
RUNNING_JOBS="$CAMPAIGN_OUTPUT_DIR/running_jobs.txt"
touch $RUNNING_JOBS

# File to log actions taken by submission script
SUBMISSION_SCRIPT_LOG="$CAMPAIGN_OUTPUT_DIR/submission_script_log.txt"
touch $SUBMISSION_SCRIPT_LOG
echo "==============================" >> $SUBMISSION_SCRIPT_LOG
echo "Begin automated simulation run" >> $SUBMISSION_SCRIPT_LOG
echo "" >> $SUBMISSION_SCRIPT_LOG


# File to log errors
ERROR_LOG="$CAMPAIGN_OUTPUT_DIR/error_log.txt"
touch $ERROR_LOG

# -----------------------------------------------------------------
# Create a unique scratch directory to run simulation campaign
# -----------------------------------------------------------------

# Define the root directory
root_dir="/scratch/$USER"

# Check if /scratch/$USER exists
if [ ! -d "$root_dir" ]; then
    xxx="000"
else
    # Find existing folders matching the pattern (three digits)
    existing_folders=($(ls -d "$root_dir"/[0-9][0-9][0-9] 2>/dev/null | grep -oE '[0-9]{3}' | sort -n))

    if [ ${#existing_folders[@]} -eq 0 ]; then
        xxx="000"
    else
        # Get the highest existing number and increment it
        highest=${existing_folders[-1]}
        highest=$((10#$highest))  # force base 10
        xxx=$(printf "%03d" $((highest + 1)))

    fi
fi

# Define the new folder path
SCRATCH_DIR="$root_dir/$xxx"


echo "SCRATCH LOCATION: $SCRATCH_DIR" >> $SUBMISSION_SCRIPT_LOG
echo "" >> $SUBMISSION_SCRIPT_LOG
echo "==============================" >> $SUBMISSION_SCRIPT_LOG

# --------------------------------------------------------
# Define functions
# --------------------------------------------------------


# Function to format the entry of logging data
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $@" >> "$SUBMISSION_SCRIPT_LOG"
}

log_Error(){
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $@" >> "$ERROR_LOG"
}


# Function to convert HH:MM:SS to seconds
convert_to_seconds() {
    local IFS=: # Setting internal field separator to colon for splitting input
    read -ra time_parts <<< "$1" # Splitting input into an array based on IFS
    local hours=${time_parts[0]}
    local minutes=${time_parts[1]}
    local seconds=${time_parts[2]}
    echo $((hours * 3600 + minutes * 60 + seconds))
}



#Evaluate the status of the cluster by seeing if it can run the "squeue" command in a timely manner
check_cluster_status() {
    timeout 5s squeue > /dev/null 2>&1
    # Set the global variable with the exit status of the timeout command
    cluster_status=$?
    
    if [ "$cluster_status" -ne 0 ]; then
        log " "
        log "Cluster appears to be struggling. Returned cluster_status = $cluster_status"
        log " "
    fi
}


# Function to check the total number of submitted jobs (running or queued) from the current run script

check_jobs() {
    check_cluster_status
    
    if [ "$cluster_status" -eq 0 ]; then
        total_jobs=0
        temp_file=$(mktemp)  # Create a temporary file

        if [ -f "$RUNNING_JOBS" ]; then
            while read -r job_id; do
                if squeue -j "$job_id" -h >/dev/null 2>&1; then
                    echo "$job_id" >> "$temp_file"  # Keep active jobs
                    ((total_jobs++))
                fi
            done < "$RUNNING_JOBS"

            mv "$temp_file" "$RUNNING_JOBS"  # Overwrite with active jobs
        fi
    fi
}


submit_jobs() {
    
    check_cluster_status
    
    if [ "$cluster_status" -ne 0 ]; then
        log " "
        log "Cluster appears to be struggling. Returned cluster_status = $cluster_status"
        log " "
        return 0
    fi
    
    log "--------------------- job submission ---------------------"
    
    # Recursively find all simulation script files
    while IFS= read -r sim_script; do
        sim_dir=$(dirname "$sim_script")
        
        # Check if the simulation has already been submitted by looking for its directory in the submitted_jobs file
        if ! grep -Fxq "$sim_dir" "$SUBMITTED_JOBS_HISTORY"; then
            
            # Submit the job and store the job ID
            JOB_ID=$(sbatch "$sim_script" | awk '{print $NF}')
                
            # Error checking on sbatch command
            if [ $? -ne 0 ]; then
                log_Error "Failed to submit job for script $sim_script"
                continue  # Skip this script and go to the next one
            else
                log "submitted job ID $JOB_ID to run the simulation $sim_dir"
            fi
            
            # Record the job as active for use in the check_jobs() function
            if [[ $JOB_ID =~ ^[0-9]+$ ]]; then
                echo "$JOB_ID" >> "$RUNNING_JOBS"
            fi
                
                
            # Record the directory as submitted
            echo "$sim_dir" >> "$SUBMITTED_JOBS_HISTORY"

            sleep 2
                
            # Get the number of submitted and running jobs
            check_jobs
                
                
            # If we've reached the maximum number of jobs, break from the loop
            if [ "$total_jobs" -ge "$MAX_JOBS" ]; then
                log "reached maximum job submission count with $total_jobs submitted jobs."
                log " "
                break
            fi
                
            if [ "$cluster_status" -ne 0 ]; then
                log " "
                log "Cluster appears to be struggling. Returned cluster_status = $cluster_status"
                log " "
                break
            fi
                
                
            
        fi
    done < <(find "$TOP_LEVEL_DIR" -type f -name "*MAIN*.sh" ! -name "*PostProcess.sh")
}


end_if_complete() {

    check_cluster_status
    
    if [ "$cluster_status" -ne 0 ]; then
        log " "
        log "Cluster appears to be struggling. Returned cluster_status = $cluster_status"
        log " "
        return 0
    fi
    
    # load Python module
    module load python3
  
    # Define the base directory where the simulations are stored
    local BASE_DIR=$SCRATCH_DIR

    # Change to the base directory
    cd "$BASE_DIR"

    # Find all directories starting with 'MAIN'
  
    log " "
    log "--------------------- Early termination evaluation ---------------------"
  
    for sim_dir in MAIN*_*; do
     
        cd "$sim_dir" || continue
    
    
        if [[ ! -f TERMINATION_STATUS.txt || $(<TERMINATION_STATUS.txt) == "n" ]]; then
    
            # Check if lock file exists
            if [[ -f PostprocessLock ]]; then
                log "Skipping $sim_dir as post-processing is actively running"
            else
                # Set up a process lock so that the post process script is not re-initiated while it is still running from previous call
                touch PostprocessLock
                chmod +x *PostProcess.sh
                ./*PostProcess.sh &
                sleep 2
            fi
        else
            log "Skipping $sim_dir as post-processing has completed"
        fi
        cd "$BASE_DIR"
    done
    
    log " "
}


# --------------------------------------------------------
# Main loop
# --------------------------------------------------------


# Create and move into scratch directory 

mkdir -p $SCRATCH_DIR
# Clean the directory to avoid contamination (this is a scratch space!)
rm -rf "$SCRATCH_DIR"/*
# copy shared ls-dyna keyword files to the base of the scratch file system
cp *.k "$SCRATCH_DIR"
cd "$SCRATCH_DIR" || { log_Error "Failed to cd into $SCRATCH_DIR"; exit 1; }

echo "$sims_left_to_run" > SimsLeft.txt



while :; do
    
    # Update job counts
    check_jobs
    
    sims_left_to_run=$(cat SimsLeft.txt)
    log " "
    log "We have $sims_left_to_run simulations of $total_sims left to run"
    log " "
    
    # Check if we have fallen below the max simulation submission limit. If we are, submit more jobs.
    if [ "$total_jobs" -lt "$MAX_JOBS" ]; then
        submit_jobs
    fi

    
    # Check to see if any simulations have finished running (rebound is complete) so that it may be terminated early
    end_if_complete
    sleep 600
    
    # Update job counts
    check_jobs
    
    
    # Check if all jobs are done by seeing if only the master run script (this script) is running and none are queued
    if [ "$total_jobs" -eq 0 ] && [ "$cluster_status" -eq 0 ]; then
    
        # Try one more time to submit any sim that hasn't yet run
        submit_jobs
        sleep 200
        check_jobs
        
        
        if [ "$total_jobs" -eq 0 ] && [ "$cluster_status" -eq 0 ]; then
            log " "
            log "No more jobs are running or queued. Exiting."
            exit 0
        fi
    fi
    
done




