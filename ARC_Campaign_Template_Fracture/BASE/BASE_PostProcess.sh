#!/bin/bash


# Always remove lock file if script exits
trap 'rm -f PostprocessLock' EXIT



# !!!!! Script is run in the SCRATCH FOLDER holding the sim output


# Define paths
OG_SIM_DIR="current_working_directory"
OG_PARENT_DIR=$(dirname "$OG_SIM_DIR")
CAMPAIGN_LOG_DIR="$OG_PARENT_DIR/Automated_Campaign_Output"
SUBMISSION_SCRIPT_LOG="$CAMPAIGN_LOG_DIR/submission_script_log.txt"
POSTPROCESS_STATUS_LOG="$CAMPAIGN_LOG_DIR/post_process_status_log.txt"
ERROR_LOG="$CAMPAIGN_LOG_DIR/error_log.txt"

# Ensure log files exist
touch $POSTPROCESS_STATUS_LOG
touch $ERROR_LOG


# Define log functions
log_PostProcess() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $@" >> "$POSTPROCESS_STATUS_LOG"
}

log_Error() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $@" >> "$ERROR_LOG"
}


# Define variables
script_path="/home/harryjang/Simulations/Campaign_Template_Fracture/Midcourse_Data_Processing"
sim_directory=$(pwd)
MinRunTime=MIN_RUN_TIME
STATUS_FILE="last_processing_time.txt"
TIME_FILE="last_nodout_time.txt"
analyze_directory="$sim_directory/Analyze"



# Helper functions
extract_binout() {
    binout_files=("$sim_directory"/binout*)
    if [[ -f "${binout_files[0]}" ]]; then
        "$script_path/ExtractBINOUT.l2a" "${binout_files[@]}"
    else
        return 1
    fi
}


run_nodout_extraction() {
    cp "$script_path/WriteLastNodoutBlock.py" "$sim_directory"
    python WriteLastNodoutBlock.py 2> >(while read -r line; do log_Error "$line"; done)
    if [[ $? -ne 0 ]]; then
        return 1
    fi
}




move_nodout_files() {
    mkdir -p "$analyze_directory"
    for nodout_file in "$sim_directory"/nodout*.txt; do
        [[ -f $nodout_file ]] || continue
        base_name=$(basename "$nodout_file" .txt)
        target_directory="$analyze_directory/$base_name"
        mkdir -p "$target_directory"
        mv "$nodout_file" "$target_directory"
    done
}




process_nodout_directories() {
    for dir in "$analyze_directory"/nodout*; do
        [[ -d $dir && "$(basename "$dir")" != "nodout_0" ]] || continue
        get_clusters_dir="$dir/GetClusters"
        mkdir -p "$get_clusters_dir"
        cp "$script_path/ExtractData.py" "$get_clusters_dir"
        # cp "$script_path/MATLAB_SlurmSub.sh" "$get_clusters_dir"
        cp "$script_path/RunScript.m" "$get_clusters_dir"
        cp "$script_path/ClusterAlgorithmParallelTermination.m" "$get_clusters_dir"

        cd "$get_clusters_dir"
        
        # Run ExtractData.py
        python ExtractData.py 2> >(while read -r line; do log_Error "$line"; done)
        if [[ $? -ne 0 ]]; then
            return 2  # Return 2 for Python script failure
        fi

        # Run MATLAB analysis
        module load MATLAB
        log_PostProcess "Begin MATLAB fragment analysis on $sim_directory"
        log_PostProcess " "
        matlab -batch RunScript
        if [[ $? -ne 0 ]]; then
            check_matlab_errors
            return 3  # Return 3 for MATLAB failure
        fi

        check_matlab_errors
        cd "$sim_directory"
    done
    return 0  # Success
}



check_matlab_errors() {
    if [[ -f error_log.txt ]]; then
        while IFS= read -r line; do
            log_Error "$line"
        done < error_log.txt
        # Clear the error file
        > error_log.txt
    fi
}





# MAIN() function

# Convert binout files to ASCII output files
if ! extract_binout; then
    log_Error "No binout files in $sim_directory."
    rm -f PostprocessLock
    exit 0
fi


# Extract the last time block from the nodout file and record the time corresponding to the last nodout recording to the file TIME_FILE
if ! run_nodout_extraction; then
    log_Error "WriteLastNodoutBlock.py script failed for $sim_directory"
    rm -f PostprocessLock
    exit 0
fi


# Read last processing time
if [[ -f "$STATUS_FILE" ]]; then
    last_processing_time_numeric=$(cat "$STATUS_FILE")
else
    echo "0.0" > "$STATUS_FILE"
    last_processing_time_numeric=0.0
fi


# Read nodout time produced by the Python script WriteLastNodoutBlock.py
if [[ ! -f "$TIME_FILE" ]]; then
    log_Error "File $TIME_FILE does not exist for $sim_directory. Exiting script."
    rm -f PostprocessLock
    exit 0
fi
# Read the number from the file
nodout_time=$(cat "$TIME_FILE")
# Convert scientific notation to a numeric value with precision
nodout_time_numeric=$(awk '{printf "%.10f\n", $1}' <<< "$nodout_time")
# Check if conversion was successful
if [[ -z "$nodout_time_numeric" || ! "$nodout_time_numeric" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    log_Error "Unable to extract a valid numeric value from $TIME_FILE"
    rm -f PostprocessLock
    exit 0
fi


# Compute time-based conditions for processing based on the current nodout time (nodout_time_numeric) and the last time at which particle fracture information has been analyzed (last_processing_time_numeric)

# Slightly increase LastProcessing_time to avoid numerical errors which may lead to incorrect conclusion that current nodout time and the previously processed nodout time are not equal when they really are
adjusted_LastProcessing_time=$(echo "$last_processing_time_numeric * (1 + $last_processing_time_numeric / 1000)" | bc -l)

condition1=$(echo "$nodout_time_numeric >= $MinRunTime" | bc -l)
condition2=$(echo "$nodout_time_numeric > $adjusted_LastProcessing_time" | bc -l)

CurrentSim=$(basename "$(pwd)")
log_PostProcess " SIMULATION: $CurrentSim"
log_PostProcess "Runtime state:     $(pwd)---Current nodout file time: $nodout_time_numeric --------------------- Minimum time to start analysis: $MinRunTime ---------------------------- Start analysis conditional: $condition1"
log_PostProcess "Runtime state:     $(pwd)---Current nodout file time: $nodout_time_numeric ---------- Previous time at which file was processed: $adjusted_LastProcessing_time ---------- Start analysis conditional: $condition2"
log_PostProcess "  "

if [ "$condition1" -eq 1 ] && [ "$condition2" -eq 1 ]; then
    echo "$nodout_time_numeric" > "$STATUS_FILE"
    move_nodout_files
    
    process_nodout_directories
    exit_status=$?

    if [[ $exit_status -eq 2 ]]; then
        log_Error "Critical error: ExtractData_V2.py script failed for $sim_directory"
        rm -f PostprocessLock
        exit 0
    elif [[ $exit_status -eq 3 ]]; then
        log_Error "Critical error: MATLAB cluster extraction script failed for $sim_directory"
        rm -f PostprocessLock
        exit 0
    elif [[ $exit_status -ne 0 ]]; then
        log_Error "Unknown error in process_nodout_directories(). Exiting."
        rm -f PostprocessLock
        exit 0
    fi
fi



# Check for termination status
if [[ -f TERMINATION_STATUS.txt && $(cat TERMINATION_STATUS.txt) == "y" ]]; then
    log_PostProcess "$(pwd) has completed. Waiting 60 minutes before termination."
    sleep 3600
    
    #Extract cluster information at end of simulation
    cd "$sim_directory"
    
    # Convert binout files to ASCII output files
    if ! extract_binout; then
        log_Error "Grace period is over, but no binout files found in $sim_directory. Will transfer previously recorded simulation data and quit."
        echo sw1 > D3KIL
        sleep 5
        exit 0
    fi
    
    # Extract the last time block from the nodout file and record the time corresponding to the last nodout recording to the file TIME_FILE
    if ! run_nodout_extraction; then
        log_Error "Grace period is over, but WriteLastNodoutBlock.py script failed for $sim_directory. Will transfer previously recorded simulation data and quit."
        echo sw1 > D3KIL
        sleep 5
        exit 0
    fi
    
    move_nodout_files
    
    process_nodout_directories
    exit_status=$?

    if [[ $exit_status -eq 2 ]]; then
        log_Error "Grace period is over, but ExtractData_V2.py script failed for $sim_directory. Will transfer previously recorded simulation data and quit."
        echo sw1 > D3KIL
        sleep 5
        exit 0
    elif [[ $exit_status -eq 3 ]]; then
        log_Error "Grace period is over, but MATLAB cluster extraction script failed for $sim_directory. Will transfer previously recorded simulation data and quit."
        echo sw1 > D3KIL
        sleep 5
        exit 0
    elif [[ $exit_status -ne 0 ]]; then
        log_Error "Grace period is over, but the function process_nodout_directories failed for $sim_directory for an UNKNOWN reason. Will transfer previously recorded simulation data and quit."
        echo sw1 > D3KIL
        sleep 5
        exit 0
    fi
    
    echo sw1 > D3KIL
    sleep 5
    log_PostProcess "$(pwd) Grace period has completed. Simulation completed."
    log_PostProcess ""
    exit 0
fi

rm -f PostprocessLock
exit 0


