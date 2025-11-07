#!/usr/bin/env bash

# Get the directory from which the script is being run
CampaignFolder=$(pwd)
BASEkFilePath="$CampaignFolder/BASE/BASE.k"
BASEshFilePath="$CampaignFolder/BASE/BASE.sh"
BASEshPostProcPath="$CampaignFolder/BASE/BASE_PostProcess.sh"
zShift_file="$CampaignFolder/Z_shift_matrix.csv"

MinRunTime=0.0005


# Find the distance which the particle should be shifted normal to the target to get it as close to the target at the
# start of the simulation as possible (to minimize run-time).
# NOTE: The file Z_shift_matrix.csv must exist in the same directory as this shell script as it stores the 
# shift values for each orientation.

find_z_shift() {
    local theta1="$1"
    local theta2="$2"
    local theta3="$3"

    # Search for the matching line in the CSV
    awk -F, -v t1="$theta1" -v t2="$theta2" -v t3="$theta3" \
    'BEGIN {found=0} 
     $2==t1 && $3==t2 && $4==t3 {printf "%.4f\n", $1; found=1; exit} 
     END {if (found==0) print "0"}' "$zShift_file"
}


create_and_populate_directory() {
    local theta1="$1"
    local theta2="$2"
    local theta3="$3"
    local dir_name="MAIN_${theta1//[.]/p}_${theta2//[.]/p}_${theta3//[.]/p}"

    mkdir -p "$CampaignFolder/$dir_name"
    
    z_shift=$(find_z_shift "$theta1" "$theta2" "$theta3")
    
    echo "$z_shift" > "$CampaignFolder/$dir_name/z_shift_value.txt"


    # Process .k file
    sed "s/THETA1/$theta1/;s/THETA2/$theta2/;s/THETA3/$theta3/;s/Z_SHIFT/$z_shift/" "$BASEkFilePath" | \
    sed "s/BASE/$dir_name/" > "$CampaignFolder/$dir_name/${dir_name}.k"


    # Process .sh file
    sed "s|BASE|$dir_name|g" "$BASEshFilePath" | \
    sed "s|current_working_directory|$CampaignFolder/$dir_name|g" > "$CampaignFolder/$dir_name/${dir_name}.sh"
    

    # Process PostProcess .sh file efficiently
    sed -e "s|BASE|$dir_name|g" \
        -e "s|current_working_directory|$CampaignFolder/$dir_name|g" \
        -e "s|MIN_RUN_TIME|$MinRunTime|g" \
        "$BASEshPostProcPath" > "$CampaignFolder/$dir_name/${dir_name}_PostProcess.sh"

}

CSV_INPUT_FILE="$CampaignFolder/theta_values.csv"  # or wherever the CSV is

if [ ! -f "$CSV_INPUT_FILE" ]; then
    echo "ERROR: CSV input file not found at $CSV_INPUT_FILE"
    exit 1
fi


while IFS=',' read -r theta1 theta2 theta3; do
    # Skip empty or malformed lines
    if [[ -n "$theta1" && -n "$theta2" && -n "$theta3" ]]; then
        create_and_populate_directory "$theta1" "$theta2" "$theta3"
    fi
done < "$CSV_INPUT_FILE"






