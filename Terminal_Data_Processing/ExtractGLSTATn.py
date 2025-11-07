import os
import csv
import re

# Define paths
output_csv = "glstat_summary.csv"  # Running script inside COLLECTED_GLSTAT

# Regex patterns to match the required data
kinetic_energy_pattern = re.compile(r"(?<!eroded )kinetic energy\.+\s+([\d\.\-+Ee]+)")
sliding_interface_pattern = re.compile(r"sliding interface energy\.+\s+([\d\.\-+Ee]+)")
internal_energy_pattern = re.compile(r"(?<!eroded )internal energy\.+\s+([\d\.\-+Ee]+)")
eroded_kinetic_pattern = re.compile(r"eroded kinetic energy\.+\s+([\d\.\-+Ee]+)")
eroded_internal_pattern = re.compile(r"eroded internal energy\.+\s+([\d\.\-+Ee]+)")
total_energy_pattern = re.compile(r"total energy\.+\s+([\d\.\-+Ee]+)")
energy_ratio_wo_erode_pattern = re.compile(r"energy ratio w/o eroded energy\.+\s+([\d\.\-+Ee]+)")

# Updated folder pattern to extract theta values from the filename
folder_pattern = re.compile(r"glstat_(-?[\d]+)_(-?[\d]+)_(-?[\d]+)")

# List to store results
results = []

# Process each glstat file in the current directory
for filename in os.listdir("."):
    if filename.startswith("glstat_"):
        file_path = filename  # Since we are running in COLLECTED_GLSTAT

        # Extract theta values from filename
        match = folder_pattern.search(filename)
        if not match:
            print(f"Warning: Could not extract theta values from {filename}")
            continue
        
        theta_1 = f"{match.group(1)}"
        theta_2 = f"{match.group(2)}"
        theta_3 = f"{match.group(3)}"

        first_kinetic_energy = None
        last_kinetic_energy = None
        last_internal_energy = None
        last_sliding_interface = None
        last_eroded_kinetic = None
        last_eroded_internal = None
        last_total_energy = None
        last_energyRatio_wo_erosion = None

        # Read file and extract values
        with open(file_path, "r") as f:
            for line in f:
                # Find first occurrence of kinetic energy
                if first_kinetic_energy is None:
                    kinetic_match = kinetic_energy_pattern.search(line)
                    if kinetic_match:
                        first_kinetic_energy = kinetic_match.group(1)

                # Find last occurrences of required values
                
                kinetic_match = kinetic_energy_pattern.search(line)
                if kinetic_match:
                    last_kinetic_energy = kinetic_match.group(1)
                    
                internal_match = internal_energy_pattern.search(line)
                if internal_match:
                    last_internal_energy = internal_match.group(1)

                sliding_match = sliding_interface_pattern.search(line)
                if sliding_match:
                    last_sliding_interface = sliding_match.group(1)

                eroded_kinetic_match = eroded_kinetic_pattern.search(line)
                if eroded_kinetic_match:
                    last_eroded_kinetic = eroded_kinetic_match.group(1)

                eroded_internal_match = eroded_internal_pattern.search(line)
                if eroded_internal_match:
                    last_eroded_internal = eroded_internal_match.group(1)

                total_energy_match = total_energy_pattern.search(line)
                if total_energy_match:
                    last_total_energy = total_energy_match.group(1)
                
                energyRatio_wo_erosion_match = energy_ratio_wo_erode_pattern.search(line)
                if energyRatio_wo_erosion_match:
                    last_energyRatio_wo_erosion = energyRatio_wo_erosion_match.group(1)

        # Ensure valid extraction of all required values
        if None in [first_kinetic_energy, last_kinetic_energy, last_internal_energy, last_sliding_interface, 
                    last_eroded_kinetic, last_eroded_internal, last_total_energy, last_energyRatio_wo_erosion]:
            print(f"Warning: Missing values in {filename}")
            continue

        # Store the extracted data
        results.append([theta_1, theta_2, theta_3, first_kinetic_energy, last_internal_energy, last_kinetic_energy, 
                        last_sliding_interface, last_eroded_kinetic, last_eroded_internal, last_total_energy, last_energyRatio_wo_erosion])

# Write results to CSV
with open(output_csv, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["theta_1", "theta_2", "theta_3","initial_ke", 
                     "final_internal_energy", "final_kinetic_energy", "final_sliding_interface_energy",
                     "final_eroded_kinetic_energy", "final_eroded_internal_energy", "final_total_energy", "final_energy_ratio_without_erosion"])
    writer.writerows(results)

print(f"CSV file created: {output_csv}")

