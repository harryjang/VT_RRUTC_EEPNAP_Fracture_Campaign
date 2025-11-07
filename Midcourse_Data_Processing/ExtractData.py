import os
import re
import csv
import shutil
import subprocess



#### Read Me ####

# For this code to run as intended, the following must be true

# 1) One folder back from the location of this python code, you must have a single keyword file corresponding to the DYNA simulation being analyzed,
#    all of the message files (this Python code assumes the DYNA simulations are run using MPP which produce 1 message file per compute core used in
#    the analysis -- these message files are only used to acces information on which elements erode during the simulation), and the binout file (this
#    code assumes MPP which compress all ASCII output into a binary file called binout#### where #### are a string of integers. This code assumes that
#    there is only one such binout file called binout0000. Modifications must be made if other binout files also have to be read).

# 2) At the location specified as l2a_executable_path, must have the executable l2a.exe (which comes with DYNA MPP executables). This executable converts
#   the above mentioned binout#### files to human-readable ASCII files.


# 3) The binout file must contain NODOUT information for the analysis.

# OTHER ASSUMPTIONS MADE BY PROGRAM

#  1) The LS-DYNA part being analyzed is part 2.
#  2) The time at which part 2 is to be analyzed corresponds to the last time reported in the NODOUT file.





# ============ Define regular expression patterns ==================

# Regular expression pattern to match "integer failed at time"
pattern = re.compile(r'(\d+)\s+failed at time')

# Using regex to filter files that start with "mes"
mes_pattern = re.compile(r'^mes.*')


# ====================== Define paths ===========================

current_directory = os.getcwd()
# Go one step back from the current script directory to search for raw simulation files
search_directory = os.path.dirname(current_directory)





# Get the current working directory
current_dir = os.getcwd()

# Define the name of the folder you want to find the parent of
folder_name = "Analyze"

# Traverse the directory tree upwards from the current directory
while True:
    # Check if the folder exists in the current directory
    if folder_name in os.listdir(current_dir):
        # Store the path to the parent directory of the Analyze folder
        root_directory = os.path.abspath(current_dir)
        break
    # Move up one level in the directory tree
    current_dir = os.path.abspath(os.path.join(current_dir, os.pardir))






output_directory = os.path.join(current_directory, 'Input_files')
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Fixed location of "l2a.exe" executable
# # # l2a_executable_path = r"C:\Users\TPL\Downloads\Peridynamics_Via_DYNA\l2a.exe"





def extract_element_lines(input_path):
    extracted_lines = []
    with open(input_path, 'r') as f:
        lines = f.readlines()
        flag = False
        for line in lines:
            # Check for the "*ELEMENT_SOLID" keyword
            if "*ELEMENT_SOLID" in line:
                flag = True
            elif flag and len(line) > 15 and "2" in line[8:16]:
                extracted_lines.append(line)
            # Reset flag if another line with '*' is found
            elif "*" in line:
                flag = False
    return extracted_lines



def save_elementData_to_csv(extracted_lines, output_path):
    with open(output_path, 'w') as csvfile:
        for line in extracted_lines:
            # Split the line into chunks of 8 characters and strip spaces
            chunks = [line[i:i+8].strip() for i in range(0, len(line), 8)]
            
            # Remove the second column
            if len(chunks) > 1:
                del chunks[1]
            
            # Remove any empty chunks (this should handle the trailing comma issue)
            chunks = [chunk for chunk in chunks if chunk]

            # Write the chunks as a comma-separated line
            csvfile.write(",".join(chunks) + "\n")



def extract_from_nodout(input_path):
    data_to_process = []
    with open(input_path, 'r') as f:
        lines = f.readlines()
        start_idx = None
        end_idx = len(lines)
        
        # Find the last occurrence of a line containing both "nodal point" and "x-coor"
        for idx, line in reversed(list(enumerate(lines))):
            if "nodal point" in line and "x-coor" in line:
                start_idx = idx
                break
        
        # Find the next occurrence of " n o d a l   p r i n t   o u t" after the identified line
        for idx, line in enumerate(lines[start_idx:], start=start_idx):
            if " n o d a l   p r i n t   o u t" in line:
                end_idx = idx
                break
        
        if start_idx is not None:
            extracted_lines = lines[start_idx + 1:end_idx]
            print(end_idx)
            # Filtering out empty or whitespace-only lines
            data_to_process = [line for line in extracted_lines if line.strip()]

    return data_to_process




def save_nodePositions_to_csv(data, output_path):
    with open(output_path, 'w') as csvfile:
        for line in data:
            first_col = line[:10].strip()
            # Extracting columns and directly filtering out any empty ones
            other_cols = [line[i:i+12].strip() for i in range(10, len(line), 12) if line[i:i+12].strip()]
            
            # Selecting the first and the last three columns
            selected_cols = [first_col] + other_cols[-3:]
            
            csvfile.write(",".join(selected_cols) + "\n")


def save_nodeV_to_csv(data, output_path):
    with open(output_path, 'w') as csvfile:
        for line in data:
            first_col = line[:10].strip()
            # Extracting columns and directly filtering out any empty ones
            other_cols = [line[i:i+12].strip() for i in range(10, len(line), 12) if line[i:i+12].strip()]
            
            # Selecting the first and the last three columns
            selected_cols = [first_col] + other_cols[-9:-6]
            
            csvfile.write(",".join(selected_cols) + "\n")



def extract_failed_elements(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            match = pattern.search(line)
            if match:
                yield int(match.group(1))



def find_file_with_extension(directory, extension):
    """
    Finds a file with the given extension in the directory.
    Returns the file path if only one is found, otherwise raises an exception.
    """
    # Adjusting for case-insensitivity on Windows
    files_with_extension = [file for file in os.listdir(directory) if file.lower().endswith(extension.lower())]
    
    if len(files_with_extension) == 1:
        return os.path.join(directory, files_with_extension[0])
    elif len(files_with_extension) > 1:
        raise Exception(f"Multiple files with the extension '{extension}' found!")
    else:
        raise Exception(f"No file with the extension '{extension}' found!")





def main():
    
    try:
        # find the .k file and get its path
        keyword_path = os.path.join(root_directory, "projectile_geometry.k")

        nodout_filename = None
        for file in os.listdir(search_directory):
            if file.lower().startswith("nodout"):
                nodout_filename = file
                break

        if nodout_filename:
            nodout_path = os.path.join(search_directory, nodout_filename)
        else:
            raise Exception("No file starting with 'nodout' found in the directory!")


    except Exception as e:
        print(f"Error: {e}")


        
    # --------------Extract element data------------------
    
    output_path = os.path.join(output_directory, 'element_data.csv')
    
    extracted_lines = extract_element_lines(keyword_path)
    save_elementData_to_csv(extracted_lines, output_path)
    print(f"Data saved to {output_path}")


    #------------- Extract node data--------------------
  
    data_to_process = extract_from_nodout(nodout_path)

    output_path = os.path.join(output_directory, 'node_positions.csv')
    save_nodePositions_to_csv(data_to_process, output_path)
    print(f"Data saved to {output_path}")
    
    output_path = os.path.join(output_directory, 'node_velocities.csv')
    save_nodeV_to_csv(data_to_process, output_path)   
    print(f"Data saved to {output_path}")


    #------------Extract failed (eroded) elements--------



    # List all files in the directory that start with "mes."
    filenames = [f for f in os.listdir(root_directory) if mes_pattern.match(f) and os.path.isfile(os.path.join(root_directory, f))]

    
    # Create or overwrite the CSV file
    output_path = os.path.join(output_directory, 'deleted_elements_data.csv')
    with open(output_path, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)

        # Extract failed times from the files
        for filename in filenames:
            full_path = os.path.join(root_directory, filename)
            for elementID in extract_failed_elements(full_path):
                csvwriter.writerow([elementID])
    print(f"Data saved to {output_path}")
    

if __name__ == "__main__":
    main()
