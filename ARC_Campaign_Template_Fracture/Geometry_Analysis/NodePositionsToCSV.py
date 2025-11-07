# Path to your input file
file_path = 'projectile_node_coordinates.txt'


# Define the path for the output CSV file
output_csv_path = 'projectile_node_coordinates.csv'

# Open the input file and the output file
with open(file_path, 'r') as infile, open(output_csv_path, 'w') as outfile:
    for line in infile:
        # Split the line into columns based on whitespace
        columns = line.split()
        
        # Check if the line has at least 4 columns to avoid index errors
        if len(columns) >= 4:
            # Extract columns 2, 3, and 4 (using 0-based indexing)
            selected_columns = [columns[1], columns[2], columns[3]]
            
            # Write the extracted columns to the output file, separated by commas
            outfile.write(','.join(selected_columns) + '\n')

print(f"Extracted columns saved to {output_csv_path}")
