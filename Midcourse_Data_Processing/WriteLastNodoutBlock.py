import re

def extract_last_non_rot_block(filename, output_filename="nodout_1.txt", time_filename="last_nodout_time.txt"):
    with open(filename, 'r', encoding='utf-8') as file:
        lines = file.readlines()
    
    blocks = []
    current_block = []
    delimiters = []

    for line in lines:
        line = line.rstrip("\n")  # Remove newline character
        
        if "time" in line:  # First delimiter
            if len(delimiters) == 3 and "rot" not in delimiters[2]:  # Save previous block if valid
                blocks.append((delimiters, current_block))
            delimiters = [line]  # Reset and store first delimiter
            current_block = []
        
        elif line == "" and len(delimiters) == 1:  # Second delimiter (blank line)
            delimiters.append(line)
        
        elif len(delimiters) == 2:  # Third delimiter
            delimiters.append(line)
        
        elif len(delimiters) == 3:  # Actual text block content
            current_block.append(line)
    
    # Final check in case the last block is valid
    if len(delimiters) == 3 and "rot" not in delimiters[2]:
        blocks.append((delimiters, current_block))
    
    if blocks:
        last_non_rot_block = blocks[-1]  # Get the last valid block
        delimiters, text_block = last_non_rot_block

        # Extract time from the first delimiter line
        time_match = re.search(r"time\s+([\d\.eE+-]+)", delimiters[0])  # Handles scientific notation too
        time_value = time_match.group(1) if time_match else "Unknown"

        # Write time to file
        with open(time_filename, 'w', encoding='utf-8') as time_file:
            time_file.write(time_value + "\n")

        # Write block to nodout_1.txt
        with open(output_filename, 'w', encoding='utf-8') as out_file:
            out_file.write("\n".join(delimiters) + "\n")  # Write delimiters with spacing
            out_file.write("\n".join(text_block) + "\n")  # Write text block
        
        print(f"Last extracted block saved to {output_filename}")
        print(f"Time value ({time_value}) saved to {time_filename}")
    else:
        print("No valid block found.")

# Example usage
filename = "nodout"
extract_last_non_rot_block(filename)
