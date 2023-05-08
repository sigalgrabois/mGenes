import json
import os

print("Welcome to the rs catalog file Parser - LevHai parser rs_catalog.py version 1.0")
input("Press any key to continue")
# Define the output filename and path
output_path = os.path.join(os.getcwd(), "config", "rs catalog.json")

# Open the input file
with open('rs catalog.txt', 'r') as f:
    # Initialize a list to hold the parsed data
    rs_data = []
    # check that the headers are correct - if not, raise an error
    correct_headers = ['rs#', 'CHROM', 'POS', 'REF', 'ALT', 'Assay_ID', 'Assay_Name', 'Amplicon START', 'Amplicon END',
                       'Amplicon Len']
    # split the headers by whitespace - only by tab and not by space
    headers = f.readline().strip().split('\t')
    for i in range(len(headers)):
        if headers != correct_headers:
            raise ValueError(f"Error: unexpected header in file. {headers[i]} instead of {correct_headers[i]}")

    # Iterate through each line in the file
    for line in f:
        # Split the line by whitespace
        fields = line.strip().split()

        # Create a dictionary to hold the data for this rs#
        rs_dict = {
            "rs#": fields[0],
            "CHROM": fields[1],
            "POS": fields[2],
            "REF": fields[3],
            "ALT": fields[4],
            "Assay_ID": fields[5],
            "Assay_Name": fields[6],
            "Amplicon START": fields[7],
            "Amplicon END": fields[8],
            "Amplicon Len": fields[9]
        }

        # Append the dictionary to the list
        rs_data.append(rs_dict)

# Convert the list of dictionaries to JSON
rs_json = json.dumps(rs_data)

# Write the JSON string to a file
with open(output_path, 'w') as out_file:
    out_file.write(rs_json)
# ask to press any key to exit
print("Done parsing rs catalog file")
input("Press any key to exit")