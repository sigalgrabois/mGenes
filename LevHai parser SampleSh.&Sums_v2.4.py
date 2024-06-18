import json
import os
import pandas as pd
# import from the same folder from LevHai parser VCF.py
# import function
from LevHai_parser_VCF import remove_suffix, replace_dash_with_underscore

"""
This script reads the data from the text files and converts it into a JSON array.

"""


def find_files():
    """
    The find_files function searches the current directory for files with a .txt extension.
        It then checks to see if the file contains 'TXA' or 'rs' in line 1 of the text file.
        If it does, it appends that filename to either TXA_files or rs_files list respectively.

    :return: A list of txa files, a list of rs files and the path to the samplesheet file
    :doc-author: Sigal and Chen
    """
    TXA_files = []
    rs_files = []
    file_path_SampleSheet = ''

    # Find the text file in the current path location
    for file in os.listdir():
        # Check if the file is a text file and contains 'TXA' in line 1
        if file.endswith('.txt') and 'TXA' in open(file).readlines()[1] and 'P1' in file:
            TXA_files.append(file)
            print("found TXA QC P1 file")
            print(TXA_files[0])
        # Check if the file is a text file and contains 'rs' in line 1
        elif file.endswith('.txt') and 'TXA' in open(file).readlines()[1] and 'P2' in file:
            TXA_files.append(file)
            print("found TXA QC P2 file")
            print(TXA_files[1])
            continue
        # Check if the file is a text file and contains 'rs' in line 1
        if file.endswith('.txt') and 'rs' and not 'TXA' in open(file).readlines()[1] and 'P1' in file:
            rs_files.append(file)
            print("found rs QC P1 file")
            print(rs_files[0])
            continue
        # Check if the file is a text file and contains 'rs' in line 1
        elif file.endswith('.txt') and 'rs' and not 'TXA' in open(file).readlines()[1] and 'P2' in file:
            rs_files.append(file)
            print("found rs QC P1 file")
            print(rs_files[1])
            continue
        # Find the SampleSheet file in the current path location
        if file.endswith('.csv') and 'Sample_ID' in open(file).readlines()[18]:
            file_path_SampleSheet = file
            print("found SampleSheet file")
            print(file_path_SampleSheet)
            continue
    # Return the lists of files and the path to the SampleSheet file
    return TXA_files, rs_files, file_path_SampleSheet


def read_convert_Data(file_path, flag, plate):
    """
    The read_convert_Data function reads the data from a text file and converts it into a JSON array.
    The function takes three arguments:
        1) The path to the input text file, which must be in tab-delimited format with column headers.
        2) A flag indicating whether the data is for TXA or rs (reign name). This determines how each row is converted.
           If 'TXA', then each row will be converted into an array of three elements: [rs_TXA, Sample ID, value].
           If 'rs', then each row will be converted into an array of three

    :param file_path: Specify the path to the input file
    :param flag: Determine which data is being read
    :param plate: Determine the output file name
    :return: A list of lists, which is the converted data
    :doc-author: Sigal and Chen
    """
    # Read the data from the text file
    with open(file_path, 'r') as f:
        data = [line.strip().split('\t') for line in f]

    # Initialize a new empty list to store the converted rows
    converted_data = []

    if flag == 'TXA':
        # Add the table headers to the converted data
        converted_data.append(['rs_TXA', 'Sample ID', 'value'])

        # Iterate through the rows and convert each row into the desired format
        for row in data[1:]:  # skip the header row
            for i in range(1, len(data[0])):  # skip the first column (rs)
                converted_data.append([row[0], data[0][i], row[i]])

        # Create the output folder if it does not exist
        output_folder = os.path.join(os.path.dirname(file_path), 'config')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Generate the output file path
        if plate == 'P1':
            output_file_path = os.path.join(output_folder, 'converted_data_TXA_P1.json')
        elif plate == 'P2':
            output_file_path = os.path.join(output_folder, 'converted_data_TXA_P2.json')
        # Part 3: Write the converted data to the output file, handling duplicates

        # Read the existing data from the output file, if it exists
        existing_data = []
        if os.path.exists(output_file_path):
            with open(output_file_path, 'r') as f:
                existing_data = json.load(f)

        # Initialize counters for new rows, duplicates, and unique Sample IDs
        new_rows = 0
        duplicates = 0
        unique_sample_ids = set(row[1] for row in existing_data if len(row) > 1)

        # Iterate through the converted data and add each row to the existing data
        # if it is not a duplicate
        for row in converted_data[1:]:  # skip the header row
            if row not in existing_data:
                existing_data.append(row)
                new_rows += 1
                unique_sample_ids.add(row[1])  # add the Sample ID to the set
            else:
                duplicates += 1

        # Add the table headers to the existing data, if they are not already present
        if ['rs_TXA', 'Sample ID', 'value'] not in existing_data:
            existing_data.insert(0, ['rs_TXA', 'Sample Name', 'Depth value'])

        # Write the updated existing data to the output file as a JSON array
        with open(output_file_path, 'w') as f:
            json.dump(existing_data, f)
    elif flag == 'rs':
        # Add the table headers to the converted data
        converted_data.append(['Reign Name', 'Sample Name', 'Depth value'])

        # Iterate through the rows and convert each row into the desired format
        for row in data[1:]:  # skip the header row
            for i in range(1, len(data[0])):  # skip the first column (rs)
                # use the functions replace underscore and replace space to replace the underscore and space in the sample name
                # and use remove suffix to remove the suffix of the sample name
                converted_data.append([row[0], data[0][i], row[i]])
                print([row[0], data[0][i], row[i]])

        # Part 2: Create the output folder and file path

        # Create the output folder if it does not exist
        output_folder = os.path.join(os.path.dirname(file_path), 'config')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Generate the output file path
        if plate == 'P1':
            output_file_path = os.path.join(output_folder, 'converted_data_rs_P1.json')
        elif plate == 'P2':
            output_file_path = os.path.join(output_folder, 'converted_data_rs_P2.json')

        # Part 3: Write the converted data to the output file, handling duplicates

        # Read the existing data from the output file, if it exists
        existing_data = []
        if os.path.exists(output_file_path):
            with open(output_file_path, 'r') as f:
                existing_data = json.load(f)

        # Initialize counters for new rows, duplicates, and unique Sample IDs
        new_rows = 0
        duplicates = 0
        unique_sample_ids = set(row[1] for row in existing_data)

        # Iterate through the converted data and add each row to the existing data
        # if it is not a duplicate
        for row in converted_data[1:]:  # skip the header row
            if row not in existing_data:
                existing_data.append(row)
                new_rows += 1
                unique_sample_ids.add(row[1])  # add the Sample ID to the set
                print(row)
            else:
                duplicates += 1

        # Add the table headers to the existing data, if they are not already present
        if ['Reign Name', 'Sample Name', 'Depth value'] not in existing_data:
            existing_data.insert(0, ['Reign Name', 'Sample Name', 'Depth value'])

        # Write the updated existing data to the output file as a JSON array
        with open(output_file_path, 'w') as f:
            json.dump(existing_data, f)

    # Print a summary of the number of new rows, duplicates, and unique Sample IDs
    print(f"{new_rows} new rows added")
    print(f"{duplicates} duplicates removed")
    print(f"{len(unique_sample_ids)} unique Sample IDs in total")


def read_convert_SampleSheet(path_Sample_Sheet):
    """
    The read_convert_SampleSheet function reads the SampleSheet.csv file and converts it to a json object
        which is then written to a json file in the config folder.

    :param path_Sample_Sheet: Specify the path to the samplesheet
    :return: A pandas dataframe
    :doc-author: Sigal
    """

    # Read the CSV file into a Pandas DataFrame
    df = pd.read_csv(path_Sample_Sheet, skiprows=18, header=0, delimiter=',')
    # create new path for the json file
    path_Sample_Sheet_temp = os.path.join(os.path.dirname(path_Sample_Sheet), 'config', 'SampleSheet.json')
    # create a copy of the dataframe
    df_temp = df.copy()
    # convert the temp dataframe to json
    df_temp.to_json(path_Sample_Sheet_temp, orient='records')
    print("SampleSheet.json file created in config folder")


def reforamt_json(path_file, flag, plate):
    """
    The reformat_json function reads the json file and rewrites it in a new format

    :param
    :return:
    :doc-author: Sigal
    """
    # the file name to be opened
    file_name = f"converted_data_{flag}_{plate}.json"
    print(file_name)
    # open the file
    with open(os.path.join(os.path.dirname(path_file), 'config', file_name), 'r') as f:
        data = json.load(f)
        headers = data[0]
        formated_data = []

        for row in data[1:]:
            # for each row change the sample name to the correct format - remove the suffix and replace the underscore and space
            row[1] = remove_suffix(row[1])
            row[1] = replace_dash_with_underscore(row[1])
            row.append(
                row[0] + "_" + row[1])  # Concatenate "Reign Name" and "Sample Name" and append it as a new column
            formated_data.append(dict(zip(headers + ["Key"], row)))  # Include the new column in the dictionary
        with open(os.path.join(os.path.dirname(path_file), 'config', file_name), 'w') as f:
            json.dump(formated_data, f)


if __name__ == "__main__":
    print("Welcome to files Parser - LevHai parser SampleSh.&Sums.py version 2.4")
    input("Press Enter to continue...")
    print("starting the program\n", "finding the files\n")
    list_path_TXA, list_path_rs, path_Sample_Sheet = find_files()
    print("files found: \n")
    print(list_path_TXA, list_path_rs, path_Sample_Sheet)
    plate = ''
    # loop to go through all the files and convert them
    # first one is for TXA
    for path_TXA in list_path_TXA:
        if 'P1' in path_TXA:
            plate = 'P1'
        elif 'P2' in path_TXA:
            plate = 'P2'
        read_convert_Data(path_TXA, 'TXA', plate)
        reforamt_json(path_TXA, 'TXA', plate)
    # second one is for rs
    for path_rs in list_path_rs:
        if 'P1' in path_rs:
            plate = 'P1'
        elif 'P2' in path_rs:
            plate = 'P2'
        read_convert_Data(path_rs, 'rs', plate)
        reforamt_json(path_rs, 'rs', plate)
    read_convert_SampleSheet(path_Sample_Sheet)
    print("Done")
    # print press any key to exit
    input("Press any key to exit")





