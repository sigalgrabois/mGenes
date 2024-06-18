# Importing the necessary libraries
import os
import shutil
from datetime import time

import pandas as pd
import subprocess


def extract_number_from_sample_id(sample_id):
    # Split the sample ID according to "_"
    parts = sample_id.split("_")

    # Return the last item from the split list
    return parts[-1]


def replace_pipe_with_slash(input_file, output_file):
    """
    The replace_pipe_with_slash function replaces the "|" character with "/" in the input file and writes the modified
    contents to the output file.
    :param input_file:
    :param output_file:
    :return: None
    """
    try:
        # Open the input file and read its contents
        with open(input_file, "r") as file:
            contents = file.read()

        # Replace "|" with "/"
        modified_contents = contents.replace("|", "/")

        # Write the modified contents to the output file
        with open(output_file, "w") as file:
            file.write(modified_contents)

        print("File successfully processed.")

    except FileNotFoundError:
        print("The input file was not found.")

    except Exception as e:
        print("An error occurred during file processing:", e)


def copy_files_from_wsl(source_path, destination_path):
    """
    Copy files from WSL to local directory using the subprocess module.

    :param source_path: The path to the file or directory to be copied from WSL.
    :param destination_path: The path to the local directory to copy the file(s) to.
    :return: None
    """
    try:
        # Use the 'cp' command to copy the file(s) from WSL to local directory
        # the file should have the same name as the original file
        subprocess.run(["wsl", "cp", "-r", source_path, destination_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error copying files: {e.cmd}")
        print(f"Output: {e.output}")
        print(f"Return code: {e.returncode}")
        print(f"Error message: {e.stderr}")
    else:
        print(f"Files successfully copied from {source_path} to {destination_path}.")


def run_sub_process(lib_plate):
    """
    The run_sub_process function runs the bcf tools commands in a subprocess.
    The function runs the commands on the converted_data_VCF_PX_Edited.json file.

    :return: The output of the subprocess
    :doc-author: Sigal
    """

    working_dir = input("Enter the working directory (wsl): ")
    # cmd01 = f"bcftools filter -i 'ID=@/mnt/c/Users/User/LevHai/Filter_ToInclude_rs_list_v4.txt' -o '/mnt/c/Users/User/LevHai/{lib_plate}/Partek_{lib_plate}_Filetr.vcf' '/mnt/c/Users/User/LevHai/{lib_plate}/Partek_{lib_plate}.vcf'"
    # cmd02 = f"bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE[\t%TGT][\t%DP]\n' --print-header '/mnt/c/Users/User/LevHai/{lib_plate}/Partek_{lib_plate}_Filetr.vcf' -o '/mnt/c/Users/User/LevHai/{lib_plate}wip.vcf'"
    # cmd01 = f"bcftools filter -i 'ID=@{working_dir}/{rs_list}' -o '{working_dir}/{lib_plate}/Partek_{lib_plate}_Filetr.vcf' '{working_dir}/{lib_plate}/Partek_{lib_plate}.vcf'"
    # cmd02 = f"bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE[\t%TGT][\t%DP]\n' --print-header '{working_dir}/{lib_plate}/Partek_{lib_plate}_Filetr.vcf' -o '{working_dir}/{lib_plate}wip.vcf'"
    # The first command gets as input tthe filtered vcf file and changes the genotype of the samples with DP<35 and saves the output to a temp file
    cmd00 = f"bcftools filter -S . -e 'FMT/DP<35' '{working_dir}/{lib_plate}/Partek_{lib_plate}_Filetr.vcf' > '{working_dir}/{lib_plate}temp.vcf'"
    cmd01 = f"bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE[\t%GT]\n' --print-header '{working_dir}/{lib_plate}temp.vcf' -o '{working_dir}/{lib_plate}Pool.vcf'"


    process = subprocess.Popen(["wsl", "bash", "-c", f"cd {working_dir} && {cmd00} && {cmd01}"],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # get the output of the subprocess
    output, error = process.communicate()
    # print the content of the output
    if error:
        print(error.decode("utf-8"))
    else:
        print(output.decode("utf-8"))
    copy_files_from_wsl(working_dir, "archive")


def run(lib_plate):
    # rs_result file:
    # Reading the rs_catalog file and converting it into a pandas dataframe
    rs_catalog = pd.read_csv("rs_list_pool.txt", header=None)
    rs_catalog.columns = ["rs"]
    # extract the lib number from the lib_plate when the lib_plate is in the format of "L02P1"
    # split the lib_plate by the "P" and take the first element
    lib_num = lib_plate.split("P")[0]

    # Reading the L6.vcf file and converting it into a pandas dataframe
    vcf_file = f'{lib_plate}Pool_edited2.vcf'
    vcf_file = pd.read_csv(vcf_file, sep="\t")
    # split the first line according to "]" and take the first element.
    # remove duplicate rows from the vcf_file - start from the second row
    vcf_file = vcf_file.iloc[:].drop_duplicates()
    # save the vcf_file as a csv file
    vcf_file.to_csv(f'{lib_plate}Pool_edited2checkk.csv', index=False)

    vcf_file.columns = vcf_file.columns.str.split("]").str[1]

    # Define the remove_suffix and replace_dash_with_underscore functions

    header_list = list(vcf_file.columns)
    for i in range(len(header_list)):
        header_list[i] = remove_suffix(header_list[i])
        header_list[i] = replace_dash_with_underscore(header_list[i])

    vcf_file.columns = header_list

    # Creating an empty dataframe to store the result
    result = pd.DataFrame(
        columns=["rs", "hom_ref", "het", "hom_alt", "total", "no_depth", "missing", "sample_ID", "Lib#", "Key"])

    # Initializing a list to store the missing rs values
    missing_rs = []

    # Looping through the rs_catalog dataframe
    for index, row in rs_catalog.iterrows():
        # Getting the current rs value
        rs = row["rs"]

        # Finding the row in the vcf file that matches the current rs value
        vcf_row = vcf_file[vcf_file["ID"] == rs]

        # If the rs value is not found in the vcf file, add it to the missing_rs list and continue to the next iteration
        if vcf_row.empty:
            missing_rs.append(rs)
            continue

        # Getting the sample values for the current rs value
        sample_values = vcf_row.iloc[:, 6:]

        # Counting the number of 0/0, 1/1, 1/0 or 0/1, and ./.' values
        num_00 = (sample_values == "0/0").sum().sum()
        num_11 = (sample_values == "1/1").sum().sum()
        num_10 = ((sample_values == "1/0") | (sample_values == "0/1")).sum().sum()
        num_total = num_00 + num_11 + num_10
        num_dot = (sample_values == "./.").sum().sum()

        # Finding the sample ID with homozygous alternate genotype
        sample_ID = ""
        for col in sample_values.columns:
            if "1/1" in sample_values[col].tolist():
                sample_ID = col
                break

        # If homozygous alternate genotype is not found, finding the sample ID with heterozygous genotype
        if not sample_ID:
            for col in sample_values.columns:
                if "1/0" in sample_values[col].tolist() or "0/1" in sample_values[col].tolist():
                    sample_ID = col
                    break

        # Adding the result to the result dataframe
        # Adding the result to the temporary DataFrame
        temp_df_data = {
            "rs": [rs],
            "hom_ref": [num_00],
            "het": [num_10],
            "hom_alt": [num_11],
            "total": [num_total],
            "no_depth": [num_dot],
            "missing": [None],
            "sample_ID": [sample_ID], "Lib#": [lib_plate]
            , "Key": [lib_plate + "_" + rs]
        }
        temp_df = pd.DataFrame(temp_df_data)
        result = pd.concat([result, temp_df], ignore_index=True)

        # Adding the missing rs values to the result dataframe
        for rs in missing_rs:
            # Looping through the missing rs values
            for rs in missing_rs:
                temp_missing_df_data = {
                    "rs": [rs],
                    "hom_ref": [None],
                    "het": [None],
                    "hom_alt": [None],
                    "total": [None],
                    "no_depth": [None],
                    "missing": ["missing"],
                    "sample_ID": [None],
                    "Key": [None]
                }

                # Creating a temporary DataFrame for missing values
                temp_missing_df = pd.DataFrame(temp_missing_df_data)

                # Concatenating the temporary DataFrame to the result DataFrame
                result = pd.concat([result, temp_missing_df], ignore_index=True)
        # remove duplicate rows
        result = result.drop_duplicates(subset=["rs", "hom_ref", "het", "hom_alt", "total", "no_depth", "missing",
                                                "sample_ID", "Lib#", "Key"], keep="first")
        # Writing the result dataframe to a csv file
        # take the sample_ID column and change each of them by using the extract function

        result.to_csv('rs_result.csv', index=False)

    # Sample_result
    # Creating an empty dataframe to store the result
    result = pd.DataFrame(
        columns=["Sample ID", "hom_ref", "het", "hom_alt", "total", "no_depth", "missing", "Lib#", "Key"])

    # Getting the sample IDs from the vcf file
    sample_IDs = vcf_file.columns[6:]

    # Looping through the sample IDs
    for sample_ID in sample_IDs:
        # Getting the sample values for the current sample ID
        sample_values = vcf_file[sample_ID]

        # Counting the number of 0/0, 1/1, 1/0 or 0/1, Total, and ./.' values
        num_00 = (sample_values == "0/0").sum()
        num_11 = (sample_values == "1/1").sum()
        num_10 = ((sample_values == "1/0") | (sample_values == "0/1")).sum()
        num_total = num_00 + num_11 + num_10
        num_dot = (sample_values == "./.").sum()

        # Adding the result to the result dataframe
        # Adding the result to the temporary DataFrame
        temp_sample_df_data = {
            "Sample ID": [sample_ID],
            "hom_ref": [num_00],
            "het": [num_10],
            "hom_alt": [num_11],
            "total": [num_total],
            "no_depth": [num_dot],
            "missing": [len(rs_catalog) - num_total],
            "Lib#": [lib_plate],
            "Key": [lib_plate + "_" + sample_ID],
        }
        # Creating a temporary DataFrame
        temp_sample_df = pd.DataFrame(temp_sample_df_data)

        # Concatenating the temporary DataFrame to the result_sample DataFrame
        result = pd.concat([result, temp_sample_df], ignore_index=True)

    # Saving the result dataframe to a csv file
    result.to_csv("sample_result.csv", index=False)


def copy_pool_vcf(lib_plate):
    # take the 4.vcf file from archive/LevHai and copy it to the working directory
    """
    The copy_wip_vcf function copies the 4.vcf file from archive/LevHai to the working directory
    and renames it according to user input.

    :return: The destination path of the copied file
    :doc-author: Sigal
    """
    # ask from user the file name
    source_file = f"archive\\LevHai\\{lib_plate}Pool.vcf"
    # destionation directory is the working directory
    destination_dir = os.getcwd()
    print(destination_dir)
    # copy the file to the working directory
    shutil.copy(source_file, destination_dir)


def processVcfFile(lib_plate):
    """
    The processVcfFile function takes in a library plate name as an argument.
    It then calls the run_sub_process function, which runs the subprocesses to create a vcf file for each sample on that library plate.
    The processVcfFile function then waits for all of those sub processes to finish before calling copy_pool_vcf, which combines all of the individual vcfs into one pool vcf file.
    Finally, it calls replace pipe with slash and reformat genotype functions to edit the pool vcf so that it can be used by other programs.

    :param lib_plate: Specify the library plate name
    :return: A file with the name of the library plate + &quot;pool_edited2
    :doc-author: Sigal
    """
    run_sub_process(lib_plate)
    # wait for the sub process to finish
    copy_pool_vcf(lib_plate)
    replace_pipe_with_slash(f"{lib_plate}Pool.vcf", f"{lib_plate}Pool_edited.vcf")
    reformat_genotype(f"{lib_plate}Pool_edited.vcf", f"{lib_plate}Pool_edited2.vcf")


def replace_dash_with_underscore(string: str) -> str:
    """
    The replace_dash_with_underscore function replaces all dashes in a string with underscores.

    :param string: str: Specify the type of data that is expected to be passed into the function
    :return: A string with all dashes replaced by underscores
    :doc-author: Sigal
    """
    return string.replace('-', '_')


def remove_suffix(string):
    """
    The remove_suffix function takes a string as input and returns the same string with its last underscore removed.
    If no underscore is found, the original string is returned.

    :param string: Pass the string to be modified into the function
    :return: A string with the last underscore and everything after it removed
    :doc-author: Sigal
    """
    last_underscore_index = string.rfind("_")
    if last_underscore_index != -1:  # check if underscore is found
        string = string[:last_underscore_index]
    return string


def remove_prefix(string):
    """
    The remove_prefix function takes a string as input and returns the same string with its first part, up to and
    including the last underscore, removed.
    If no underscore is found, the original string is returned.

    :param string: Pass the string to be modified into the function
    :return: A string with the first part up to the last underscore removed
    :doc-author: Sigal
    """
    last_underscore_index = string.rfind("_")
    if last_underscore_index != -1:  # check if underscore is found
        string = string[last_underscore_index + 1:]
    return string


def reformat_genotype(input_file, output_file):
    # open the input file -L14P1Pool_edited.vcf and create the output file -L14P1Pool_edited_reformat.vcf
    # now when there is * in the [5]ALT column, change the genotype in the whole row according to this dictoinary
    # 2/2 -> 1/1
    # 1/2 -> 1/1
    # 2/1 -> 1/1
    # 0/2 -> 0/1
    # 2/0 -> 1/0

    # Define the dictionary mapping original values to modified values
    replacements = {
        "2/2": "1/1",
        "1/2": "1/1",
        "2/1": "1/1",
        "0/2": "0/1",
        "2/0": "1/0",
    }
    # list to store the file contents in after modification
    modified_contents = []
    # Open the input file and read its contents
    with open(input_file, "r") as file:
        contents = file.read()
    # loop through the contents and check if the [5]ALT column contains a * character
    for line in contents.split("\n"):
        if len(line.split("\t")) >= 5 and line.split("\t")[4].__contains__("*"):
            # if * is found, replace the genotype with the modified value
            for original, modified in replacements.items():
                line = line.replace(original, modified)
            modified_contents.append(line)
            print(modified_contents[-1])
        else:
            # if * is not found, add the line to the modified contents list without modification
            modified_contents.append(line)
            print(modified_contents[-1])
    # write the modified contents to the output file
    with open(output_file, "w") as file:
        file.write("\n".join(modified_contents))


def reformat_SampleID(filename):
    # Open the file
    df = pd.read_csv(filename)

    # Check if the file name contains "sample_result"
    if "sample_result" in filename:
        # Check if "Sample ID" column exists in the DataFrame
        if "Sample ID" in df.columns:
            # Apply remove_prefix function only to non-null values in the "Sample ID" column
            df["Sample ID"] = df["Sample ID"].apply(lambda x: remove_prefix(x) if pd.notnull(x) else x)
            df["Key"] = df["Lib#"] + "_" + df["Sample ID"]
            df.to_csv(filename, index=False)
        else:
            print("No 'Sample ID' column found in the DataFrame.")
    else:
        # Check if "sample_ID" column exists in the DataFrame
        if "sample_ID" in df.columns:
            # Apply remove_prefix function only to non-null values in the "sample_ID" column
            df["sample_ID"] = df["sample_ID"].apply(lambda x: remove_prefix(x) if pd.notnull(x) else x)
            df.to_csv(filename, index=False)
        else:
            print("No 'sample_ID' column found in the DataFrame.")


if __name__ == "__main__":
    print("Welcome to the Sample Pool MyGenes Pipeline version 2.4")
    lib_plate = input("Enter the library and plate number of vcf file: ")
    input("make sure do have an archive folder with LevHai folder inside it, press enter to continue\n")
    processVcfFile(lib_plate)
    run(lib_plate)
    # press enter to exit
    reformat_SampleID("rs_result.csv")
    reformat_SampleID("sample_result.csv")
    input("Press Enter to exit")
