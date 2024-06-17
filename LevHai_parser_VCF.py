import json
import os
import shutil
import pandas as pd
import re
import functools
from tqdm import tqdm
import subprocess


class Sample:
    """
    This class represents a Sample - each sample has its sample ID and a list of the rs in the sample.
    In addition, sample has a list of the rs that are not in the VCF file and a flag indicating if the depth is above 25.

    """

    def __init__(self, sample_ID, rs_list):
        self.Sample_ID = sample_ID
        self.rs_list = rs_list
        self.rs_missing = []
        self.flag = []


# this function inits the rs catalog - the rs catalog from the txt file - rs catalog.txt
# gets the path and returns the members of the rs catalog
def init_catalog(path_rs_catalog):
    """
    The init_catalog function takes in a path to the rs_catalog.txt file and returns lists in order to init the Catalog class:
        1) rs_list - list of all the rs numbers in the catalog
        2) CHROM_list - list of all chromosome values for each SNP (rs number)
        3) POS_list - list of all position values for each SNP (rs number)
        4) REF_list - list of reference allele values for each SNP (rs number). This is what is used as a baseline when comparing genotypes.  For example, if you have an A/A genotype, then your gen

    :param path_rs_catalog: Specify the path to the rs_catalog
    :return: A list of the rs in the first column
    :doc-author: Sigal
    """
    # read the data from the text file
    data = pd.read_table(path_rs_catalog, sep='\t')
    # create a list of the values in each column
    rs_list = [row[0] for index, row in data.iterrows()]
    CHROM_list = [row[1] for index, row in data.iterrows()]
    POS_list = [row[2] for index, row in data.iterrows()]
    REF_list = [row[3] for index, row in data.iterrows()]
    ALT_list = [row[4] for index, row in data.iterrows()]
    Assay_ID_list = [row[5] for index, row in data.iterrows()]
    Assay_Name_list = [row[6] for index, row in data.iterrows()]
    Amplicon_START_list = [row[7] for index, row in data.iterrows()]
    Amplicon_END_list = [row[8] for index, row in data.iterrows()]
    Amplicon_Len_list = [row[9] for index, row in data.iterrows()]
    return rs_list, CHROM_list, POS_list, REF_list, ALT_list, Assay_ID_list, Assay_Name_list, Amplicon_START_list, Amplicon_END_list, Amplicon_Len_list


# this class holds the rs catalog - the rs catalog from the txt file - rs catalog.txt
# uses the init_catalog function to init the members of the rs catalog
class Rs_catalog:
    # constructor
    def __init__(self, path_rs_catalog):
        self.rs, self.CHROM, self.POS, self.REF, self.ALT, self.Assay_ID, self.Assay_Name, self.Amplicon_START, self.Amplicon_END, self.Amplicon_Len = init_catalog(
            path_rs_catalog)

    def find_index(self, rs):
        """
        The find_index function takes a list of strings and returns the index of the string that is passed in as an argument.
            If there are multiple instances of the same string, it will return only one index.

        :param self: Refer to the instance of the class
        :param rs: Find the index of a specific value in the list
        :return: The index of the element in a list
        :doc-author: Sigal
        """
        return self.rs.index(rs)


def make_new_df(SampleName, SampleGT_list, SampleDP_list, df_original):
    """
    The make_new_df function takes the SampleName, SampleGT_list and SampleDP_list lists as input.
    It also takes the original data frame as an input.
    The function creates a new data frame with columns: CHROM, POS, ID, REF, ALT TYPE and 3 columns for each sample:
    SampleName (the name of the sample), SampleGT (the genotype value) and SampleDP (the depth value).
    The number of rows in this new data frame is equal to the number of samples * number of snps in original df.

    :param SampleName: Get the samplename from the original data frame
    :param SampleGT_list: Get the gt value from the original data frame
    :param SampleDP_list: Get the dp value from the original data frame
    :param df_original: Get the values from the original data frame
    :return: A data frame with the following columns:
    :doc-author: Trelent
    """
    # create a new data frame with the following columns
    df_new = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'TYPE', 'SampleName', 'SampleGT', 'SampleDP'])
    # snp list is the value from the column ID in the original data frame the 3rd column
    snp_list = df_original.iloc[:, 2].values.tolist()
    # make a list of the GT values in the original data frame - ['SampleGT'] run in a nested loop - the outer loop is
    # on the number of samples and the inner loop is on the number of snps each row in the new data frame is a snp in
    # a sample so the number of rows is the number of samples * the number of snps in the original data frame
    for i in range(len(SampleName)):
        for j in range(len(snp_list)):
            # make a new row in the new data frame
            # the GT and DP values are the values from the column SampleValue[i] in the original data frame
            GT = SampleGT_list[i]
            DP = SampleDP_list[i]
            # get the value from the column in GT place and in the j row
            df_row = pd.DataFrame(
                {'CHROM': df_original.iloc[j, 0], 'POS': df_original.iloc[j, 1], 'ID': df_original.iloc[j, 2],
                 'REF': df_original.iloc[j, 3], 'ALT': df_original.iloc[j, 4], 'TYPE': df_original.iloc[j, 5],
                 'SampleName': replace_dash_with_underscore(remove_suffix(SampleName[i])),
                 'SampleGT': df_original.loc[j, GT], 'SampleDP': df_original.loc[j, DP]}, index=[0])
            # add the new row to the new data frame using concat
            print(df_row)
            df_new = pd.concat([df_new, df_row], ignore_index=True)
    return df_new


# this function gets a path of the VCF file and a flag that indicates if the file is VCF or CSV and returns a new data frame - with 9 columns:
# ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'TYPE', 'SampleName', 'SampleGT', 'SampleDP']
# and converts the data frame to a json file
def read_convert_Data(file_path):
    """
    The read_convert_Data function reads the VCF file into a pandas DataFrame.
        It then creates two lists, one for SampleGT and one for SampleDP.
        The function extracts the sample name from the SampleGT list - it is the string before the first ':'.
        Then it removes from SampleName all of [] and any numbers inside of them.
         Finally, it converts this data frame to a json file.

    :param file_path: Get the file path of the vcf file
    :return: A data frame and a string (plate name)
    :doc-author: Sigal
    """
    # Read the VCF file into a pandas DataFrame
    global plate
    df_original = pd.read_table(file_path, header=0, delimiter="\t")
    # list holds each column name
    first_row = df_original.columns.values.tolist()
    # get the number of columns in the data frame
    # start the data frame from the 6th row
    first_row = first_row[6:]
    # split the list to two lists - one for SampleGT and one for SampleDP
    SampleGT = []
    SampleDP = []
    # split by the 2 last characters of the string - if it is GT then add it to SampleGT list and if it is DP then add it to SampleDP list
    for i in range(len(first_row)):
        if first_row[i][-2:] == 'GT':
            SampleGT.append(first_row[i])
        elif first_row[i][-2:] == 'DP':
            SampleDP.append(first_row[i])
    # extract the sample name from the SampleGT list - it is the string before the first ':'
    SampleName = []
    for i in range(len(SampleGT)):
        SampleName.append(SampleGT[i].split(':')[0])
    # remove from SampleName all the [] and the number inside the []
    for i in range(len(SampleName)):
        # use regex to find the pattern of [digit] and replace it with ''
        SampleName[i] = re.sub(r'\[\d+\]', '', SampleName[i])
        print(SampleName[i])

    new = make_new_df(SampleName, SampleGT, SampleDP, df_original)

    # check if the file is from plate 1 or plate 2 and convert the data frame to a json file

    if 'P1' in file_path:
        # take the new data frame and convert it to a json file - the first row is the header and the rest of the
        # rows are the values the keys of the json file are the header and the values are the rows
        # change the CHROM column to the string chr + the value in the column - for example chr1
        new['CHROM'] = 'chr' + new['CHROM'].astype(str)
        new.to_json("config/converted_data_VCF_P1.json", orient="records", indent=4)
        plate = 'P1'

    elif 'P2' in file_path:
        # take the new data frame and convert it to a json file - the first row is the header and the rest of the
        # rows are the values the keys of the json file are the header and the values are the rows
        new['CHROM'] = 'chr' + new['CHROM'].astype(str)
        new.to_json("config/converted_data_VCF_P2.json", orient="records", indent=4)
        plate = 'P2'
    return new, plate


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
    The remove_prefix function takes a string as input and returns the same string with its first underscore removed.
    If no underscore is found, the original string is returned.

    :param string: Pass the string to be modified into the function
    :return: A string with the first underscore and everything before it removed
    :doc-author: Sigal
    """
    first_underscore_index = string.find("_")
    if first_underscore_index != -1:  # check if underscore is found
        string = string[first_underscore_index + 1:]
    return string


# this function gets a path of the sample sheet and returns a list of the sample IDs
def create_sample_ID_list(vcf_df):
    """
    The create_sample_ID_list function takes in a vcf_df and returns a list of the sample IDs.
    The function first creates an empty list called sample_ID_list. The function then uses the set() method to create
    a unique set of values from the 'SampleName' column in vcf_df, which contains all of the sample IDs. The function
    then converts this unique set into a list using .tolist(). Finally, it returns this new list.

    :param vcf_df: Create a list of the sample ids from the vcf_df,
    :return: A list of the sample ids from the vcf_df
    :doc-author: Sigal
    """
    # create a list of the sample IDs from the vcf_df, the sample IDs are in the column 'SampleName'

    sample_ID_list = list(set(vcf_df['SampleName'].values.tolist()))
    print(len(sample_ID_list))
    # return a list of the sample IDs
    return sample_ID_list


def create_rs_list(sample_ID):
    """
    The create_rs_list function takes in a sample ID and returns a list of rs numbers that are associated with the
    sample ID. The function first checks if the sample ID has P2 or P3 in it, then reads the corresponding json file.
    The function then creates an empty set to store all of the rs numbers for that specific sample. It loops through each
    entry in data and modifies each entry's SampleName by removing its last 3 characters (the '_Px' part). If this modified
    SampleName matches with our inputted SampleID, we add its corresponding rs number to our set.

    :param sample_ID: Identify the sample id of the user
    :return: A list of rs numbers for the sample id
    :doc-author: Sigal
    """
    # read the json file
    # check if the sample ID has P1 of P2 in the string
    # if it does then read the json file of the P1 or P2
    if 'P1' in sample_ID:
        try:
            with open('config/converted_data_VCF_P1.json') as json_file:
                data = json.load(json_file)
        except FileNotFoundError:
            print("Error: JSON file not found")
            return []
        except json.JSONDecodeError:
            print("Error: invalid JSON format P1")
            return []
    elif 'P2' in sample_ID:
        try:
            with open('config/converted_data_VCF_P2.json') as json_file:
                data = json.load(json_file)
        except FileNotFoundError:
            print("Error: JSON file not found")
            return []
        except json.JSONDecodeError:
            print("Error: invalid JSON format P2")
            return []

    # create a set of the rs in the sample
    rs_list = set()
    for i in range(len(data)):
        # modify sample name to remove the last 3 characters
        sample_name = data[i]['SampleName']

        # check if modified sample name matches the input sample ID
        if sample_name == sample_ID:
            rs_list.add(data[i]['ID'])

    return rs_list


def open_VCF_file(file_path):
    # open the converted_data_VCF_P1.json file
    """
    The open_VCF_file function opens the converted_data_VCF_P2.json file and returns a dictionary of the data.

    :param file_path: Specify the file that is being opened
    :return: A list of dictionaries
    :doc-author: Sigal
    """
    with open(file_path, 'r') as json_file:
        data = json.load(json_file)
    return data


# this function gets a sample ID, rs, catalog flag and the data(json file content)
# and adds the rs to the data - depending on the flag it adds the rs as below or above the threshold (25)
def add_rs_to_data(sample_name, rs, catalog, flag, data, result):
    """
    The add_rs_to_data function takes in a sample name, rs number, catalog object, flag (True or False), data list and result.
    It finds the index of the rs in the catalog list. It then checks if the flag is True or False - if it is True then
    the rs is above threshold and adds it to data as a SNP with SampleGT = REF/REF else SampleGT = ./.. The function returns
    the updated data.

    :param sample_name: Specify the sample name
    :param rs: Find the index of the rs in the catalog list
    :param catalog: Find the index of the rs in the catalog list
    :param flag: Determine if the rs is above or below the threshold
    :param data: Store the data that will be written to the vcf file
    :param result: Store the result of the calculation
    :return: A list of dictionaries
    :doc-author: Sigal
    """
    # find the index of the rs in the catalog list
    index = catalog.find_index(rs)
    # check if the flag is True or False - if it is True then the rs is above the threshold
    if flag:
        sample_gt = f"{catalog.REF[index]}/{catalog.REF[index]}"
    else:
        sample_gt = "./."
    # add the rs to the data
    data.append({"CHROM": catalog.CHROM[index], "POS": catalog.POS[index], "ID": rs, "REF": catalog.REF[index],
                 "ALT": catalog.ALT[index], "TYPE": "SNP",
                 "SampleName": sample_name,
                 "SampleGT": sample_gt, "SampleDP": f"({result})"})
    return data


def check_rs_in_sample(sample, catalog_list):
    """
    The check_rs_in_sample function takes in a sample and a list of rs from the catalog.
    It returns a list of the rs that are not in the sample but are in the catalog.

    :param sample: Access the rs_list attribute of the sample object
    :param catalog_list: Compare the rs_list in the sample to
    :return: A list of the rs that are not in the sample but are in the catalog
    :doc-author: Sigal
    """
    # return a list of the rs that are not in the sample but are in the catalog
    rs_from_catalog = []
    for i in range(len(catalog_list)):
        if catalog_list[i] not in sample.rs_list:
            rs_from_catalog.append(catalog_list[i])
    return rs_from_catalog


@functools.lru_cache(maxsize=None)
def check_depth(sample, rs, plate):
    """
    The check_depth function takes a sample object, an rs value, and the plate number as arguments.
    It then opens the converted_data_rs json file for that plate and finds all of the rows where
    the first column is equal to one of the rs values in sample.rs_missing. It then applies two functions
    to each value in column 1: replace dash with underscore and remove suffix (see above). Then it finds
    the row where both columns match up with their respective arguments to this function (sample name &amp; rs). If a matching row was found, return its depth value; otherwise return -2.

    :param sample: Identify the sample
    :param rs: Find the row in the data frame where the rs value matches
    :param plate: Determine which converted_data_rs json file to open
    :return: The depth of a given sample and rs value
    :doc-author: Sigal
    """
    # open the converted data rs json file as a data frame - check which plate it is and open the correct file
    if plate == 'P1':
        df = pd.read_json('config/converted_data_rs_P1.json')
    elif plate == 'P2':
        df = pd.read_json('config/converted_data_rs_P2.json')
    # take from the df only the data with the rs in the first column that are in sample.rs_missing
    df = df[df.iloc[:, 0].isin(sample.rs_missing)]
    # then apply to all of the values the function remove_suffix
    # Find the row where the sample name and rs value match the given arguments to the function
    row = df[(df.iloc[:, 1] == sample.Sample_ID) & (df.iloc[:, 0] == rs)]
    # If a matching row was found, return the value of the depth column
    if len(row) > 0:
        # return the value of the depth column
        return int(row.iloc[0, 2])

    # Otherwise, return -1
    return -1


# this function gets a sample object and the catalog and adds the rs that are in the catalog but not in the sample
def add_missing_rs(sample, catalog):
    """
    The add_missing_rs function takes in a sample object and a catalog object as parameters.
    It opens the rs_summary_reign file and loops over the sample.missing_rs list to create a list of the rs content in the sample.
    The function then checks if each missing rs is above or below 25 depth, if it is above 25 depth it adds true to sample.flag
    and adds that missing rs to converted_data_VCF_P2 file with its corresponding information, otherwise it sets flag = false
    and does not add that missing rs to converted data VCF P2.

    :param sample: Access the sample
    :param catalog: Find the index of the rs in the catalog
    :return: The sample object
    :doc-author: Sigal
    """
    # # open the rs_summary_reign file
    # loop over the sample.missing_rs list and create a list of the rs content in the sample
    rs_content = []
    # loop over the sample.missing_rs list and check the depth of each rs
    for i in range(len(sample.rs_missing)):
        # check the depth of the rs
        depth = check_depth(sample.Sample_ID, sample.rs_missing[i])
        # if the depth is -1 then add the rs to the sample.missing_rs list
        index = catalog.find_index(sample.rs_missing[i])
        if depth == -1:
            print("rs not found in the summary reign file")
            sample.flag = -1
        # if the depth is not -1 then check if it is above the threshold of 25 and if it is above the threshold add true to sample.flag
        # and add the rs to the sample.missing_rs list
        else:
            if int(depth) > 25:
                sample.flag = True
                sample_gt = f"{catalog.REF[index]}/{catalog.REF[index]}"

            else:
                sample.flag = False
                sample_gt = "./."
            rs_content.append(
                {"CHROM": catalog.CHROM[index], "POS": catalog.POS[index], "ID": sample.rs_missing[i],
                 "REF": catalog.REF[index],
                 "ALT": catalog.ALT[index], "TYPE": "SNP", "SampleName": sample.Sample_ID,
                 "SampleGT": sample_gt, "SampleDP": f"({depth})"})
    # add the rs content to the converted_data_VCF_P1.json file and save it open the file with no with statement and with append mode
    with open('config/converted_data_VCF_P1.json', 'a') as outfile:
        json.dump(rs_content, outfile, indent=4)


def find_files():
    """
    The find_files function searches the current directory for files ending in .vcf and containing either 'P2' or 'P3'.
    It then returns a list of these files.

    :return: A list of file names
    :doc-author: Sigal
    """
    vcf_files = []

    # Find the text file in the current path location
    for file in os.listdir():
        # if file contains Pool continue
        if file.__contains__('Pool'):
            continue
        if file.endswith('wip.vcf') and 'P1' in file:
            vcf_files.append(file)
            print("found P1 vcf file")
            print(vcf_files)
            continue

        if file.endswith('wip.vcf') and 'P2' in file:
            vcf_files.append(file)
            print("found P2 vcf file")
            print(vcf_files)
            continue
    return vcf_files


def create_edited_vcf_file(data, plate):
    """
    The create_edited_vcf_file function creates a file converted_data_VCF_PX_Edited.json and saves the data list in it.
    The function uses on the data list - on each of the SampleName keys, replace dash with underscore and remove suffix.
    It also adds a new column to the converted_data_VC_Edited dataset as a concatenation of ID and SampleName.

    :param data: Pass the data list to the function
    :param plate: Determine which file to open and write the data list in
    :return: The data list
    :doc-author: Sigal
    """
    # Add a new column to the converted_data_VC_Edited dataset
    for item in data:
        item["Key"] = item["ID"] + "_" + item["SampleName"]

    # Create a file converted_data_VCF_PX_Edited.json and save the data list in it
    if plate == 'P1':
        with open('config/converted_data_VCF_P1_Edited.json', 'w') as outfile:
            # Use on the data list - on each of the SampleName keys the function replace_dash_with_underscore and
            # remove_suffix
            json.dump(data, outfile, indent=4)
    elif plate == 'P2':
        with open('config/converted_data_VCF_P2_Edited.json', 'w') as outfile:
            # Use on the data list - on each of the SampleName keys the function replace_dash_with_underscore and
            # remove_suffix
            json.dump(data, outfile, indent=4)

    return data


def create_final_result(plate):
    """
    The create_final_result function takes the converted_data_VCF file and pivots it to get the desired column format.
    It also removes duplicate IDs for each SampleName, as well as saves the final result to a CSV file.

    :param plate: Specify which plate to use
    :return: The pivoted table as a pandas dataframe
    :doc-author: Sigal
    """
    # Read the input file into a pandas DataFrame
    if plate == 'P1':
        input_df = pd.read_json("config/converted_data_VCF_P1_Edited.json")
    elif plate == 'P2':
        input_df = pd.read_json("config/converted_data_VCF_P2_Edited.json")

    # Check for duplicates of ID for each SampleName
    for sample in input_df["SampleName"].unique():
        sample_df = input_df[input_df["SampleName"] == sample]
        if sample_df["ID"].duplicated().any():
            print(f"Duplicate IDs found for SampleName '{sample}', removing duplicates...")
            input_df.drop(sample_df[sample_df["ID"].duplicated()].index, inplace=True)

    # Pivot the input table to get the specified column format with only SampleGT
    pivoted_table = input_df.pivot_table(index=["CHROM", "POS", "ID", "REF", "ALT", "TYPE"], columns="SampleName",
                                         values="SampleGT", aggfunc=lambda x: ';'.join(x))

    # Reset the index to get the pivoted columns as regular columns
    pivoted_table = pivoted_table.reset_index()
    # remove the column TYPE from the pivoted table
    pivoted_table = pivoted_table.drop(columns=["TYPE"])
    # run over the columns headers and if the column header has _ then use the function remove prefix to remove the prefix
    # from the column header
    for column in pivoted_table.columns:
        if "_" in column:
            pivoted_table = pivoted_table.rename(columns={column: remove_prefix(remove_prefix(column))})

    # Save the pivoted table to a CSV file with progress bar
    if plate == 'P1':
        with tqdm(total=len(pivoted_table.index), desc="Exporting to CSV") as pbar:
            pivoted_table.to_csv("converted_data_VCF_final_P1.csv", index=False, chunksize=1000)
            pbar.update(len(pivoted_table.index))
    elif plate == 'P2':
        with tqdm(total=len(pivoted_table.index), desc="Exporting to CSV") as pbar:
            pivoted_table.to_csv("converted_data_VCF_final_P2.csv", index=False, chunksize=1000)
            pbar.update(len(pivoted_table.index))

    # Print the pivoted table
    print(pivoted_table)


def run_sub_process(lib_plate):
    """
    The run_sub_process function runs the bcf tools commands in a subprocess.
    The function runs the commands on the converted_data_VCF_PX_Edited.json file.

    :return: The output of the subprocess
    :doc-author: Sigal
    """

    working_dir = input("Enter the working directory (wsl): ")
    rs_list = input("Enter the rs list file name: ")
    DP_value = input("Enter the DP value: (default is 35)")
    if DP_value == "":
        DP_value = "35"

    # cmd01 = f"bcftools filter -i 'ID=@/mnt/c/Users/User/LevHai/Filter_ToInclude_rs_list_v4.txt' -o '/mnt/c/Users/User/LevHai/{lib_plate}/Partek_{lib_plate}_Filetr.vcf' '/mnt/c/Users/User/LevHai/{lib_plate}/Partek_{lib_plate}.vcf'"
    # cmd02 = f"bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE[\t%TGT][\t%DP]\n' --print-header '/mnt/c/Users/User/LevHai/{lib_plate}/Partek_{lib_plate}_Filetr.vcf' -o '/mnt/c/Users/User/LevHai/{lib_plate}wip.vcf'"

    cmd01 = f"bcftools filter -i 'ID=@{working_dir}/{rs_list}' -o '{working_dir}/{lib_plate}/Partek_{lib_plate}_Filetr.vcf' '{working_dir}/{lib_plate}/Partek_{lib_plate}.vcf'"
    cmd02 = f"bcftools filter -S . -e 'FMT/DP<{DP_value}' '{working_dir}/{lib_plate}/Partek_{lib_plate}_Filetr.vcf' > '{working_dir}/{lib_plate}temp.vcf'"
    cmd03 = f"bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE[\t%TGT][\t%DP]\n' --print-header '{working_dir}/{lib_plate}temp.vcf' -o '{working_dir}/{lib_plate}wip.vcf'"

    process = subprocess.Popen(["wsl", "bash", "-c", f"cd {working_dir} && {cmd01} && {cmd02} && {cmd03}"],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # get the output of the subprocess
    output, error = process.communicate()
    # print the content of the output
    if error:
        print(error.decode("utf-8"))
    else:
        print(output.decode("utf-8"))
    copy_files_from_wsl(working_dir, "archive")


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


def copy_wip_vcf(lib_plate):
    # take the 4.vcf file from archive/LevHai and copy it to the working directory
    """
    The copy_wip_vcf function copies the 4.vcf file from archive/LevHai to the working directory
    and renames it according to user input.

    :return: The destination path of the copied file
    :doc-author: Sigal
    """
    # ask from user the file name
    source_file = f"archive\\LevHai\\{lib_plate}wip.vcf"
    # destionation directory is the working directory
    destination_dir = os.getcwd()
    print(destination_dir)
    # copy the file to the working directory
    shutil.copy(source_file, destination_dir)


def check_valid_input(lib_plate):
    if lib_plate[0] != 'L':
        print("Invalid input, please try again")
        lib_plate = input("Please enter the libary and plate you would like to run(e.g. L07P1): ")
        lib_plate = check_valid_input(lib_plate)
    for i in range(1, lib_plate.find('P')):
        if lib_plate[i].isdigit() == False:
            print("Invalid input, please try again")
            lib_plate = input("Please enter the libary and plate you would like to run(e.g. L07P1): ")
            lib_plate = check_valid_input(lib_plate)
    for i in range(lib_plate.find('P') + 1, len(lib_plate)):
        if lib_plate[i].isdigit() == False:
            print("Invalid input, please try again")
            lib_plate = input("Please enter the libary and plate you would like to run(e.g. L07P1): ")
            lib_plate = check_valid_input(lib_plate)
    return lib_plate


def replace_pipe_with_slash(input_file, output_file):
    """
    The replace_pipe_with_slash function takes two arguments:
        - input_file: The name of the file to read from.
        - output_file: The name of the file to write to.

    :param input_file: Specify the path to the input file
    :param output_file: Specify the name of the file to which
    :return: None
    :doc-author: Sigal
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


def edit_vcf_stars(input_file, output_file):
    try:
        # Open the input file and read its contents
        with open(input_file, "r") as file:
            contents = file.read()

        # Define the dictionary mapping original values to modified values
        replacements = {
            "C/*": "C/C",
            "G/*": "G/G",
            "T/*": "T/T",
            "A/*": "A/A",
            "*/C": "C/C",
            "*/G": "G/G",
            "*/T": "T/T",
            "*/A": "A/A"
        }

        # Replace each instance of the original with the value in the dictionary
        for key, value in replacements.items():
            contents = contents.replace(key, value)

        # Write the modified contents to the output file
        with open(output_file, "w") as file:
            # ensure that the file is empty before writing to it - delete all the content
            file.truncate(0)
            file.write(contents)

        print("File successfully processed.")

    except FileNotFoundError:
        print("The input file was not found.")

    except Exception as e:
        print("An error occurred during file processing:", e)


def runPlate():
    # input for the library and plate number of the first plate
    lib_plate = input("Please enter the libary and plate you would like to run(e.g. L07P1): ")
    # check if the input is valid
    lib_plate = check_valid_input(lib_plate)
    # run the sub process - run on wsl the bcftools commands
    run_sub_process(lib_plate)
    # copy the wip.vcf file to the working directory
    copy_wip_vcf(lib_plate)
    # input for the library and plate number of the second plate
    # input for the library and plate number of the second plate
    result = input("Would you like to edit the wip file? (y/n): ")
    while result != 'y' and result != 'n':
        print("Invalid input, please try again")
        result = input("Would you like to edit the wip file? (y/n): ")

    if result == 'y':  # if the answer is yes, replace the pipe with splash and the *
        replace_pipe_with_slash(f"{lib_plate}wip.vcf", f"{lib_plate}wip1.vcf")
        edit_vcf_stars(f"{lib_plate}wip1.vcf", f"{lib_plate}wip.vcf")
        print("The file was edited successfully")


if __name__ == '__main__':
    rs_catalog = Rs_catalog('rs catalog.txt')
    print("Welcome to files Parser - LevHai_parser_VCF.py version 2.4\n")
    print("This program converts vcf files by your choice - please follow the instructions\n")
    result = input("Would you like to run the program on 1 plate or on 2 plates?\n")
    if result == "1":
        runPlate()
    if result == "2":
        runPlate()
        runPlate()
    # look for the files in the working directory and create a list of the files
    vcf_files_found = find_files()
    for file in tqdm(vcf_files_found, desc="Processing VCF files"):
        vcf_table, plate = read_convert_Data(file)
        Sample_ID_list = create_sample_ID_list(vcf_table)
        sample_list = []
        # open the json file with the converted data according to the plate number
        if plate == 'P1':
            data = open_VCF_file('config/converted_data_VCF_P1.json')

        elif plate == 'P2':
            data = open_VCF_file('config/converted_data_VCF_P2.json')
        #  loop over the samples and create a rs list for each sample
        # create rs_missing list for each sample
        # for the missing rs - check the depth
        # if rs not found in the region summary file - the depth is -1 ( added to log file, and genotype set as ./.)
        # depth < 25 - genotype set as ./.
        # depth >= 25 - genotype set as homozygous reference
        # output is a json file with the converted data and csv with final results
        for i in tqdm(range(len(Sample_ID_list)), desc=f"Processing {os.path.basename(file)}"):
            rs_list = create_rs_list(Sample_ID_list[i])
            sample_list.append(Sample(Sample_ID_list[i], rs_list))
            rs_missing_list = check_rs_in_sample(sample_list[i], rs_catalog.rs)
            sample_list[i].rs_missing = rs_missing_list
            for rs in sample_list[i].rs_missing:
                result = check_depth(sample_list[i], rs, plate)
                if result == -1:
                    with open('rs_not_found.txt', 'a') as outfile:
                        outfile.write(f"{rs} - not found\n")
                    # open a file - for the rs that not found in the vcf file - save the rs and the sample id
                    print("rs not found in region summary file, added to the file as ./.")
                    data = add_rs_to_data(sample_list[i].Sample_ID, rs, rs_catalog, False, data, result)
                elif int(result) < 25:
                    print("depth is less then 25")
                    data = add_rs_to_data(sample_list[i].Sample_ID, rs, rs_catalog, False, data, result)
                else:
                    print("depth is more then 25")
                    data = add_rs_to_data(sample_list[i].Sample_ID, rs, rs_catalog, True, data, result)
            create_edited_vcf_file(data, plate)
            create_final_result(plate)
    # press any key to exit
    input("Press any key to exit")
