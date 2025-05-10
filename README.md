# LevHai mGenes - Parser version 2.4

## About
LevHai Parser version 2.4 processes initial variant analysis output files, converting them into comprehensive CSV reports with genotype data for each sample.

## Components

The LevHai Parser contains the following scripts:

1. **LevHai parser SampleSh.&Sums_v2.4.py** - Generates JSON files from the sample sheet and region summaries.
2. **LevHai_parser_VCF.py** - Processes VCF files and generates the final reports.
3. **LevHai parser rs_catalog_v2.4.py** - Generates a JSON file from the rs catalog file.
4. **edit_catalog.py** - Edits the rs catalog file.
5. **sample_poolv2.4.py** - Creates JSON files from the input files for the pool.

## Usage Guide

### 1. LevHai Parser SampleSh.&Sums

This script finds the following files in the working directory and parses them into JSON files:

**Input Files:**
- SampleSheet.csv (for the 2 plates)
- RegionReportSummary_P1_rs.txt
- RegionReportSummary_P2_rs.txt
- RegionReportSummary_P1_TXA.txt
- RegionReportSummary_P2_TXA.txt

**Output Files:**
- SampleSheet.json
- converted_data_P1_rs.json
- converted_data_P2_rs.json
- converted_data_P1_TXA.json
- converted_data_P2_TXA.json

### 2. LevHai Parser VCF

This script finds the following files in the working directory and parses them into JSON files:

**Input Files:**
- rs catalog.txt (updated catalog)
- VCF files from Partek server.

**Output Files:**
- rs_not_found.txt
- config/converted_data_VCF_P1_Edited.json
- config/converted_data_VCF_P2_Edited.json
- converted_data_VCF_final_P1.csv (in the working directory) - the final file for plate 1
- converted_data_VCF_final_P2.csv (in the working directory) - the final file for plate 2

### 3. LevHai Catalog

This script finds the file rs catalog.txt (updated catalog) and creates a JSON file from it.

**Input File:**
- rs catalog.txt (updated catalog)

**Output File:**
- config/rs_catalog.json

### 4. Sample Pool

This script receives the library and plate number of the VCF file from the user (should exist in the LevHai directory).

**Input:**
- VCF files - should be located in the LevHai folder.

**Output:**
- sample_result.csv
- rs_result.csv

## Notes

- Input files containing PHI (Protected Health Information) are not included in this repository.
- Please ensure the appropriate files are placed in the working directory before running each script.



## Requirements

- Python 3.x
- pandas
- tqdm
- functools
- json
- other dependencies as specified in the code files
