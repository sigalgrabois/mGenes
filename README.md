# LevH mGenes - Parser version 2.4

## About the Project

LevHai Parser version 2.4 processes initial variant analysis output files, converting them into comprehensive CSV reports with genotype data for each sample.

## How to Use - Developer

The LevHai Parser contains the following scripts:

1. **LevHai parser SampleSh.&Sums_v2.4.py** - Generates JSON files from the sample sheet and region summaries.
2. **LevHai parser rs_catalog_v2.4.py** - Generates a JSON file from the rs catalog file.
3. **LevHai_parser_VCF.py** - Processes VCF files and generates the final reports.
4. **edit_catalog.py** - Edits the rs catalog file.
5. **sample_poolv2.4.py** - Creates JSON files from the input files for the pool.

### Running Details

#### LevHai Parser SampleSh.&Sums:

This script finds the following files in the working directory and parses them into JSON files:

**Input Files:**
1. SampleSheet.csv (for the 2 plates)
2. RegionReportSummary_P1_rs.txt
3. RegionReportSummary_P2_rs.txt
4. RegionReportSummary_P1_TXA.txt
5. RegionReportSummary_P2_TXA.txt

**Output Files:**
1. SampleSheet.json
2. converted_data_P1_rs.json
3. converted_data_P2_rs.json
4. converted_data_P1_TXA.json
5. converted_data_P2_TXA.json

#### LevHai Parser VCF:

This script finds the following files in the working directory and parses them into JSON files:

**Input Files:**
1. rs catalog.txt (updated catalog)
2. config/converted_data_VCF_P1.json (created by LevHai Parser SampleSh.&Sums_v2.4.py)
3. config/converted_data_VCF_P2.json (created by LevHai Parser SampleSh.&Sums_v2.4.py)

**Output Files:**
1. rs_not_found.txt
2. config/converted_data_VCF_P1_Edited.json
3. config/converted_data_VCF_P2_Edited.json
4. converted_data_VCF_final_P1.csv (in the working directory) - the final file for plate 1
5. converted_data_VCF_final_P2.csv (in the working directory) - the final file for plate 2

#### LevHai Catalog:

This script finds the file rs catalog.txt (updated catalog) and creates a JSON file from it.

**Input File:**
- rs catalog.txt (updated catalog)

**Output File:**
- config/rs_catalog.json

#### Sample Pool:

This script receives the library and plate number of the VCF file from the user (should exist in the LevHai directory).

**Input:**
1. VCF files - should be located in the LevHai folder.

**Output:**
1. sample_result.csv
2. rs_result.csv
