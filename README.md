# LevH mGenes - Parser version 2.4

This is a README file of LevHai Parser version 2.4. <br>
The LevHai Parser takes LevHai output files from Partek analysis as input, and creates 2 csv final report files as output. <br>
Each csv file is from a plate with the following columns:

1. CHROM
2. POS
3. ID
4. REF
5. ALT

The rest of the columns are the samples from the plate, each sample is given a genotype, per snp. <br>
## How to use - developer

The LevHai Parser contains 4 different scripts:

1. LevHai parser rs_catalog_v2.4.py - creates a json file from the rs catalog file. <br>
2. LevHai parser SampleSh.&Sums_v2.3.py - creates a json file from the sample sheet and region summaries. <br>
3. LevHai_parser_VCF.py - creates a json file from vcf files and procesess the final reports. <br>
4. sample_poolv2.3.py - create json files from the input files for the pool. <br>

## How to use - user <br>
There are 4 zip files in the ParserLevH2.4 folder, each unzipped file is a script that can be run by double clicking on it. <br>

The exe files are:
1. VCFparser2.4.exe - creates a json file from vcf files and procesess the final reports. <br>
2. SummariesSamples2.4.exe - creates a json file from the sample sheet and region summaries. <br>
3. SamplePool2.4.exe - create json files from the input files for the pool. <br>
4. rsCatalogParser2.3.exe - creates a json file from the rs catalog file. <br>


### running details: <br>
- LevHai Parser SampleSh.&Sums:

finds the following files in the working directory, and parses them into json files: <br>
<br>
input files:
1. SampleSheet.csv (for the 2 plates)
2. RegionReportSummary_P1_rs.txt
3. RegionReportSummary_P2_rs.txt
4. RegionReportSummary_P1_TXA.txt
5. RegionReportSummary_P2_TXA.txt
after running the script, the following files will be created in config folder (if not exist, will be created):
<br>
output files:
1. SampleSheet.json
2. converted_data_P1_rs.json
3. converted_data_P2_rs.json
4. converted_data_P1_TXA.json
5. converted_data_P2_TXA.json

- LevHai Parser VCF:
finds the following files in the working directory, and parses them into json files: <br>
<br>
input files:
1. rs catalog.txt (updated catalog)
2. config/converted_data_VCF_P1.json ( created by LevHai Parser SampleSh.&Sums.py)
3. config/converted_data_VCF_P2.json ( created by LevHai Parser SampleSh.&Sums.py)
after running the script, the following files will be created in config folder (if not exist, will be created):
<br>
output files:
1. rs_not_found.txt
2. config/converted_data_VCF_P1_Edited.json
3. config/converted_data_VCF_P2_Edited.json
4. converted_data_VCF_final_P1.csv (in the working directory) - the final file for plate 1
5. converted_data_VCF_final_P2.csv (in the working directory) - the final file for plate 2

- LevHai Catalog
finds the file rs catalog.txt (updated catalog) and creates a json file from it.<br>
the json file will be created in config folder (if not exist, will be created).
input file:
rs catalog.txt (updated catalog)
output file:
config/rs_catalog.json

- sample_poolv2.3
This script get from the user as input the library and plate number of the vcf file(should exist in the wsl LevHai directory).
After getting this it preprocesses the files, then runs the program to get rs_result and sample result of the initial vcf files.
input:
1. vcf files - should be located in the LevHai wsl folder.
output:
1. sample_result.csv
2. rs_result.csv

