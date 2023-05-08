# LevH - Parser version 1.0

This is a README file of LevHai Parser version 1.0. <br>
The LevHai Parser takes LevHai output files from running and analyzing 2 plates and parses them into
readable formats. <br>
The final output is 2 csv files, each is from a plate with the following columns:

1. Sample Name
2. rs#
3. Chromosome
4. Position
5. Reference Allele
6. Alternative Allele
7. Depth
8. Allele Frequency
9. Genotype

## How to use

The LevHai Parser contains 3 different scripts:

1. LevHai Parser SampleSh.&Sums.py
2. LevHai Parser VCF.py
3. LevHai Catalog.py

### LevHai Parser SampleSh.&Sums.py

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

### LevHai Parser VCF.py
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

### LevHai Catalog.py
finds the file rs catalog.txt (updated catalog) and creates a json file from it.<br>
the json file will be created in config folder (if not exist, will be created).
input file:
rs catalog.txt (updated catalog)
output file:
config/rs_catalog.json



## change log

[sigal 01/02/2023 12:43]
LevHai Parser v1.1.1 - fixed line 56 "SampleName[i] = SampleName[i]" deleted this part, the current code only remove the
parentheses and keep the number within, result in wrong sample number.

- now working and deleting the parentheses and the number within.

[chen 24/01/2023 16:17]
LevHai Parser v1.1.1.1 - new code to include catalogue of rs#, every execution need to repopulate the data.
LevHai Parser v1.1.1 - line 56 "SampleName[i] = SampleName[i]" deleted this part, the current code only remove the
parentheses and keep the number within, result in wrong sample number.
LevHai Parser v1.1 - update to include script_log.txt file (need to include error msgs).

- The code for sctipr_log was not saved...

[sigal 26/02/2023 14:00]
LevHai Parser v1.1.1 - added the sample class in order to add the missing rs from the catalog

[sigal 09/03/2023 14:00]
LevHai Parser v1.1.1 - renaming the file to LevHai Parser VCF.py, fixing the running time of the check depth function
added catalog and sample class to the code
LevHai Parser v1.1 - renaming the file to LevHai Parser SampleSh.&Sums.py
LevHai Parser v1.1.1.1 - renaming the file to LevHai Catalog.py

[sigal 11/03/2023 14:00]
LevHai Parser v1.1.1 - renaming the file to LevHai Parser VCF.py, fixing the running time of the check depth function
added catalog and sample class to the code
LevHai Parser v1.1 - renaming the file to LevHai Parser SampleSh.&Sums.py
LevHai Parser v1.1.1.1 - renaming the file to LevHai Catalog.py

[sigal 01/04/2023 13:00]
LevHai parser SampleSh.&Sums.py -  added a reformat function to create the json file for the VCF parser - as list of jsons
                                -  use the remove suffix function to remove the suffix from the sample name
                                -  use function to replace the '-' with '_' in the sample name
LevHai parser rs_catalog.py     - added at the end the option to press any key to exit the program




