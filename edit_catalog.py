import pandas as pd
from decimal import Decimal as decimal


# Function to find the minimum value while ignoring 'inf' and converting strings to Decimal
def custom_min(values):
    filtered_values = [decimal(val) for val in values if val != 'inf']
    if filtered_values:
        return min(filtered_values)
    else:
        return decimal('inf')

def scientific_to_float(scientific_notation):
    try:
        # Use Decimal to convert the scientific notation string to a float
        float_value = decimal(scientific_notation)
        print(f"Converted {scientific_notation} to {float_value}")
        return float_value
    except ValueError:
        # Handle the case where the input is not a valid scientific notation
        return None

# Define the file path
file_path = "25 Sequences Annotations.csv"

# Read the CSV file into a DataFrame
df = pd.read_csv(file_path)

# print the first 5 rows of the DataFrame
print(df.head())

# create a list of names from the df
names_list = df['Name'].tolist()

# Take the FREQ column and make a list out of this
freq_list = df['FREQ'].tolist()

# Initialize an empty list to store the extracted numeric values
numeric_values = []

# Extract numeric values from each element in freq_list using regular expressions
source_list = []
frequencies_list = []
for item in freq_list:
    # split each item by ":" the first item is the source of the data, the second is the frequency
    source, frequency = item.split(":")
    # append the source to the source list
    source_list.append(source)
    # append the frequency to the frequencies list
    frequencies_list.append(frequency)

# Split frequencies_list into lists of numbers (split by ",")
frequencies_list = [item.split(",") for item in frequencies_list]

# Iterate through each list in frequencies_list and make a list of the indexes of the strings containing scientific notation
scientific_notation_indexes = []
for i in range(len(frequencies_list)):
    for j in range(len(frequencies_list[i])):
        if "e" in frequencies_list[i][j]:
            scientific_notation_indexes.append((i, j))

# use function scientific_to_float to convert the scientific notation to float
for i, j in scientific_notation_indexes:
    frequencies_list[i][j] = scientific_to_float(frequencies_list[i][j])

# Iterate through each list in frequencies_list and convert each string to a float, check if the type is decimal and if so, continue and dont convert
for i in range(len(frequencies_list)):
    for j in range(len(frequencies_list[i])):
        if type(frequencies_list[i][j]) == str and not ".":
            frequencies_list[i][j] = float(frequencies_list[i][j])

# replace the "." with value of infinty in the frequencies list so that the minor frequency will be the only value
for i in range(len(frequencies_list)):
    for j in range(len(frequencies_list[i])):
        if frequencies_list[i][j] == ".":
            frequencies_list[i][j] = float("inf")

print(frequencies_list)

# Find the minimum for each sublist in data
minor_frequencies = [custom_min(sublist) for sublist in frequencies_list]

print(minor_frequencies)

# create a txt file in the format of name(from name list) + "_" + source + "_" + minor_frequency all in one file one after the other
with open("minor_frequencies.txt", "w") as f:
    f.write("MAF\n")
    for i in range(len(names_list)):
        f.write(names_list[i] + "_" + source_list[i] + "_" + str(minor_frequencies[i]) + "\n")



