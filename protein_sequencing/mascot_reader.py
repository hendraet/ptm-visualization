import pandas as pd

# Read CSV file into a DataFrame
input_data = pd.read_csv("data/mascot/AD_Insol_01_F045338.csv", skiprows=63)

input_data = input_data[input_data['pep_var_mod'].notna() & ~input_data['pep_var_mod'].str.contains('Oxidation', na=False)]

for index, row in input_data.iterrows():
    print(row)
