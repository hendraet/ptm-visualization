from collections import defaultdict
import pandas as pd
import numpy as np

file = "data/chris/cleavage_plot/PPc_COMPLETE_cutoff_0-05FDR_reformat_XX.csv"

def extract_numeric_value(value):
    return int(value[1][1:])

df = pd.read_csv(file)

new_data_2 = {}

# new_data_2['ID'] = df.iloc[:, 0].tolist()
# new_data_2['Neuropathology'] = df.iloc[:, 1].tolist()

ids = [np.nan]
ids.extend(df.iloc[:, 0].tolist())
new_data_2['ID'] = ids
neuropathologies = [np.nan]
neuropathologies.extend(df.iloc[:, 1].tolist())
new_data_2['Neuropathology'] = neuropathologies

new_data = {}
for col in df.columns[6:]:
    second_row_value = df.loc[0, col]    
    sights = second_row_value.split(';')
    for sight in sights:
        mod_pos = sight.split('@')[1]
        mod_type = sight.split('(')[1][0]
        mod_name = sight.split('(')[0].strip()
        mod_iso = col.split('.')[0]
        mod_origin = second_row_value.replace(" ", "")
        data = [mod_name, mod_type+mod_pos]
        data.extend(df[col][1:].tolist())
        new_data[mod_name+"."+mod_type+mod_pos+".("+mod_origin+")"+"."+mod_iso] = data

# new_data = {}
# for col in df.columns[6:]:
#     second_row_value = df.loc[0, col]    
#     sights = second_row_value.split(';')
#     for sight in sights:
#         mod_pos = sight.split('@')[1]
#         mod_type = sight.split('@')[0]
#         data = [mod_pos]
#         data.extend(df[col][1:].tolist())
#         new_data[mod_type+"."+mod_pos] = data

temp_dict = defaultdict(list)

sorted_dict = dict(sorted(new_data.items(), key=lambda item: extract_numeric_value(item[1])))
for key, value in sorted_dict.items():
    prefix = key.split('.(')[0]
    if prefix not in temp_dict:
        temp_dict[prefix] = []
    temp_dict[prefix].append(value)

result = defaultdict(list)
for prefix, lists in temp_dict.items():
    combined_list = lists[0][:2]  # Start with the first two common values
    for i in range(2, len(lists[0])):
        if any(lst[i] == '1' for lst in lists):
            combined_list.append(1)
        else:
            combined_list.append(0)
    result[prefix] = combined_list

new_data_2.update(result)
new_df = pd.DataFrame(new_data_2)
new_df.to_csv("data/chris/cleavage_plot/PPc_COMPLETE_cutoff_0-05FDR_reformat_XX_tarik.csv", index=False)
