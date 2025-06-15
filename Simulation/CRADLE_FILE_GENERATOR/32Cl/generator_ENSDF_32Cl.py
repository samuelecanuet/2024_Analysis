import pandas as pd
import shutil
import numpy as np
import gspread
from google.oauth2.service_account import Credentials

Mass = {}
sheet_id = "1UPYVJ3ssUV1XORsqBXmlmHLym6_y6D54FvzY5YOheK8"
sheet_name = "32Cl_AdoptedLevels_ENSDF_2024"


hbar_eV = 6.582119569e-16  # eV.s

url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"
df = pd.read_csv(url)

# with pd.option_context('display.max_rows', None, 'display.max_columns', None):
#     print(df)

def ConvertionEnergytoTime(value):
    if (value == 0):
        return 0.
    return hbar_eV / (1e3*value) * np.log(2)  # Convert MeV to eV and then to seconds

def LABtoCM(value, type):
    if (type == "p"):
        return value*Mass["32S"]/Mass["31P"] 
    elif (type == "a"):
        return value*Mass["32S"]/Mass["28Si"]

########################### GENREAL INFO ###########################
Mass["32Cl"] = float(df["Energy Error.1"][df[df["31P levels Energy"] == "32Cl"].index[0]+3]) #µu
Mass["32S"] = float(df["Energy Error.1"][df[df["31P levels Energy"] == "32S"].index[0]+3]) #µu
Mass["31P"] = float(df["Energy Error.1"][df[df["31P levels Energy"] == "31P"].index[0]+3]) #µu
Mass["28Si"] = float(df["Energy Error.1"][df[df["31P levels Energy"] == "28Si"].index[0]+3]) #µu
####################################################################

sheet_name = "32Cl_DeducedProtonsAlpha"
url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"
df = pd.read_csv(url)
########################### 32Cl --------> 32S + e + nu ###########################

# copy and past without base in the filename the file z17.a32_base
shutil.copy("z17.a32_base", "z17.a32_ENSDF")

p_betaplus = df["BR (absolute)"][df["BR (absolute)"].shape[0]-1]
# p_KEC = df["BR (absolute) ecKp"][df["BR (absolute) ecKp"].shape[0]-1]
# p_LEC = df["BR (absolute) ecLp"][df["BR (absolute) ecLp"].shape[0]-1]
# p_MEC = df["BR (absolute) ecMp"][df["BR (absolute) ecMp"].shape[0]-1]


# Build z17.a32 only for beta p discrading gammas
with open("z17.a32_ENSDF", "a") as f:
    # write total decay mode probability
    f.write(f"                               BetaPlus            0  -       {p_betaplus}\n")
    
    for index, row in df.iterrows():
        if (row["BR (absolute)"] != 0):
            if (index == df.shape[0]-1):
                # last row is the total decay mode probability
                continue
            f.write(f"                               BetaPlus            {row['Levels Energy']}  -       {row['BR (absolute)']}     {row['Qβ']}\n")




########################### 32S* --------> 31P + p ###########################
shutil.copy("z16.a32_base", "z16.a32_ENSDF")

counter_peak = 0
list_peak = []
with open("z16.a32_ENSDF", "a") as f:
    for index, row in df.iterrows(): 
        if (row["BR (absolute)"] != 0):
            if (index == df.shape[0]-1):
                # last row is the total decay mode probability
                continue
            if (row['BR'] > 0):
                f.write(f"P                {row['Levels Energy']}    -  {ConvertionEnergytoTime(0.0 if row['Width'] == '-' else float(row['Width']))}\n")

            # GS
                e = LABtoCM(row['Proton Energy'], "p")
                f.write(f"                 Proton         0.0    -  {row['BR']}        {e}\n")
                list_peak.append([counter_peak, row['Proton Energy'], row['Energy Error']])


# ordering the map with the proton energy
# list_peak.sort(key=lambda x: x[1])

# with open("32Cl_protons.txt", "w") as f:
#     for i, peak in enumerate(list_peak):
#         f.write(f"{i+1} \t {peak[1]} \t {peak[2]}\n")

########################### 32S* --------> 28Si + a ###########################
counter_peak = 0
# list_peak = []
print(df)
with open("z16.a32_ENSDF", "a") as f:
    for index, row in df.iterrows(): 
        if (row["BR (absolute)"] != 0):
            if (index == df.shape[0]-1):
                # last row is the total decay mode probability
                continue
            if (row['BR.1'] > 0):
                f.write(f"P                {row['Levels Energy']}    -  {ConvertionEnergytoTime(0.0 if row['Width'] == '-' else float(row['Width']))}\n")

            # GS
                e = LABtoCM(row['Alpha Energy'], "a")
                f.write(f"                 Alpha         0.0    -  {row['BR.1']}        {e}\n")
                list_peak.append([counter_peak, row['Alpha Energy'], row['Energy Error.1']])


# ordering the map with the proton energy
list_peak.sort(key=lambda x: x[1])

with open("32Cl_decalyed.txt", "w") as f:
    for i, peak in enumerate(list_peak):
        f.write(f"{i+1} \t {peak[1]} \t {peak[2]}\n")