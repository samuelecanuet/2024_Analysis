import pandas as pd
import shutil
import numpy as np
import gspread
from google.oauth2.service_account import Credentials

Mass = {}
sheet_id = "1UPYVJ3ssUV1XORsqBXmlmHLym6_y6D54FvzY5YOheK8"
sheet_name = "33Ar_AdoptedLevels_ENSDF_2024"


hbar_eV = 6.582119569e-16  # eV.s

url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"
df = pd.read_csv(url)

# with pd.option_context('display.max_rows', None, 'display.max_columns', None):
#     print(df)

def ConvertionEnergytoTime(value):
    if (value == 0):
        return 0.
    return hbar_eV / (1e3*value) * np.log(2)  # Convert MeV to eV and then to seconds

def LABtoCM(value):
    return value*Mass["33Cl"]/Mass["32S"] 

########################### GENREAL INFO ###########################
Mass["33Ar"] = float(df["Energy Error.1"][df[df["32S levels Energy"] == "33Ar"].index[0]+3]) #µu
Mass["33Cl"] = float(df["Energy Error.1"][df[df["32S levels Energy"] == "33Cl"].index[0]+3]) #µu
Mass["32S"] = float(df["Energy Error.1"][df[df["32S levels Energy"] == "32S"].index[0]+3]) #µu
Energy_32S_1 = df["32S levels Energy"][1]
Energy_32S_2 = df["32S levels Energy"][2]

####################################################################

sheet_name = "33Ar_AllDeducedProtons"
url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"
df = pd.read_csv(url)

########################### 33Ar --------> 33Cl + e + nu ###########################

# copy and past without base in the filename the file z18.a32_base
shutil.copy("z18.a33_base", "z18.a33_ENSDF")

p_betaplus = df["BR (absolute) p"][df["BR (absolute) p"].shape[0]-1]
# p_KEC = df["BR (absolute) ecKp"][df["BR (absolute) ecKp"].shape[0]-1]
# p_LEC = df["BR (absolute) ecLp"][df["BR (absolute) ecLp"].shape[0]-1]
# p_MEC = df["BR (absolute) ecMp"][df["BR (absolute) ecMp"].shape[0]-1]


# Build z18.a32 only for beta p discrading gammas
with open("z18.a33_ENSDF", "a") as f:
    # write total decay mode probability
    f.write(f"                               BetaPlus            0  -       {p_betaplus}\n")
    # f.write(f"                               KshellEC            0  -       {p_KEC}\n")
    # f.write(f"                               LshellEC            0  -       {p_LEC}\n")
    # f.write(f"                               MshellEC            0  -       {p_MEC}\n")

    
    for index, row in df.iterrows():
        if (row["BR (absolute) p"] != 0):
            if (index == df.shape[0]-1):
                # last row is the total decay mode probability
                continue
            f.write(f"                               BetaPlus            {row['Levels Energy']}  -       {row['BR (absolute) βp']}     {row['Qβ']}\n")
    
    # for index, row in df.iterrows():
    #     if (row["BR (absolute) p"] != 0):
    #         if (index == df.shape[0]-1):
    #             # last row is the total decay mode probability
    #             continue
    #         f.write(f"                               KshellEC            {row['Levels Energy']}  -       {row['BR (absolute) ecKp']}     {row['Qβ']}\n")

    # for index, row in df.iterrows():
    #     if (row["BR (absolute) p"] != 0):
    #         if (index == df.shape[0]-1):
    #             # last row is the total decay mode probability
    #             continue
    #         f.write(f"                               LshellEC            {row['Levels Energy']}  -       {row['BR (absolute) ecLp']}     {row['Qβ']}\n")
    
    # for index, row in df.iterrows():
    #     if (row["BR (absolute) p"] != 0):
    #         if (index == df.shape[0]-1):
    #             # last row is the total decay mode probability
    #             continue
    #         f.write(f"                               MshellEC            {row['Levels Energy']}  -       {row['BR (absolute) ecMp']}     {row['Qβ']}\n")




digit = 2
########################### 33Cl* --------> 32S + p ###########################
shutil.copy("z17.a33_base", "z17.a33_ENSDF")

counter_peak = 0
list_peak = []
with open("z17.a33_ENSDF", "a") as f:
    for index, row in df.iterrows(): 
        if (row["BR (absolute) p"] != 0):
            if (index == df.shape[0]-1):
                # last row is the total decay mode probability
                continue

            f.write(f"P                {row['Levels Energy']}    -  {ConvertionEnergytoTime(0.0 if str(float(row['Width'])) == 'nan' else float(row['Width']))}\n")

            if (row['BR'] > 0):
            # GS
                f.write(f"                 Proton         0.0    -  {row['BR']}        {round(LABtoCM(row['Ground State Energy']), digit)}\n")
                list_peak.append([counter_peak, row['Ground State Energy'], row['Energy Error'], row['BR'], row['BR Error']])

            # First excited state
            if (row['BR.1'] > 0):
                f.write(f"                 Proton         {Energy_32S_1}    -  {row['BR.1']}        {round(LABtoCM(row['1st Excited State Energy']), digit)}\n")
                list_peak.append([counter_peak, row['1st Excited State Energy'], row['Energy Error.1'], row['BR.1'], row['BR Error.1']])

            # Second excited state
            if (row['BR.2'] > 0):
                f.write(f"                 Proton         {Energy_32S_2}    -  {row['BR.2']}        {round(LABtoCM(row['2nd Excited State Energy']), digit)}\n")
                list_peak.append([counter_peak, row['2nd Excited State Energy'], row['Energy Error.2'], row['BR.2'], row['BR Error.2']])

# ordering the map with the proton energy
list_peak.sort(key=lambda x: x[1])

with open("33Ar_protons.txt", "w") as f:
    for i, peak in enumerate(list_peak):
        f.write(f"{i+1} \t {peak[1]} \t {peak[2]} \t {peak[3]} \t {peak[4]}\n")


