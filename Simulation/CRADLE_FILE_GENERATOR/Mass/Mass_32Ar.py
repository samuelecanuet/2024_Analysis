import pandas as pd
import shutil
import numpy as np
import os
from oauth2client.service_account import ServiceAccountCredentials
import gspread

ONLY_FERMI = True

amu = 931494.10242  # keV/c^2
mass = 31997637.824 # µ amu

Mass = {}
hbar_eV = 6.582119569e-16  # eV.s
def ConvertionEnergytoTime(value):
    if (value == 0):
        return 0.
    return hbar_eV / (1e3*value) * np.log(2)  # Convert MeV to eV and then to seconds

def LABtoCM(value):
    return value*Mass["32Cl"]/Mass["31S"] 

#   -4   14   18   32 Ar    x   -2200.352       1.770      7700.0089     0.0553  B- -24190#       400#       31 997637.824       1.900
########################### 32Ar --------> other masses ###########################
for m in (list(range(-100, 110, 10))+[-500, 500, -1000, 1000]): # keV

    print(f"Generating files for mass {m} keV")
    if not os.path.exists(f"{m}"):
        os.mkdir(f"{m}")

    ######## AME DATA ###########
    shutil.copy(f"AMEdata.txt", f"{m}/AMEdata.txt")

    mass_keV = mass * 1e-6 * amu
    new_mass_keV = mass_keV + m
    new_mass = new_mass_keV * 1e6 / amu

    userline = f"  -4   14   18   32 Ar    x   -2200.352       1.770      7700.0089     0.0553  B- -24190#       400#       {new_mass}       1.900"
    # open file and replace the line with % with the user_line
    with open(f"{m}/AMEdata.txt", "r") as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if lines[i].startswith("  -4   14   18   32 Ar"):
                lines[i] = userline + "\n"
                break
    with open(f"{m}/AMEdata.txt", "w") as f:
        f.writelines(lines)
    f.close()

    ######## RADIATION ##########
    # google sheet modification
    sheet_name = "32Ar_AdoptedLevels_ENSDF_2024"
    sheet_id = "158YsHZMtVH4RmywMT9NbkNgxWDt-E9jrhj40iqikF8I"
    url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"

    # changing mass 32Ar and Mass Excess 32Ar in the google sheet

    # rewriting the google sheet with the new mass in cell H22
    scope = ["https://spreadsheets.google.com/feeds", "https://www.googleapis.com/auth/drive"]
    creds = ServiceAccountCredentials.from_json_keyfile_name("../gsheet_credentials.json", scope)
    client = gspread.authorize(creds)

    # Open the Google Sheet (by name or by ID)
    google_spreadsheet = client.open_by_key(sheet_id)     # or client.open(google_spreadsheet_name)
    google_spreadsheet_url = f"https://docs.google.com/spreadsheets/d/{google_spreadsheet.id}"

    # Select the worksheet by name
    worksheet = google_spreadsheet.worksheet(sheet_name)
    # Update the cell with the new mass
    worksheet.update_acell("I24", new_mass) # Mass 32Ar
    worksheet.update_acell("I23", -2200.352 + m) # Mass Excess 32Ar

    ## generating radiationdata file
    sheet_name = "32Ar_AdoptedLevels_ENSDF_2024"
    url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"
    df = pd.read_csv(url)
    ########################### GENREAL INFO ###########################
    Mass["32Ar"] = float(df["Energy Error.1"][df[df["31S levels Energy"] == "32Ar"].index[0]+3]) #µu
    Mass["32Cl"] = float(df["Energy Error.1"][df[df["31S levels Energy"] == "32Cl"].index[0]+3]) #µu
    Mass["31S"] = float(df["Energy Error.1"][df[df["31S levels Energy"] == "31S"].index[0]+3]) #µu
    Energy_31S_1 = df["31S levels Energy"][1]
    Energy_31S_2 = df["31S levels Energy"][2]

    # copy and past without base in the filename the file z18.a32_base
    shutil.copy("../32Ar/z18.a32_base", f"{m}/z18.a32")
    sheet_name = "32Ar_AllDeducedProtons"
    url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"
    df = pd.read_csv(url)
    p_betaplus = df["BR (absolute) p"][df["BR (absolute) p"].shape[0]-1]
    p_KEC = df["BR (absolute) ecKp"][df["BR (absolute) ecKp"].shape[0]-1]
    p_LEC = df["BR (absolute) ecLp"][df["BR (absolute) ecLp"].shape[0]-1]
    p_MEC = df["BR (absolute) ecMp"][df["BR (absolute) ecMp"].shape[0]-1]


    # Build z18.a32 only for beta p discrading gammas
    with open(f"{m}/z18.a32", "a") as f:
        # write total decay mode probability
        f.write(f"                               BetaPlus            0  -       {p_betaplus}\n")
        f.write(f"                               KshellEC            0  -       {p_KEC}\n")
        f.write(f"                               LshellEC            0  -       {p_LEC}\n")
        f.write(f"                               MshellEC            0  -       {p_MEC}\n")

        
        for index, row in df.iterrows():
            if (row["BR (absolute) p"] != 0):
                if (index == df.shape[0]-1):
                    # last row is the total decay mode probability
                    continue
                if (ONLY_FERMI == True and float(row['Levels Energy']) != 5046.3):
                    # only write the decay to the ground state
                    continue
                f.write(f"                               BetaPlus            {row['Levels Energy']}  -       {row['BR (absolute) βp']}     {row['Qβ']}\n")
        
        for index, row in df.iterrows():
            if (row["BR (absolute) p"] != 0):
                if (index == df.shape[0]-1):
                    # last row is the total decay mode probability
                    continue
                if (ONLY_FERMI == True and float(row['Levels Energy']) != 5046.3):
                    # only write the decay to the ground state
                    continue
                f.write(f"                               KshellEC            {row['Levels Energy']}  -       {row['BR (absolute) ecKp']}     {row['Qβ']}\n")

        for index, row in df.iterrows():
            if (row["BR (absolute) p"] != 0):
                if (index == df.shape[0]-1):
                    # last row is the total decay mode probability
                    continue
                if (ONLY_FERMI == True and float(row['Levels Energy']) != 5046.3):
                    # only write the decay to the ground state
                    continue
                f.write(f"                               LshellEC            {row['Levels Energy']}  -       {row['BR (absolute) ecLp']}     {row['Qβ']}\n")
        
        for index, row in df.iterrows():
            if (row["BR (absolute) p"] != 0):
                if (index == df.shape[0]-1):
                    # last row is the total decay mode probability
                    continue
                if (ONLY_FERMI == True and float(row['Levels Energy']) != 5046.3):
                    # only write the decay to the ground state
                    continue
                f.write(f"                               MshellEC            {row['Levels Energy']}  -       {row['BR (absolute) ecMp']}     {row['Qβ']}\n")

    digit = 2
    ########################### 32Cl* --------> 31S + p ###########################
    shutil.copy("../32Ar/z17.a32_base", f"{m}/z17.a32")

    counter_peak = 0
    list_peak = []
    with open(f"{m}/z17.a32", "a") as f:
        for index, row in df.iterrows(): 
            if (row["BR (absolute) p"] != 0):
                if (index == df.shape[0]-1):
                    # last row is the total decay mode probability
                    continue
                if (ONLY_FERMI == True and float(row['Levels Energy']) != 5046.3):
                #     # only write the decay to the ground state
                    continue
                f.write(f"P                {row['Levels Energy']}    -  {ConvertionEnergytoTime(0.0 if row['Width'] == '-' else float(row['Width']))}\n")
                if (row['BR'] > 0):
                # GS
                    f.write(f"                 Proton         0.0    -  {row['BR']}        {round(LABtoCM(row['Ground State Energy']), digit)}\n")

                if (ONLY_FERMI == True):
                    # only write the decay to the ground state
                    continue
                # First excited state
                if (row['BR.1'] > 0):
                    f.write(f"                 Proton         {Energy_31S_1}    -  {row['BR.1']}        {round(LABtoCM(row['1st Excited State Energy']), digit)}\n")

                # Second excited state
                if (row['BR.2'] > 0):
                    f.write(f"                 Proton         {Energy_31S_2}    -  {row['BR.2']}        {round(LABtoCM(row['2nd Excited State Energy']), digit)}\n")






    