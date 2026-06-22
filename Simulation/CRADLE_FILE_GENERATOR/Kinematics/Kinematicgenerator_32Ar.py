import pandas as pd
import shutil
import numpy as np
import os
from oauth2client.service_account import ServiceAccountCredentials
import gspread

ONLY_FERMI = True

amu = 931494.10242  # keV/c^2
mass_Ar = 31997637.824 # µ amu
mass_Cl = 31985684.605 # µ amu
mass_S = 30979557.002 # µ amu

sep_Ar = -2200.352 # keV
sep_Cl = -13334.706 # keV
sep_S = -19042.531 # keV

Mass = {}
hbar_eV = 6.582119569e-16  # eV.s
def ConvertionEnergytoTime(value):
    if (value == 0):
        return 0.
    return hbar_eV / (1e3*value) * np.log(2)  # Convert MeV to eV and then to seconds

def LABtoCM(value):
    return value*Mass["32Cl"]/Mass["31S"] 

############# INITIALIZATION USING GOOGLESHEET TABLE ################
sigma_Ar = 1.9
sigma_Cl = 0.603
sigma_S = 0.246
sigma_Ex_Ar = 0.3
DATA = [

]
#reading the google sheet with the new masses
sheet_name = "sheet"
sheet_id = "1Qk-XyQm8JEIIqQkvgnSpS0mGkMPsbPCtC6rTGDcpofw"
scope = ["https://spreadsheets.google.com/feeds", "https://www.googleapis.com/auth/drive"]
creds = ServiceAccountCredentials.from_json_keyfile_name("../gsheet_credentials.json", scope)
client = gspread.authorize(creds)
google_spreadsheet = client.open_by_key(sheet_id)     # or client.open(google_spreadsheet_name)
google_spreadsheet_url = f"https://docs.google.com/spreadsheets/d/{google_spreadsheet.id}"
worksheet = google_spreadsheet.worksheet(sheet_name)

for row in worksheet.get_all_values()[1:]:  # Skip header
    name = row[0]
    sigma_sign_Ar = float(row[1])
    sigma_sign_Cl = float(row[2])
    sigma_sign_S = float(row[3])
    sigma_sign_Ex = float(row[4])
    DATA.append({"Name": name, "Ar": sigma_sign_Ar*sigma_Ar, "Cl": sigma_sign_Cl*sigma_Cl, "S": sigma_sign_S*sigma_S, "Ex": sigma_sign_Ex*sigma_Ex_Ar})


########################### 32Ar --------> other masses ###########################
for name, sigma_Ar, sigma_Cl, sigma_S, sigma_Ex in [(d["Name"], d["Ar"], d["Cl"], d["S"], d["Ex"]) for d in DATA]:
    print(f"Generating files for case {name} keV")

    # if "Ar" in name or "S" in name or "Ex" in name:
    #     continue    

    if not os.path.exists(f"{name}"):
        os.mkdir(f"{name}")

    ######## AME DATA ###########
    shutil.copy(f"AMEdata.txt", f"{name}/AMEdata.txt")

    new_mass_Ar = mass_Ar + sigma_Ar
    new_sep_Ar = sep_Ar + sigma_Ar*amu * 1e-6

    new_mass_Cl = mass_Cl + sigma_Cl
    new_sep_Cl = sep_Cl + sigma_Cl*amu * 1e-6

    new_mass_S = mass_S + sigma_S
    new_sep_S = sep_S + sigma_S*amu * 1e-6

    new_Ex = 5046.3 + sigma_Ex

    userlineAr = f"  -4   14   18   32 Ar    x   -2200.352       1.770      7700.0089     0.0553  B- -24190#       400#       {new_mass_Ar}       1.900"
    userlineCl = f"  -2   15   17   32 Cl       -13334.706       0.562      8072.4058     0.0176  B- -11134.3536     1.8568   {new_mass_Cl}       0.603"
    userlineS = f"  -1   15   16   31 S        -19042.531       0.229      8281.8013     0.0074  B- -12007.9790     3.4541   {new_mass_S}       0.246"
    # open file and replace the line with % with the user_line
    with open(f"{name}/AMEdata.txt", "r") as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if lines[i].startswith("  -4   14   18   32 Ar"):
                lines[i] = userlineAr + "\n"
            if lines[i].startswith("  -2   15   17   32 Cl"):
                lines[i] = userlineCl + "\n"
            if lines[i].startswith("  -1   15   16   31 S"):
                lines[i] = userlineS + "\n"
    with open(f"{name}/AMEdata.txt", "w") as f:
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
    worksheet.update_acell("I24", new_mass_Ar) # Mass 32Ar
    worksheet.update_acell("I23", new_sep_Ar) # Mass Excess 32Ar
    worksheet.update_acell("I12", new_mass_Cl) # Mass 32Cl
    worksheet.update_acell("I11", new_sep_Cl) # Mass Excess 32Cl
    worksheet.update_acell("I18", new_mass_S) # Mass 31S
    worksheet.update_acell("I17", new_sep_S) # Mass Excess 31S
    worksheet.update_acell("B29", new_Ex) # Ex 32Cl

    ## generating radiationdata file
    sheet_name = "32Ar_AdoptedLevels_ENSDF_2024"
    url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"
    df = pd.read_csv(url)
    ########################### GENREAL INFO ###########################
    Mass["32Ar"] = new_mass_Ar  
    Mass["32Cl"] = new_mass_Cl
    Mass["31S"] = new_mass_S
    Energy_31S_1 = df["31S levels Energy"][1]
    Energy_31S_2 = df["31S levels Energy"][2]

    # copy and past without base in the filename the file z18.a32_base
    shutil.copy("../32Ar/z18.a32_base", f"{name}/z18.a32")
    sheet_name = "32Ar_AllDeducedProtons"
    url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"
    df = pd.read_csv(url)
    p_betaplus = df["BR (absolute) p"][df["BR (absolute) p"].shape[0]-1]
    p_KEC = df["BR (absolute) ecKp"][df["BR (absolute) ecKp"].shape[0]-1]
    p_LEC = df["BR (absolute) ecLp"][df["BR (absolute) ecLp"].shape[0]-1]
    p_MEC = df["BR (absolute) ecMp"][df["BR (absolute) ecMp"].shape[0]-1]


    # Build z18.a32 only for beta p discrading gammas
    with open(f"{name}/z18.a32", "a") as f:
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
                if (ONLY_FERMI == True and float(row['Levels Energy']) != new_Ex):
                    # only write the decay to the ground state
                    continue
                f.write(f"                               BetaPlus            {row['Levels Energy']}  -       {row['BR (absolute) βp']}     {row['Qβ']}\n")
        
        for index, row in df.iterrows():
            if (row["BR (absolute) p"] != 0):
                if (index == df.shape[0]-1):
                    # last row is the total decay mode probability
                    continue
                if (ONLY_FERMI == True and float(row['Levels Energy']) != new_Ex):
                    # only write the decay to the ground state
                    continue
                f.write(f"                               KshellEC            {row['Levels Energy']}  -       {row['BR (absolute) ecKp']}     {row['Qβ']}\n")

        for index, row in df.iterrows():
            if (row["BR (absolute) p"] != 0):
                if (index == df.shape[0]-1):
                    # last row is the total decay mode probability
                    continue
                if (ONLY_FERMI == True and float(row['Levels Energy']) != new_Ex):
                    # only write the decay to the ground state
                    continue
                f.write(f"                               LshellEC            {row['Levels Energy']}  -       {row['BR (absolute) ecLp']}     {row['Qβ']}\n")
        
        for index, row in df.iterrows():
            if (row["BR (absolute) p"] != 0):
                if (index == df.shape[0]-1):
                    # last row is the total decay mode probability
                    continue
                if (ONLY_FERMI == True and float(row['Levels Energy']) != new_Ex):
                    # only write the decay to the ground state
                    continue
                f.write(f"                               MshellEC            {row['Levels Energy']}  -       {row['BR (absolute) ecMp']}     {row['Qβ']}\n")

    digit = 2
    ########################### 32Cl* --------> 31S + p ###########################
    shutil.copy("../32Ar/z17.a32_base", f"{name}/z17.a32")

    counter_peak = 0
    list_peak = []
    with open(f"{name}/z17.a32", "a") as f:
        for index, row in df.iterrows(): 
            if (row["BR (absolute) p"] != 0):
                if (index == df.shape[0]-1):
                    # last row is the total decay mode probability
                    continue
                if (ONLY_FERMI == True and float(row['Levels Energy']) != new_Ex):
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
