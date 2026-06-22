import pandas as pd
import shutil
import numpy as np
import os

########################### 32Ar --------> 32Cl + e + nu ###########################
for Qbeta in range(500, 10000+100, 100):
    # copy and past without base in the filename the file z18.a32_base
    os.mkdir(f"{Qbeta}")
    shutil.copy(f"z18.a32_base", f"{Qbeta}/z18.a32")

    # Build z18.a32 only for beta p discrading gammas
    with open(f"{Qbeta}/z18.a32", "a") as f:
        # write total decay mode probability
        f.write("                               BetaPlus            0  -       1\n")
        f.write(f"                               BetaPlus            {11134.3-Qbeta}  -       100     {Qbeta}\n")
    f.close()