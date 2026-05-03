import os
import time
import subprocess

directories = ["Bfield", "Catcher_thickness", "DL", "a_Calibration", "Sampling", "Silicon_Position", "Sensibility_B"]

# inpath = '/run/media/local1/DATANEX/Samuel-G4/Positions/Detectors/AsyEfficiency/'


import shlex

remote_base = "/data333/lecanuet/Result/Systematics/"
remote_host = "lecanuet@borlin333"
# inpath = "/run/media/local1/Disque_Dur/2025_DATA/SIMULATED_DATA/Systematics/"
inpath = "/run/media/local1/DATANEX/Samuel-G4/Systematics/"

while True:
    for d in directories:
        print(f"----- Checking directory: {d}")
        local_dir = os.path.join(inpath, d)
        remote_dir = f"{remote_base}{d}/"

        os.makedirs(local_dir, exist_ok=True)

        # List remote files matching "*new*"
        cmd_list = [
            "ssh",
            remote_host,
            f'find {shlex.quote(remote_dir)} -maxdepth 1 -type f -name "*new*"'
        ]

        result = subprocess.run(cmd_list, stdout=subprocess.PIPE)
        

        if result.returncode != 0:
            print(f"Error while listing remote files for {d}:")
            print(result.stderr)
            continue

        remote_files = result.stdout.strip().splitlines()

        for remote_file in remote_files:
            # print(f"Remote file found: {remote_file}")
            fname = str(os.path.basename(remote_file))[2:-1]
            # print(f"Processing remote file: {fname}")

            # Example:
            # xxx_analysed_new.root  ->  xxx_result.root
            if "_analysed_new.root" in fname:
                result_name = fname[:-18] + "_result.root"
            else:
                # Adapt here if you also want other "*new*" patterns
                continue

            print("### Result file to check/sync: ", local_dir+"/" + result_name)
            local_result = os.path.join(local_dir, result_name)

            # Skip rsync if corresponding result file already exists locally
            if os.path.exists(local_result):
                # print(f"Result file {result_name} already exists locally. Skipping sync.")
                continue
            else:
                print(f"Result file {result_name} not found locally. Syncing remote file...")

            # Sync the remote file
            cmd_rsync = [
                "rsync", "-av",
                f"{remote_host}:{remote_dir+fname}",
                local_dir + "/"
            ]

            # print(f"Running command: {' '.join(cmd_rsync)}")

            rs = subprocess.run(cmd_rsync)
            if rs.returncode != 0:
                print(f"Error while syncing {remote_file}")
        print("\n")
        # Run the analyser
        for f in os.listdir(local_dir):
            if f.endswith("_analysed_new.root"):
                if os.path.exists(f"{local_dir}/{f[:-18]}_result.root"):
                    continue
                os.system(f'cd Grouper/; ./Analyser_Sim "{local_dir}/{f}"; cd -')

    print("Waiting for the next cycle...")
    time.sleep(900)



# #!/bin/bash
# cd Grouper
# # for all the file in dir /run/media/local1/DATANEX/Samuel-G4/PositionCatcherAngle/ run ReaderNew $f
# for f in /run/media/local1/DATANEX/Samuel-G4/RelativePosition/*_analysed_new.root;
# do
#     ./Analyser_Sim "$f"
# done
