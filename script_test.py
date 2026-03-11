import os
import time

while True:
    for f in os.listdir('/run/media/local1/DATANEX/Samuel-G4/Positions/Detectors/'):
        if f.endswith('_analysed_new.root'):
            if os.path.exists(f'/run/media/local1/DATANEX/Samuel-G4/Positions/Detectors/{f[:-18]}_result.root'):
                continue
            os.system(f'cd Grouper/; ./Analyser_Sim /run/media/local1/DATANEX/Samuel-G4/Positions/Detectors/{f}; cd -')
    print(" Waiting for the next cycle...")
    time.sleep(900)



# #!/bin/bash
# cd Grouper
# # for all the file in dir /run/media/local1/DATANEX/Samuel-G4/PositionCatcherAngle/ run ReaderNew $f
# for f in /run/media/local1/DATANEX/Samuel-G4/RelativePosition/*_analysed_new.root;
# do
#     ./Analyser_Sim "$f"
# done
