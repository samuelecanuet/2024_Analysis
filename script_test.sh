#!/bin/bash
cd Grouper
# for all the file in dir /run/media/local1/DATANEX/Samuel-G4/PositionCatcherAngle/ run ReaderNew $f
for f in /run/media/local1/DATANEX/Samuel-G4/RelativePosition/*_analysed_new.root;
do
    ./Analyser_Sim "$f"
done
