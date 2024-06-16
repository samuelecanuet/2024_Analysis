for run in 57, 58, 59, 60, 61, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 77, 79, 112, 113, 114, 115, 116, 117, 118

do
    cd Analysis/ReadFaster
    FWisardReadRun "../../../../../../../mnt/hgfs/shared-2/2024_Wisard/run_%3R_multifast_32Ar.fast/run_%3R_*_*.fast" -r"$run"
    cd -
    cd Grouper
    Grouper $run
    Matcher $run
    cd -
done

