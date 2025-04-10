#include "MCMC_SiPM.hh"

#include <iostream>
#include <sys/mman.h>
#include <sys/shm.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>


int main(int argc, char **argv)
{
    InitDetectors("../../Grouper/Config_Files/sample.pid");

    double Threshold = stod(argv[1]);
    double Offset_calibration = stod(argv[2]);
    double Coefficients1_calibration = stod(argv[3]);
    double Coefficients1_resolution = stod(argv[4]);
    double Coefficients2_resolution = stod(argv[5]);

    string THREAD = argv[6];

    vector<double> Par = {Offset_calibration, Coefficients1_calibration, 0, Coefficients1_resolution, Coefficients2_resolution, Threshold, 0};

    ///////////////////////////////  INITIALISATION ///////////////////////////////
    InitWindows(0, "../");
    InitSiliconCalibration("../");
    InitHistograms(0);

    ////////////////////////  EXPERIMETAL TREE FOR EACH PEAK //////////////////////
    InitTree("READ");    

    ////////////////////////  SIMULATION TREE FOR EACH PEAK //////////////////////
    InitSimulatedTree();

    current_detector = 7;
    const double *bestPar = Par.data();
    VERBOSE = 1;

    Start("Detector: SiPM " + to_string(current_detector));
    double chi2 = Chi2TreeHist_conv(bestPar);

    const char* shared_mem_name = THREAD.c_str();
    const int SIZE = sizeof(double);

    int shm_fd = shm_open(shared_mem_name, O_CREAT | O_RDWR, 0666);
    ftruncate(shm_fd, SIZE);
    void* ptr = mmap(0, SIZE, PROT_WRITE, MAP_SHARED, shm_fd, 0);

    if (ptr == MAP_FAILED) {
        std::cerr << "Memory mapping failed!" << std::endl;
        return 1;
    }

    double value = chi2; 
    memcpy(ptr, &value, sizeof(double));
    munmap(ptr, SIZE);
    close(shm_fd);


    cout << "Chi2: " << chi2 << endl;
    return 0;

}