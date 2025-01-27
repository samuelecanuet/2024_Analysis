#include "Calibration.hh"

int main(int argc, char *argv[])
{
    int thread = atoi(argv[9]);

    WRITTING = true;

    

    ROOT_file = new TFile("MCP_008_4T_0001.root", "READ");
    H_precorrrected = (TH2D*)ROOT_file->Get("h_Image");
    H_precorrrected->SetMaximum(MAXI);
    H_precorrrected->GetXaxis()->SetRangeUser(-n, n);
    H_precorrrected->GetYaxis()->SetRangeUser(-n, n);

    if (WRITTING)
    {
        FINAL_file = new TFile("Calibration.root", "RECREATE");
        H_precorrrected->Write();
    }


    //1D
    // MyMinimizer();

    //2D
    // MyMinimizer2D();

    ///////////////////////// REDO 2D GRID FIT /////////////////////////
    // ### Loading guessed centers from guess_centers.txt
    LoadCenters();
    // ### Fitting the grid to get centers for calibration purpose
    FullFunctionToMinimize2D();

    ///////////////////////// REDO POLY FIT /////////////////////////
    // Reconstruction
    // ### Fitting polynomial fucntion X and Y from centers point written in guess_centers.txt
    FittingReconstruction(); 
    // ### Redo the grid fit in the extended region of interest to get the MCP resolution (not used)
    // FinalFunctionToMinimize2D();
    //Measurement
    TFile *file_measurement = new TFile("../../../../../../mnt/hgfs/shared-2/2024_DATA/MCP_DATA/RAW/MCP_010_4T.root", "READ");
    H_measurement = (TH2D*)file_measurement->Get("anode_image_func");
    FINAL_file->cd();
    // ### Regenerate the experimental calibrated data with f_X and f_Y computed earlier in FittingReconstruction()
    H_measurement->Write();
    Measurement2D();
    /////////////////////////////////////////////////////////////////

    double par[8] = {atof(argv[1]), atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]), atof(argv[8])};

    fSaved = new TFile("Calibration_Saved.root", "READ");
    H_measurement_reconstructed = (TH2D*)fSaved->Get("H_measurement_reconstructed");
    double chi2 = MeasurementFunctionToMinimize2D(par);
    

    ofstream file;
    file.open("tmp/"+to_string(thread) + ".txt");
    file << chi2;
    file.close();

    if (WRITTING)
    {
        PrintResults();
        FINAL_file->Close();
    }

    return 0;
}
