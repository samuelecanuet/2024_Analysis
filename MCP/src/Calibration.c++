#include "Calibration.hh"




int main(int argc, char *argv[])
{
    string Run = "";
    if (argc > 1)
    {
        Run = argv[1];
    }
    
    FLAG2025 = true;
    InitDetectors("Config_Files/sample.pid");

    year = to_string(YEAR);
    InitYearConfiguration(year, Run);

    WRITTING = true;

    ROOT_file = MyTFile(Calibration_Filename.c_str(), "READ");
    H_precorrrected = (TH2D*)ROOT_file->Get("h_Image");
    H_precorrrected->SetMaximum(MAXI);
    H_precorrrected->GetXaxis()->SetRangeUser(-n, n);
    H_precorrrected->GetYaxis()->SetRangeUser(-n, n);

    if (WRITTING)
    {
        FINAL_file = MyTFile((DIR_ROOT_DATA_MCP_CALIBRATED +  Measurement_Filename.substr(0, Measurement_Filename.find("_grouped.root")) + "_calibrated.root").c_str(), "RECREATE");
        H_precorrrected->Write();
    }


    //1D
    // MyMinimizer();

    //2D
    // MyMinimizer2D();

    

    ///////////////////////// REDO 2D GRID FIT /////////////////////////
    // ### Loading guessed centers from guess_centers.txt
    LoadCenters(year);
    // ### Fitting the grid to get centers/corner for calibration purpose
    // FullFunctionToMinimize2D_CORNER();

    ///////////////////////// REDO POLY FIT /////////////////////////
    // Reconstruction
    // ### Fitting polynomial fucntion X and Y from centers point written in out_centers.txt
    // FittingReconstruction_CORNER(); 
    // ### Redo the grid fit in the extended region of interest to get the MCP resolution (not used)
    // FinalFunctionToMinimize2D();
    
    // TFile *f = MyTFile((DIR_ROOT_DATA_MCP_CALIBRATED + "run_006_MCP_Scan_calibrated.root").c_str(), "READ");
    // TCanvas *c = (TCanvas*)f->Get("Measurement_2D_View");
    // if (c == nullptr)
    // {
    //     cout << "Error: Canvas is nullptr" << endl;
    //     return -1;
    // }
    // H_reconstruction_interpolated = (TH2D*)c->GetPrimitive("H_measurement_reconstructed");
    // if (H_reconstruction_interpolated == nullptr)
    // {
    //     cout << "Error: H_reconstruction_interpolated is nullptr" << endl;
    //     return -1;
    // }
    
    // TFile *fout = MyTFile((DIR_ROOT_DATA_MCP_CALIBRATED + "run_001_MCP_Scan_calibrated_secondstep.root").c_str(), "RECREATE");    
    // // FittingResolution_CORNER(H_reconstruction_interpolated);

    // TH2D* H = SecondStepFit_CORNER(H_reconstruction_interpolated, fout);
    // fout ->cd();
    // H->Write();
    // TFile *fint = MyTFile("MSecondStepCalibrated_006.root", "READ");
    // TH2D* H = (TH2D*)fint->Get("H_reconstruction_interpolated_second_measurement");
    // TFile *fout = MyTFile("MSecondStepCalibrated_006_resolution.root", "RECREATE");
    // FittingResolution_CORNER(H);
    // fout->Close();
    //Measurement
    // TFile *file_measurement = MyTFile((DIR_ROOT_DATA_MCP_GROUPED+Measurement_Filename).c_str(), "READ");
    // H_measurement = (TH2D*)file_measurement->Get("anode_image_func");
    // H_measurement = (TH2D*)file_measurement->Get("h_Image");

    // FINAL_file->cd();
    // // ### Regenerate the experimental calibrated data with f_X and f_Y computed earlier in FittingReconstruction()
    // H_measurement->Write();
    // Measurement2D();
    string filenamein = "run_030_MCP_32Ar_Beam_1T_calibrated_secondstep.root";
    string filenamefout = filenamein.substr(0, filenamein.find(".root")) + "_resolution2.root";
    TFile *fin = MyTFile(filenamein.c_str(), "READ");
    TH2D *H_measurement_reconstructed = (TH2D*)fin->Get("H_reconstruction_interpolated_second_measurement");
    TFile *fout = MyTFile(filenamefout.c_str(), "RECREATE");
    FittingResolution_CORNER(H_measurement_reconstructed);
    fout->Close();

    // /////////////////////////////////////////////////////////////////

    // fSaved = new TFile("Calibration_Saved.root", "READ");
    // H_measurement_reconstructed = (TH2D*)fSaved->Get("H_measurement_reconstructed");
    // double chi2 = MeasurementFunctionToMinimize2D_CORNER();
    

    // ofstream file;
    // file.open("tmp/"+to_string(thread) + ".txt");
    // file << chi2;
    // file.close();

    if (WRITTING)
    {
        PrintResults();
        FINAL_file->Close();
    }

    return 0;
}
