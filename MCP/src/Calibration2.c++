#include "Calibration2.hh"

int main(int argc, char *argv[])
{
    ROOT_file = new TFile("../../MCP_program/MCP_010_4T.root", "READ");
    H_precorrrected = (TH2D*)ROOT_file->Get("h_Image_corr");


    ROOT_file = new TFile("MCP_008_4T_0001.root", "READ");
    H_precorrrected = (TH2D*)ROOT_file->Get("h_Image");
    H_precorrrected->SetMaximum(MAXI);
    H_precorrrected->GetXaxis()->SetRangeUser(-n, n);
    H_precorrrected->GetYaxis()->SetRangeUser(-n, n);
    // Tree = (TTree *)ROOT_file->Get("Tree_Group");
    // Reader = new TTreeReader(Tree);
    // TTreeReaderArray<double> Channel(*Reader, "Channel");
    // TTreeReaderArray<int> Label(*Reader, "Label");


    FINAL_file = new TFile("Calibration.root", "RECREATE");
    H_precorrrected->Write();
    

    //1D
   // MyMinimizer();

    //2D
    // MyMinimizer2D();

    // 2D Full
    LoadCenters();
    FittingReconstruction();
    FullFunctionToMinimize2D();

    // Reconstruction
    // FinalFunctionToMinimize2D();
    
    //Measurement
    // TFile *file_measurement = new TFile("../../../../../../mnt/hgfs/shared-2/2024_DATA/MCP_DATA/RAW/MCP_010_4T.root", "READ");
    // H_measurement = (TH2D*)file_measurement->Get("anode_image_func");
    
    FINAL_file->cd();
    // H_measurement->Write();
    // Measurement2D();
    // MeasurementFunctionToMinimize2D();

    // testdiag();

    FINAL_file->Close();

    return 0;
}
