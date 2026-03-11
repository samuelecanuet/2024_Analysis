#include "../Grouper/include/Detectors.hh"
#include "../Grouper/include/Utilities.hh"

#include <iostream>
using namespace std;
#include "TH1D.h"
#include "TGraph2D.h"


int ApplySecondStepCalibration(string GroupedFilename="")
{   


    GroupedFilename = "run_034_MCP_32Ar_BeamScan_2T_grouped.root";
    FLAG2025 = true;  
    InitDetectors("../Grouper/Config_Files/sample.pid");
    // First step interpolation
    TFile *fcalib = MyTFile((DIR_ROOT_DATA_MCP_CALIBRATED + "run_001_StableBeamScan_calibrated.root").c_str(), "READ");
    TCanvas *c_calibx = (TCanvas*)fcalib->Get("c_fit_X");
    TGraph2D *G2D_X;
    for (int i = 0 ; i < c_calibx->GetListOfPrimitives()->GetEntries() ; i++)
    {
        TObject *obj = c_calibx->GetListOfPrimitives()->At(i);
        if (obj->InheritsFrom("TGraph2D"))
        {
            G2D_X = (TGraph2D*)obj;
            break;
        }
    }

    TCanvas *c_caliby = (TCanvas*)fcalib->Get("c_fit_Y");
    TGraph2D *G2D_Y;
    for (int i = 0 ; i < c_caliby->GetListOfPrimitives()->GetEntries() ; i++)
    {
        TObject *obj = c_caliby->GetListOfPrimitives()->At(i);
        if (obj->InheritsFrom("TGraph2D"))
        {
            G2D_Y = (TGraph2D*)obj;
            break;
        }
    }

    //SEcond step interpolation
    TGraph2D *G_SecondInterpolation_X;
    TGraph2D *G_SecondInterpolation_Y;
    TFile *fcalibs = MyTFile((DIR_ROOT_DATA_MCP_CALIBRATED + "run_001_StableBeamScan_secondstep_calibrated(resolutionbeforesecond).root").c_str(), "READ");
    for (int i = 0 ; i < fcalibs->GetListOfKeys()->GetEntries() ; i++)
    {
        TKey *key = (TKey*)fcalibs->GetListOfKeys()->At(i);
        TObject *obj = key->ReadObj();
        if (string(obj->GetName()) == "Graph2D")
        {
            G_SecondInterpolation_X = (TGraph2D*)obj;
            key = (TKey*)fcalibs->GetListOfKeys()->At(i+1);
            obj = key->ReadObj();
            G_SecondInterpolation_Y = (TGraph2D*)obj;
            break;
        }
    }

    // check that the graphs were found
    if (G2D_X == nullptr || G2D_Y == nullptr || G_SecondInterpolation_X == nullptr || G_SecondInterpolation_Y == nullptr)
    {
        cout << "Error: One or more TGraph2D objects not found." << endl;
        return -1;
    }

    TH2D* H_reconstruction_interpolated_second_measurement = new TH2D("H_reconstruction_interpolated_second_measurement", "H_reconstruction_interpolated_second_measurement", 800, -8, 8, 800, -8, 8);

    string filenameout = GroupedFilename.substr(0, GroupedFilename.find("_grouped.root")) + "_calibrated_secondstep.root";
    TFile *fout2 = new TFile(filenameout.c_str(), "RECREATE");
    TTree * tree_final = new TTree("treeMCP", "treeMCP");
    double X_final, Y_final, TIME_final;
    tree_final->Branch("X", &X_final, "X/D");
    tree_final->Branch("Y", &Y_final, "Y/D");
    tree_final->Branch("Time", &TIME_final, "Time/D");

    TFile * fin = MyTFile((DIR_ROOT_DATA_MCP_GROUPED + GroupedFilename).c_str(), "READ");
    TTree *tree_meas = (TTree *)fin->Get("treeMCP");
    TTreeReader *Reader_meas = new TTreeReader(tree_meas);
    TTreeReaderValue<double> *X0_meas = new TTreeReaderValue<double>(*Reader_meas, "X");
    TTreeReaderValue<double> *Y0_meas = new TTreeReaderValue<double>(*Reader_meas, "Y");
    TTreeReaderValue<double> *TIME0_meas = new TTreeReaderValue<double>(*Reader_meas, "Time");
    cout << "Starting to process events..." << endl;
    while (Reader_meas->Next())
    {
        double x = **X0_meas;
        double y = **Y0_meas;
        double x_int = G2D_X->Interpolate(x, y);
        double y_int = G2D_Y->Interpolate(y, x);


        double x_int_second = G_SecondInterpolation_X->Interpolate(x_int, y_int);
        double y_int_second = G_SecondInterpolation_Y->Interpolate(x_int, y_int);

        if (abs(x_int) >= 6.68 || abs(y_int) >= 6.68)
            continue;
        
        H_reconstruction_interpolated_second_measurement->Fill(x_int_second, y_int_second);

        X_final = x_int_second;
        Y_final = y_int_second;
        TIME_final = **TIME0_meas;
        tree_final->Fill();
    }

    cout << "Finished processing events." << endl;

    fin->Close();
    fout2->cd();
    
    
    H_reconstruction_interpolated_second_measurement->Write();
    tree_final->Write();
    fout2->Close();

    return 0;
}