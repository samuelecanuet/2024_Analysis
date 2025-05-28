#include "Defect.hh"
#include "Math/MinimizerOptions.h"
int main()
{

    string DET = "D5.1";
    ///////////////////////////////////  FILES //////////////////////////////////
    GROUPED_File["32Ar"] = new TFile((DIR_ROOT_DATA_GROUPED+"merge3.root").c_str(), "READ");
    GROUPED_File["4alpha"] = new TFile((DIR_ROOT_DATA_GROUPED + "run_127_data_4alpha_grouped.root").c_str(), "READ");

    SIMULATED_File["18N"] = new TFile((DIR_ROOT_DATA_SIMULATED + "18N__CS0_CSP0_CV1_CVP1.root").c_str(), "READ");
    SIMULATED_File["148Gd"] = new TFile((DIR_ROOT_DATA_SIMULATED + "148Gd_700nm_width_analysed.root").c_str(), "READ");
    SIMULATED_File["239Pu"] = new TFile((DIR_ROOT_DATA_SIMULATED + "239Pu_700nm_width_analysed.root").c_str(), "READ");
    SIMULATED_File["241Am"] = new TFile((DIR_DATA_HDD + "../SIMULATED_DATA/03-09/" + "241Am_inter0025_a1_b0_analysed.root").c_str(), "READ");
    SIMULATED_File["244Cm"] = new TFile((DIR_DATA_HDD + "../SIMULATED_DATA/03-09/" + "244Cm_inter0025_a1_b0_analysed.root").c_str(), "READ");

    ///////////////////////////////////  OUTPUT ///////////////////////////////////
    DEFECT_File = new TFile((DIR_ROOT_DATA_MATCHED + "defect2.root").c_str(), "RECREATE");
    DEFECT_File->cd();
    // WriteTime(GROUPED_File[], MATCHED_File);

    ///////////////////////////////////////////////////////////////////////////////

    TH1D *halphas = (TH1D *)GROUPED_File["4alpha"]->Get(("Strip/Strip_Channel/"+DET+"/Channel_Cleaned_"+DET).c_str());

    ///////////////////////////////////////////////////////////////////////////////

    ////////////////////////////// FULL PEAKS //////////////////////////////
    halphas->GetXaxis()->SetRangeUser(WindowExp["148Gd"].first, WindowExp["244Cm"].second);
    TF1 *f1 = new TF1("f1", fullgaussBortelsCollaers, 0, 100000, 22);
    f1->SetNpx(1000000);

    f1->SetParameter(0, 2e5);                       //Gd
    f1->SetParLimits(0, 1.5e5, 2.5e5);
    f1->SetParameter(6, 2e5*0.11/0.70);             //Pu1
    // f1->SetParLimits(6, 1.7e4, 1.9e4);
    f1->SetParameter(8, 1.8e4*0.17/0.11);             //Pu2
    // f1->SetParLimits(8, 1.7e4*0.17/0.11, 1.9e4*0.17/0.11);
    f1->SetParameter(10, 2e5);                       //Pu3
    // f1->SetParLimits(10, 100, 400000);
    f1->SetParameter(12, 1.5e5*0.0166/0.84);          //Am1
    f1->SetParLimits(12, 2e3, 8e3);
    f1->SetParameter(14, 1.5e5*0.13/0.84);            //Am2
    f1->SetParLimits(14, 2e4, 8e4);
    f1->SetParameter(16, 1.5e5);                      //Am3
    // f1->SetParLimits(16, 100, 400000);
    f1->SetParameter(18, 2e5*0.23/0.77);            //Cm1
    // f1->SetParLimits(18, 100, 400000);
    f1->SetParameter(20, 2e5);                      //Cm2
    // f1->SetParLimits(20, 100, 400000);
    f1->FixParameter(1, Peaks_Exp["148Gd"][0]);
    f1->SetParLimits(1, Peaks_Exp["148Gd"][0]-200, Peaks_Exp["148Gd"][0]+200);
    f1->FixParameter(7, Peaks_Exp["239Pu"][0]);
    f1->SetParLimits(7, Peaks_Exp["239Pu"][0]-200, Peaks_Exp["239Pu"][0]+200);
    f1->FixParameter(9, Peaks_Exp["239Pu"][1]);
    f1->SetParLimits(9, Peaks_Exp["239Pu"][1]-200, Peaks_Exp["239Pu"][1]+200);
    f1->FixParameter(11, Peaks_Exp["239Pu"][2]);
    f1->SetParLimits(11, Peaks_Exp["239Pu"][2]-200, Peaks_Exp["239Pu"][2]+200);
    f1->FixParameter(13, Peaks_Exp["241Am"][0]);
    f1->SetParLimits(13, Peaks_Exp["241Am"][0]-200, Peaks_Exp["241Am"][0]+200);
    f1->FixParameter(15, Peaks_Exp["241Am"][1]);
    f1->SetParLimits(15, Peaks_Exp["241Am"][1]-200, Peaks_Exp["241Am"][1]+200);
    f1->FixParameter(17, Peaks_Exp["241Am"][2]);
    f1->SetParLimits(17, Peaks_Exp["241Am"][2]-200, Peaks_Exp["241Am"][2]+200);
    f1->FixParameter(19, Peaks_Exp["244Cm"][0]);
    f1->SetParLimits(19, Peaks_Exp["244Cm"][0]-200, Peaks_Exp["244Cm"][0]+200);
    f1->FixParameter(21, Peaks_Exp["244Cm"][1]);
    f1->SetParLimits(21, Peaks_Exp["244Cm"][1]-200, Peaks_Exp["244Cm"][1]+200);

    f1->SetParameter(2, 83);
    f1->SetParLimits(2, 40, 90);
    f1->SetParameter(3, 0.02);
    f1->SetParLimits(3, 0.0, 0.05);
    f1->SetParameter(4, 70);
    f1->SetParLimits(4, 10, 300);
    f1->SetParameter(5, 2000);
    f1->SetParLimits(5, 500, 5000);

    cout << "Fitting" << endl;
    halphas->Fit("f1", "RL", "", WindowExp["148Gd"].first, WindowExp["244Cm"].second);
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    halphas->Draw("HIST");
    f1->Draw("SAME");
    c1->Write();
    delete c1;


    // BRANCHING RATIO CHECK //
    // 239Pu
    TF1 *f_Pu239_1 = new TF1("f_Pu239_1", gaussBortelsCollaers, 0, 100000, 6);
    f_Pu239_1->SetParameter(0, f1->GetParameter(6));
    f_Pu239_1->SetParameter(1, f1->GetParameter(7));
    f_Pu239_1->SetParameter(2, f1->GetParameter(2));
    f_Pu239_1->SetParameter(3, f1->GetParameter(3));
    f_Pu239_1->SetParameter(4, f1->GetParameter(4));
    f_Pu239_1->SetParameter(5, f1->GetParameter(5));
    double Pu239_1 = f_Pu239_1->Integral(0, 100000);

    TF1 *f_Pu239_2 = new TF1("f_Pu239_2", gaussBortelsCollaers, 0, 100000, 6);
    f_Pu239_2->SetParameter(0, f1->GetParameter(8));
    f_Pu239_2->SetParameter(1, f1->GetParameter(9));
    f_Pu239_2->SetParameter(2, f1->GetParameter(2));
    f_Pu239_2->SetParameter(3, f1->GetParameter(3));
    f_Pu239_2->SetParameter(4, f1->GetParameter(4));
    f_Pu239_2->SetParameter(5, f1->GetParameter(5));
    double Pu239_2 = f_Pu239_2->Integral(0, 100000);

    TF1 *fPu239_3 = new TF1("fPu239_3", gaussBortelsCollaers, 0, 100000, 6);
    fPu239_3->SetParameter(0, f1->GetParameter(10));
    fPu239_3->SetParameter(1, f1->GetParameter(11));
    fPu239_3->SetParameter(2, f1->GetParameter(2));
    fPu239_3->SetParameter(3, f1->GetParameter(3));
    fPu239_3->SetParameter(4, f1->GetParameter(4));
    fPu239_3->SetParameter(5, f1->GetParameter(5));
    double Pu239_3 = fPu239_3->Integral(0, 100000);

    double BrPu239_1 = Pu239_1/Pu239_2;
    double BrPu239_1_err = sqrt( pow(Pu239_1/pow(Pu239_2, 2) * sqrt(Pu239_2), 2 ) + pow(1/Pu239_2 * sqrt(Pu239_1), 2 ));
    double BrPu239_2 = Pu239_3/Pu239_2;
    double BrPu239_2_err = sqrt( pow(Pu239_3/pow(Pu239_2, 2) * sqrt(Pu239_2), 2 ) + pow(1/Pu239_2 * sqrt(Pu239_3), 2 ));

    // 241Am
    TF1 *f_Am241_1 = new TF1("f_Am241_1", gaussBortelsCollaers, 0, 100000, 6);
    f_Am241_1->SetParameter(0, f1->GetParameter(12));
    f_Am241_1->SetParameter(1, f1->GetParameter(13));
    f_Am241_1->SetParameter(2, f1->GetParameter(2));
    f_Am241_1->SetParameter(3, f1->GetParameter(3));
    f_Am241_1->SetParameter(4, f1->GetParameter(4));
    f_Am241_1->SetParameter(5, f1->GetParameter(5));
    double Am241_1 = f_Am241_1->Integral(0, 100000);

    TF1 *f_Am241_2 = new TF1("f_Am241_2", gaussBortelsCollaers, 0, 100000, 6);
    f_Am241_2->SetParameter(0, f1->GetParameter(14));
    f_Am241_2->SetParameter(1, f1->GetParameter(15));
    f_Am241_2->SetParameter(2, f1->GetParameter(2));
    f_Am241_2->SetParameter(3, f1->GetParameter(3));
    f_Am241_2->SetParameter(4, f1->GetParameter(4));
    f_Am241_2->SetParameter(5, f1->GetParameter(5));
    double Am241_2 = f_Am241_2->Integral(0, 100000);

    TF1 *f_Am241_3 = new TF1("f_Am241_3", gaussBortelsCollaers, 0, 100000, 6);
    f_Am241_3->SetParameter(0, f1->GetParameter(16));
    f_Am241_3->SetParameter(1, f1->GetParameter(17));
    f_Am241_3->SetParameter(2, f1->GetParameter(2));
    f_Am241_3->SetParameter(3, f1->GetParameter(3));
    f_Am241_3->SetParameter(4, f1->GetParameter(4));
    f_Am241_3->SetParameter(5, f1->GetParameter(5));
    double Am241_3 = f_Am241_3->Integral(0, 100000);

    double BrAm241_1 = Am241_1/Am241_3;
    double BrAm241_1_err = sqrt( pow(Am241_1/pow(Am241_3, 2) * sqrt(Am241_3), 2 ) + pow(1/Am241_3 * sqrt(Am241_1), 2 ));
    double BrAm241_2 = Am241_2/Am241_3;
    double BrAm241_2_err = sqrt( pow(Am241_2/pow(Am241_3, 2) * sqrt(Am241_3), 2 ) + pow(1/Am241_3 * sqrt(Am241_2), 2 ));

    // 244Cm
    TF1 *f_Cm244_1 = new TF1("f_Cm244_1", gaussBortelsCollaers, 0, 100000, 6);
    f_Cm244_1->SetParameter(0, f1->GetParameter(18));
    f_Cm244_1->SetParameter(1, f1->GetParameter(19));
    f_Cm244_1->SetParameter(2, f1->GetParameter(2));
    f_Cm244_1->SetParameter(3, f1->GetParameter(3));
    f_Cm244_1->SetParameter(4, f1->GetParameter(4));
    f_Cm244_1->SetParameter(5, f1->GetParameter(5));
    double Cm244_1 = f_Cm244_1->Integral(0, 100000);

    double BrCm244_1 = Cm244_1/Am241_3;
    double BrCm244_1_err = sqrt( pow(Cm244_1/pow(Am241_3, 2) * sqrt(Am241_3), 2 ) + pow(1/Am241_3 * sqrt(Cm244_1), 2 ));


    cout << "----- 239Pu ----- " << endl;
    cout << "Pu1 : " << BrPu239_1 << " +/- " << BrPu239_1_err << " : " << 11.94/70.77 << " +/- " << sqrt(pow(11.94/pow(70.77,2) * 0.14, 2) + pow(1/70.77 * 0.07, 2)) << endl;
    cout << "Pu2 : " << BrPu239_2 << " +/- " << BrPu239_2_err << " : " << 17.11/70.77 << " +/- " << sqrt(pow(17.11/pow(70.77,2) * 0.14, 2) + pow(1/70.77 * 0.14, 2)) << endl;

    cout << "----- 241Am ----- " << endl;
    cout << "Am1 : " << BrAm241_1 << " +/- " << BrAm241_1_err << " : " << 1.66/84.8 << " +/- " << sqrt(pow(1.66/pow(84.8,2) * 0.5, 2) + pow(1/84.8 * 0.2, 2)) << endl;
    cout << "Am2 : " << BrAm241_2 << " +/- " << BrAm241_2_err << " : " << 13.1/84.8 << " +/- " << sqrt(pow(13.1/pow(84.8,2) * 0.5, 2) + pow(1/84.8 * 0.3, 2)) << endl;

    cout << "----- 244Cm ----- " << endl;
    cout << "Cm1 : " << BrCm244_1 << " +/- " << BrCm244_1_err << " : " << 23.1/76.9 << " +/- " << sqrt(pow(23.1/pow(76.9,2) * 0.1, 2) + pow(1/76.9 * 0.1, 2)) << endl;


    TGraphErrors *G_AlphaBR_Pu = new TGraphErrors();
    G_AlphaBR_Pu->SetPoint(0, 1, BrPu239_1);
    G_AlphaBR_Pu->SetPointError(0, 0, BrPu239_1_err);
    G_AlphaBR_Pu->SetPoint(1, 2, BrPu239_2);
    G_AlphaBR_Pu->SetPointError(1, 0, BrPu239_2_err);
    TGraphErrors *G_AlphaBR_Pu_Lit = new TGraphErrors();
    G_AlphaBR_Pu_Lit->SetPoint(0, 1, 11.94/70.77);
    G_AlphaBR_Pu_Lit->SetPointError(0, 0, sqrt(pow(11.94/pow(70.77,2) * 0.14, 2) + pow(1/70.77 * 0.07, 2)));
    G_AlphaBR_Pu_Lit->SetPoint(1, 2, 17.11/70.77);
    G_AlphaBR_Pu_Lit->SetPointError(1, 0, sqrt(pow(17.11/pow(70.77,2) * 0.14, 2) + pow(1/70.77 * 0.14, 2)));

    TGraphErrors *G_AlphaBr_Am = new TGraphErrors();
    G_AlphaBr_Am->SetPoint(0, 3, BrAm241_1);
    G_AlphaBr_Am->SetPointError(0, 0, BrAm241_1_err);
    G_AlphaBr_Am->SetPoint(1, 4, BrAm241_2);
    G_AlphaBr_Am->SetPointError(1, 0, BrAm241_2_err);
    TGraphErrors *G_AlphaBr_Am_Lit = new TGraphErrors();
    G_AlphaBr_Am_Lit->SetPoint(0, 3, 1.66/84.8);
    G_AlphaBr_Am_Lit->SetPointError(0, 0, sqrt(pow(1.66/pow(84.8,2) * 0.5, 2) + pow(1/84.8 * 0.2, 2)));
    G_AlphaBr_Am_Lit->SetPoint(1, 4, 13.1/84.8);
    G_AlphaBr_Am_Lit->SetPointError(1, 0, sqrt(pow(13.1/pow(84.8,2) * 0.5, 2) + pow(1/84.8 * 0.3, 2)));

    TGraphErrors *G_AlphaBr_Cm = new TGraphErrors();
    G_AlphaBr_Cm->SetPoint(0, 5, BrCm244_1);
    G_AlphaBr_Cm->SetPointError(0, 0, BrCm244_1_err);
    TGraphErrors *G_AlphaBr_Cm_Lit = new TGraphErrors();
    G_AlphaBr_Cm_Lit->SetPoint(0, 5, 23.1/76.9);
    G_AlphaBr_Cm_Lit->SetPointError(0, 0, sqrt(pow(23.1/pow(76.9,2) * 0.1, 2) + pow(1/76.9 * 0.1, 2)));

    TCanvas *cBR = new TCanvas("cBR", "cBR", 800, 600);
    cBR->cd();
    // G_AlphaBR_Pu->Draw("AP");
    G_AlphaBR_Pu->SetMarkerStyle(20);
    G_AlphaBR_Pu->SetMarkerSize(1);
    G_AlphaBR_Pu->SetMarkerColor(kRed);

    // G_AlphaBR_Pu_Lit->Draw("P SAME");
    G_AlphaBR_Pu_Lit->SetMarkerStyle(54);
    G_AlphaBR_Pu_Lit->SetMarkerSize(1);
    G_AlphaBR_Pu_Lit->SetMarkerColor(kRed);

    // G_AlphaBr_Am->Draw("P SAME");
    G_AlphaBr_Am->SetMarkerStyle(20);
    G_AlphaBr_Am->SetMarkerSize(1);
    G_AlphaBr_Am->SetMarkerColor(kBlue);

    // G_AlphaBr_Am_Lit->Draw("P SAME");
    G_AlphaBr_Am_Lit->SetMarkerStyle(54);
    G_AlphaBr_Am_Lit->SetMarkerSize(1);
    G_AlphaBr_Am_Lit->SetMarkerColor(kBlue);

    // G_AlphaBr_Cm->Draw("P SAME");
    G_AlphaBr_Cm->SetMarkerStyle(20);
    G_AlphaBr_Cm->SetMarkerSize(1);
    G_AlphaBr_Cm->SetMarkerColor(kGreen);

    // G_AlphaBr_Cm_Lit->Draw("P SAME");
    G_AlphaBr_Cm_Lit->SetMarkerStyle(54);
    G_AlphaBr_Cm_Lit->SetMarkerSize(1);
    G_AlphaBr_Cm_Lit->SetMarkerColor(kGreen);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(G_AlphaBR_Pu);
    mg->Add(G_AlphaBR_Pu_Lit);
    mg->Add(G_AlphaBr_Am);
    mg->Add(G_AlphaBr_Am_Lit);
    mg->Add(G_AlphaBr_Cm);
    mg->Add(G_AlphaBr_Cm_Lit);
    mg->SetTitle("Alpha Branching Ratios");
    mg->GetXaxis()->SetTitle("Peak number");
    mg->GetYaxis()->SetTitle("BR_i/BR_main");
    mg->Draw("AP");


    //Add Legend
    TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
    legend->AddEntry(G_AlphaBR_Pu, "239Pu Fit", "p");
    legend->AddEntry(G_AlphaBr_Am, "241Am Fit", "p");
    legend->AddEntry(G_AlphaBr_Cm, "244Cm Fit", "p");
    legend->AddEntry(G_AlphaBR_Pu_Lit, "239Pu Literature", "p");
    legend->AddEntry(G_AlphaBr_Am_Lit, "241Am Literature", "p");
    legend->AddEntry(G_AlphaBr_Cm_Lit, "244Cm Literature", "p");
    legend->Draw("SAME");
    cBR->Write();




    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////   CALIBRATION   /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////   
    //// SOURCES ////
    /////////////////
    TH1D* H_Sim_Gd = (TH1D*)SIMULATED_File["148Gd"]->Get(("Silicon_Detector_Energy_Deposit_"+DET+"_All").c_str());
    H_Sim_Gd->Reset();
    for (int i = 1; i <= 8; i++)
    {
        H_Sim_Gd->Add((TH1D*)SIMULATED_File["148Gd"]->Get(("Silicon_Detector_Energy_Deposit_D"+to_string(i)+".1_All").c_str()));
    }
    TH1D* H_Sim_Pu = (TH1D*)SIMULATED_File["239Pu"]->Get(("Silicon_Detector_Energy_Deposit_"+DET+"_All").c_str());
    H_Sim_Pu->Reset();
    for (int i = 1; i <= 8; i++)
    {
        H_Sim_Pu->Add((TH1D*)SIMULATED_File["239Pu"]->Get(("Silicon_Detector_Energy_Deposit_D"+to_string(i)+".1_All").c_str()));
    }

    TH1D* H_Sim_Am = (TH1D*)SIMULATED_File["241Am"]->Get(("Silicon_Detector_Energy_Deposit_"+DET+"_All").c_str());
    H_Sim_Am->Reset();
    for (int i = 1; i <= 8; i++)
    {
        H_Sim_Am->Add((TH1D*)SIMULATED_File["241Am"]->Get(("Silicon_Detector_Energy_Deposit_D"+to_string(i)+".1_All").c_str()));
    }

    TH1D* H_Sim_Cm = (TH1D*)SIMULATED_File["244Cm"]->Get(("Silicon_Detector_Energy_Deposit_"+DET+"_All").c_str());
    H_Sim_Cm->Reset();
    for (int i = 1; i <= 8; i++)
    {
        H_Sim_Cm->Add((TH1D*)SIMULATED_File["244Cm"]->Get(("Silicon_Detector_Energy_Deposit_D"+to_string(i)+".1_All").c_str()));
    }


    H_Sim_Gd->GetXaxis()->SetRangeUser(Peaks_Sim["148Gd"][0]-50, Peaks_Sim["148Gd"][0]+50);
    H_Sim_Pu->GetXaxis()->SetRangeUser(Peaks_Sim["239Pu"][0]-50, Peaks_Sim["239Pu"][2]+50);
    H_Sim_Am->GetXaxis()->SetRangeUser(Peaks_Sim["241Am"][0]-50, Peaks_Sim["241Am"][2]+50);
    H_Sim_Cm->GetXaxis()->SetRangeUser(Peaks_Sim["244Cm"][0]-50, Peaks_Sim["244Cm"][1]+50);

    double peak_uncertainty = 10;
    TF1 *fGd = new TF1("fGd", gauss, 0, 10000, 3);
    fGd->SetNpx(1000000);
    fGd->SetParLimits(0, 1, 1e7);
    fGd->SetParameter(0, 5.6e3);
    fGd->SetParameter(1, Peaks_Sim["148Gd"][0]);
    fGd->SetParLimits(1, Peaks_Sim["148Gd"][0]-peak_uncertainty, Peaks_Sim["148Gd"][0]+peak_uncertainty);
    fGd->SetParLimits(2, 2, 5);
    fGd->SetParameter(2, 3.7);
    H_Sim_Gd->Fit("fGd", "LR", "", Peaks_Sim["148Gd"][0]-peak_uncertainty, Peaks_Sim["148Gd"][0]+peak_uncertainty);

    TF1 *fPu = new TF1("fPu", triplegauss, 0, 10000, 7);
    fPu->SetNpx(1000000);   
    fPu->SetParLimits(0, 1, 5000); 
    fPu->SetParLimits(3, 1, 5000); 
    fPu->SetParLimits(5, 1, 1e6); 
    fPu->SetParameter(1, Peaks_Sim["239Pu"][0]);
    fPu->SetParLimits(1, Peaks_Sim["239Pu"][0]-peak_uncertainty, Peaks_Sim["239Pu"][0]+peak_uncertainty);
    fPu->SetParameter(4, Peaks_Sim["239Pu"][1]);
    fPu->SetParLimits(4, Peaks_Sim["239Pu"][1]-peak_uncertainty, Peaks_Sim["239Pu"][1]+peak_uncertainty);
    fPu->SetParameter(6, Peaks_Sim["239Pu"][2]);
    fPu->SetParLimits(6, Peaks_Sim["239Pu"][2]-peak_uncertainty, Peaks_Sim["239Pu"][2]+peak_uncertainty);
    fPu->SetParLimits(2, 1, 50);
    H_Sim_Pu->Fit("fPu", "RL", "", Peaks_Sim["239Pu"][0]-peak_uncertainty, Peaks_Sim["239Pu"][2]+peak_uncertainty);

    TF1 *fAm = new TF1("fAm", triplegauss, 0, 10000, 7);
    fAm->SetNpx(1000000);
    fAm->SetParLimits(0, 0.1, 1e5);
    fAm->SetParLimits(3, 0.1, 1e5);
    fAm->SetParLimits(5, 0.1, 1e5);
    fAm->SetParameter(1, Peaks_Sim["241Am"][0]);
    fAm->SetParLimits(1, Peaks_Sim["241Am"][0]-peak_uncertainty, Peaks_Sim["241Am"][0]+peak_uncertainty);
    fAm->SetParameter(4, Peaks_Sim["241Am"][1]);
    fAm->SetParLimits(4, Peaks_Sim["241Am"][1]-peak_uncertainty, Peaks_Sim["241Am"][1]+peak_uncertainty);
    fAm->SetParameter(6, Peaks_Sim["241Am"][2]);
    fAm->SetParLimits(6, Peaks_Sim["241Am"][2]-peak_uncertainty, Peaks_Sim["241Am"][2]+peak_uncertainty);
    fAm->SetParLimits(2, 1, 50);
    H_Sim_Am->Fit("fAm", "RL", "", Peaks_Sim["241Am"][0]-peak_uncertainty, Peaks_Sim["241Am"][2]+peak_uncertainty);

    TF1 *fCm = new TF1("fCm", doublegauss, 0, 10000, 5);
    fCm->SetNpx(1000000);
    fCm->SetParLimits(0, 0.1, 1e6); 
    fCm->SetParLimits(3, 0.1, 1e6); 
    fCm->SetParameter(1, Peaks_Sim["244Cm"][0]);
    fCm->SetParLimits(1, Peaks_Sim["244Cm"][0]-peak_uncertainty, Peaks_Sim["244Cm"][0]+peak_uncertainty);
    fCm->SetParameter(4, Peaks_Sim["244Cm"][1]);
    fCm->SetParLimits(4, Peaks_Sim["244Cm"][1]-peak_uncertainty, Peaks_Sim["244Cm"][1]+peak_uncertainty);
    fCm->SetParLimits(2, 1, 50);
    H_Sim_Cm->Fit("fCm", "RL", "", Peaks_Sim["244Cm"][0]-peak_uncertainty, Peaks_Sim["244Cm"][1]+peak_uncertainty);

    TCanvas *cCalibration = new TCanvas("cCalibration", "cCalibration", 800, 600);
    cCalibration->Divide(2, 2);
    cCalibration->cd(1);
    H_Sim_Gd->Draw("HIST");
    fGd->Draw("SAME");
    cCalibration->cd(2);
    H_Sim_Pu->Draw("HIST");
    fPu->Draw("SAME");
    cCalibration->cd(3);
    H_Sim_Am->Draw("HIST");
    fAm->Draw("SAME");
    cCalibration->cd(4);
    H_Sim_Cm->Draw("HIST");
    fCm->Draw("SAME");
    cCalibration->Write();

    TGraphErrors *G_Calibration = new TGraphErrors();
    G_Calibration->SetPoint(0, f1->GetParameter(1), fGd->GetParameter(1));
    G_Calibration->SetPointError(0, f1->GetParError(1), fGd->GetParError(1));
    // G_Calibration->SetPoint(1, f1->GetParameter(7), fPu->GetParameter(1));
    // G_Calibration->SetPointError(1, f1->GetParError(7), fPu->GetParError(1));
    // G_Calibration->SetPoint(2, f1->GetParameter(9), fPu->GetParameter(4));
    // G_Calibration->SetPointError(2, f1->GetParError(9), fPu->GetParError(4));
    // G_Calibration->SetPoint(3, f1->GetParameter(11), fPu->GetParameter(6));
    // G_Calibration->SetPointError(3, f1->GetParError(11), fPu->GetParError(6));
    G_Calibration->SetPoint(4-3, f1->GetParameter(13), fAm->GetParameter(1));
    G_Calibration->SetPointError(4-3, f1->GetParError(13), fAm->GetParError(1));
    G_Calibration->SetPoint(5-3, f1->GetParameter(15), fAm->GetParameter(4));
    G_Calibration->SetPointError(5-3, f1->GetParError(15), fAm->GetParError(4));
    G_Calibration->SetPoint(6-3, f1->GetParameter(17), fAm->GetParameter(6));
    G_Calibration->SetPointError(6-3, f1->GetParError(17), fAm->GetParError(6));
    G_Calibration->SetPoint(7-3, f1->GetParameter(19), fCm->GetParameter(1));
    G_Calibration->SetPointError(7-3, f1->GetParError(19), fCm->GetParError(1));
    G_Calibration->SetPoint(8-3, f1->GetParameter(21), fCm->GetParameter(4));
    G_Calibration->SetPointError(8-3, f1->GetParError(21), fCm->GetParError(4));


    // /////////////
    // //// 18N ////
    // /////////////
    // TH1D* H_Exp_18N = (TH1D*)GROUPED_File["32Ar"]->Get(("Strip/Strip_Channel/"+DET+"/Channel_Cleaned_"+DET).c_str());
    // H_Exp_18N->GetXaxis()->SetRangeUser(0, 20000);

    // peak_uncertainty = 500;
    // TF1 *f2BKG = new TF1("BKG_function", "( gaus(0) + gaus(3) + gaus(10) + [6]*exp([7] * x) ) * 0.5*(1+TMath::Erf((x-[8])/[9]))", 0, 20000);
    // f2BKG->SetNpx(1000000);
    // f2BKG->SetParameters(1e5, Peaks_Exp["18N"][0], 100, 5000, Peaks_Exp["18N"][1], 100, 1e5, -0.003, Peaks_Exp["18N"][0]-100, 5000);
    // f2BKG->SetParLimits(0, 5, 3000);
    // f2BKG->SetParLimits(1, Peaks_Exp["18N"][0]-peak_uncertainty, Peaks_Exp["18N"][0]+peak_uncertainty);
    // f2BKG->SetParLimits(2, 1, 500);
    // f2BKG->SetParLimits(3, 5, 3000);
    // f2BKG->SetParLimits(4, Peaks_Exp["18N"][1]-peak_uncertainty, Peaks_Exp["18N"][1]+peak_uncertainty);
    // f2BKG->SetParLimits(5, 1, 500);
    // f2BKG->SetParLimits(6, 1, 1e7);
    // f2BKG->SetParLimits(7, -0.001, 0);
    // f2BKG->SetParLimits(8, 4000, Peaks_Exp["18N"][0]);
    // f2BKG->SetParLimits(9, 150, 10000);
    // f2BKG->SetParameter(10, 5000);
    // f2BKG->SetParameter(11, 8000);
    // f2BKG->SetParLimits(12, 1, 100);

    // H_Exp_18N->Fit("BKG_function", "RL", "", 0, 20000);

    // TCanvas *c18N = new TCanvas("c18N", "c18N", 800, 600);
    // c18N->cd();
    // H_Exp_18N->Draw("HIST");
    // f2BKG->Draw("SAME");
    // c18N->Write();
    

    // TH1D* H_Sim_18N = (TH1D*)SIMULATED_File["18N"]->Get((DET+"_single").c_str());
    // H_Sim_18N->Add((TH1D*)SIMULATED_File["18N"]->Get("D6.1_single"));
    // H_Sim_18N->Add((TH1D*)SIMULATED_File["18N"]->Get("D7.1_single"));
    // H_Sim_18N->Add((TH1D*)SIMULATED_File["18N"]->Get("D8.1_single"));
    // H_Sim_18N->Rebin(10);

    // H_Sim_18N->GetXaxis()->SetRangeUser(Peaks_Sim["18N"][0]-200, Peaks_Sim["18N"][1]+200);

    // peak_uncertainty = 50;
    // TF1 *f18N = new TF1("f18N", "gaus(0) + gaus(3)", 0, 10000);
    // f18N->SetNpx(1000000);
    // f18N->SetParLimits(0, 10, 1e4);
    // f18N->SetParLimits(3, 10, 1e3);
    // f18N->SetParameter(1, Peaks_Sim["18N"][0]);
    // f18N->SetParLimits(1, Peaks_Sim["18N"][0]-peak_uncertainty, Peaks_Sim["18N"][0]+peak_uncertainty);
    // f18N->SetParameter(4, Peaks_Sim["18N"][1]);
    // f18N->SetParLimits(4, Peaks_Sim["18N"][1]-peak_uncertainty, Peaks_Sim["18N"][1]+peak_uncertainty);
    // f18N->SetParLimits(2, 10, 100);
    // f18N->SetParLimits(5, 10, 100);
    // H_Sim_18N->Fit("f18N", "RL", "", 0, 1400);

    // H_Sim_18N->GetXaxis()->SetRangeUser(Peaks_Sim["18N"][0]-100, Peaks_Sim["18N"][0]+100);
    // f18N->SetParameter(1, H_Sim_18N->GetMean());
    // f18N->SetParError(1, H_Sim_18N->GetRMS());
    // H_Sim_18N->GetXaxis()->SetRangeUser(Peaks_Sim["18N"][1]-100, Peaks_Sim["18N"][1]+100);
    // f18N->SetParameter(4, H_Sim_18N->GetMean());
    // f18N->SetParError(4, H_Sim_18N->GetRMS());




    // G_Calibration->SetPoint(1, f2BKG->GetParameter(1), f18N->GetParameter(1));
    // G_Calibration->SetPointError(1, f2BKG->GetParError(1), f18N->GetParError(1));
    // G_Calibration->SetPoint(2, f2BKG->GetParameter(4), f18N->GetParameter(4));
    // G_Calibration->SetPointError(2, f2BKG->GetParError(4), f18N->GetParError(4));




    // H_Sim_18N->Fit("f18N", "RL", "", 0, 1400);

    // TCanvas *cCalibration_18N = new TCanvas("cCalibration_18N", "cCalibration_18N", 800, 600);
    // cCalibration_18N->cd();
    // H_Sim_18N->Draw("HIST");
    // f18N->Draw("SAME");
    // cCalibration_18N->Write();

    // G_Calibration->SetPoint(9, f2BKG->GetParameter(1), f18N->GetParameter(1));
    // G_Calibration->SetPointError(9, f2BKG->GetParError(1), f18N->GetParError(1));
    // G_Calibration->SetPoint(10, f2BKG->GetParameter(4), f18N->GetParameter(4));
    // G_Calibration->SetPointError(10, f2BKG->GetParError(4), f18N->GetParError(4));



    ///////////////////////////////////////
    ///////// FITTING CALIBRATION //////////
    ///////////////////////////////////////

    TCanvas *cCalibrationGraph = new TCanvas("cCalibrationGraph", "cCalibrationGraph", 800, 600);
    cCalibrationGraph->Divide(1, 2);
    cCalibrationGraph->cd(1);
    G_Calibration->Draw("AP");
    G_Calibration->SetName("Calibration");
    G_Calibration->SetMarkerStyle(20);
    G_Calibration->SetMarkerSize(1);
    G_Calibration->SetMarkerColor(kBlack);
    G_Calibration->SetTitle("Calibration");
    TF1 *fit = new TF1("fit", "[0] + [1] *x", 0, 10000);
    fit->SetParameters(-20, 0.07);
    fit->SetParLimits(0, -1000, 1000);
    fit->SetParLimits(1, 0, 0.1);
    G_Calibration->Fit("fit");
    G_Calibration->GetFunction("fit")->Draw("SAME");
    cCalibrationGraph->cd(2);
    /// residus
    TGraphErrors *G_Residus = new TGraphErrors();
    for (int i = 0; i < G_Calibration->GetN(); i++)
    {
        double x, y;
        G_Calibration->GetPoint(i, x, y);
        G_Residus->SetPoint(i, x, y - G_Calibration->GetFunction("fit")->Eval(x));
        G_Residus->SetPointError(i, G_Calibration->GetErrorX(i), G_Calibration->GetErrorY(i));
    }
    G_Residus->Draw("AP");
    G_Residus->SetMarkerStyle(20);
    G_Residus->SetMarkerSize(1);
    G_Residus->SetMarkerColor(kBlack);
    TF1 *basleine = new TF1("basleine", "0", 0, 100000);
    basleine->Draw("SAME");
    basleine->SetLineColor(kRed);
    cCalibrationGraph->Write();  


    TTree* tree = (TTree*)GROUPED_File["4alpha"]->Get(("CLEANED_Tree_"+DET).c_str());
    TTreeReader *Reader = new TTreeReader(tree);
    TTreeReaderValue<double> Channel_Cleaned = { *Reader, "Channel" };


    TH1D* H_Calib = new TH1D("H_Calib", "H_Calib", 10000, 0, 10000);
    while(Reader->Next())
    {
        double Channel = *Channel_Cleaned;
        double Energy = fit->Eval(Channel);
        H_Calib->Fill(Energy);
    }


    H_Calib->Write();

    TCanvas *cCalib = new TCanvas("cCalib", "cCalib", 800, 600);
    cCalib->cd();
    H_Calib->Draw("HIST");
    cCalib->Write();

    TCanvas *cCalibGd = new TCanvas("cCalibGd", "cCalibGd", 800, 600);

    double sigma_g = 2.5;

    TF1 *fsource148Gd = new TF1("fsource148Gd", "1-exp(-0.12*x)", 0, 200);
    fsource148Gd->SetNpx(2000);
    TH1D* hsource148Gd = (TH1D*)fsource148Gd->GetHistogram();
    hsource148Gd->Write();

    TH1D* H_Sim_Gd_conv = new TH1D("H_Sim_Gd", "H_Sim_Gd", 10000, 0, 10000);
    for (int i = 0; i < H_Sim_Gd->GetEntries(); i++)
    {
        double value = H_Sim_Gd->GetRandom();
        uniform_real_distribution<double> distribution1(0, 1);
        double random = distribution1(generator);
        // get value in the pdf of fsource148Gd
        double x = 0;
        double y = 0;
        while (y < random)
        {
            x += 0.5;
            y = fsource148Gd->Eval(x);
        }
        normal_distribution<double> distribution(value-x, sigma_g);
        H_Sim_Gd_conv->Fill(distribution(generator));
    }

    TH1D* H_Sim_Am_conv = new TH1D("H_Sim_Am", "H_Sim_Am", 10000, 0, 10000);
    for (int i = 0; i < H_Sim_Am->GetEntries(); i++)
    {
        
        double value = H_Sim_Am->GetRandom();
        uniform_real_distribution<double> distribution1(0, 1);
        double random = distribution1(generator);
        double x = 0;
        double y = 0;
        while (y < random)
        {
            x += 0.5;
            y = fsource148Gd->Eval(x);
        }
        normal_distribution<double> distribution(value-x, sigma_g);
        H_Sim_Am_conv->Fill(distribution(generator));
    }

    TH1D* H_Sim_Cm_conv = new TH1D("H_Sim_Cm", "H_Sim_Cm", 10000, 0, 10000);
    for (int i = 0; i < H_Sim_Cm->GetEntries(); i++)
    {
        double value = H_Sim_Cm->GetRandom();
        uniform_real_distribution<double> distribution1(0, 1);
        double random = distribution1(generator);
        double x = 0;
        double y = 0;
        while (y < random)
        {
            x += 0.5;
            y = fsource148Gd->Eval(x);
        }
        normal_distribution<double> distribution(value-x, sigma_g);
        H_Sim_Cm_conv->Fill(distribution(generator));
    }

    TH1D* H_Sim_Pu_conv = new TH1D("H_Sim_Pu", "H_Sim_Pu", 10000, 0, 10000);
    for (int i = 0; i < H_Sim_Pu->GetEntries(); i++)
    {
        double value = H_Sim_Pu->GetRandom();
        uniform_real_distribution<double> distribution1(0, 1);
        double random = distribution1(generator);
        // get value in the pdf of fsource148Gd
        double x = 0;
        double y = 0;
        while (y < random)
        {
            x += 0.5;
            y = fsource148Gd->Eval(x);
        }
        normal_distribution<double> distribution(value-x, sigma_g);
        H_Sim_Pu_conv->Fill(distribution(generator));
    }


    cCalibGd->cd();
    

    H_Sim_Am_conv->GetXaxis()->SetRangeUser(Peaks_Sim["241Am"][0]-50, Peaks_Sim["241Am"][2]+50);
    H_Calib->GetXaxis()->SetRangeUser(Peaks_Sim["241Am"][0]-50, Peaks_Sim["241Am"][2]+50);
    H_Sim_Am_conv->Scale(H_Calib->GetMaximum()/H_Sim_Am_conv->GetMaximum());

    H_Sim_Cm_conv->GetXaxis()->SetRangeUser(Peaks_Sim["244Cm"][0]-50, Peaks_Sim["244Cm"][1]+50);
    H_Calib->GetXaxis()->SetRangeUser(Peaks_Sim["244Cm"][0]-50, Peaks_Sim["244Cm"][1]+50);
    H_Sim_Cm_conv->Scale(H_Calib->GetMaximum()/H_Sim_Cm_conv->GetMaximum());

    H_Sim_Pu_conv->GetXaxis()->SetRangeUser(Peaks_Sim["239Pu"][0]-50, Peaks_Sim["239Pu"][2]+50);
    H_Calib->GetXaxis()->SetRangeUser(Peaks_Sim["239Pu"][0]-50, Peaks_Sim["239Pu"][2]+50);
    H_Sim_Pu_conv->Scale(H_Calib->GetMaximum()/H_Sim_Pu_conv->GetMaximum());

    H_Sim_Gd_conv->GetXaxis()->SetRangeUser(Peaks_Sim["148Gd"][0]-50, Peaks_Sim["148Gd"][0]+50);
    H_Calib->GetXaxis()->SetRangeUser(Peaks_Sim["148Gd"][0]-50, Peaks_Sim["148Gd"][0]+50);
    H_Sim_Gd_conv->Scale(H_Calib->GetMaximum()/H_Sim_Gd_conv->GetMaximum());

    H_Calib->GetXaxis()->SetTitle("Energy [keV");
    H_Calib->GetYaxis()->SetTitle("Counts / keV");
    H_Calib->Draw("HIST");
    H_Sim_Gd_conv->Draw("HIST SAME");
    H_Sim_Am_conv->Draw("HIST SAME");
    H_Sim_Cm_conv->Draw("HIST SAME");
    H_Sim_Pu_conv->Draw("HIST SAME");
    H_Sim_Gd_conv->SetLineColor(kCyan);
    H_Sim_Am_conv->SetLineColor(kRed);
    H_Sim_Cm_conv->SetLineColor(kGreen);
    H_Sim_Pu_conv->SetLineColor(kBlack);

    TLegend *legendCalib = new TLegend(0.1, 0.7, 0.48, 0.9);
    legendCalib->AddEntry(H_Sim_Gd_conv, "148Gd", "l");
    legendCalib->AddEntry(H_Sim_Pu_conv, "239Pu", "l");
    legendCalib->AddEntry(H_Sim_Am_conv, "241Am", "l");
    legendCalib->AddEntry(H_Sim_Cm_conv, "244Cm", "l");
    legendCalib->Draw("SAME");
    cCalibGd->Write();

    



    // TF1 *calibration_p = new TF1("calibration_p", "[0] + [1] *x + [2] *x*x", 0, 10000);
    // calibration_p->SetParameters(-20.213525, 72.591352, -0.006507);
    // TH1D* H_ExpCalib_18N = new TH1D("H_ExpCalib_18N", "H_ExpCalib_18N", 10000, 0, 10000);


    // TCanvas *cCalib_18N = new TCanvas("cCalib_18N", "cCalib_18N", 800, 600);
    // cCalib_18N->cd();

    // TH1D* H_Calib_18N = new TH1D("H_Calib_18N", "H_Calib_18N", 10000, 0, 10000);
    // for (int i = 0; i < 1000000; i++)
    // {
    //     double value = H_Sim_18N->GetRandom();
    //     // double x = (-fit->GetParameter(1)+ sqrt(pow(fit->GetParameter(1), 2) - 4*fit->GetParameter(2) * (fit->GetParameter(0)-value) ))/(2*fit->GetParameter(2));
    //     double x = (value - fit->GetParameter(0))/fit->GetParameter(1);
    //     double energy = calibration_p->Eval(x/1000);
    //     H_Calib_18N->Fill(energy);
    // }

    // tree = (TTree*)GROUPED_File["32Ar"]->Get(("CLEANED_Tree_"+DET).c_str());
    // Reader = new TTreeReader(tree);
    // TTreeReaderValue<double> Channel_Cleaned1 = { *Reader, "Channel" };
    // Reader->Restart();
    // while (Reader->Next())
    // {
    //     double Channel = *Channel_Cleaned1;
    //     double Energy = calibration_p->Eval(Channel/1000);
    //     H_ExpCalib_18N->Fill(Energy);
    // }

    // H_ExpCalib_18N->Draw("HIST");
    // H_Calib_18N->Scale(7500/1600);
    // H_Calib_18N->Draw("HIST SAME");
    // H_Calib_18N->SetLineColor(kRed);

    // cCalib_18N->Write();


    DEFECT_File->Close();
}

//     /////////////////////////////     Am     ////////////////////////////////
//     min = 74000;
//     max = 78000;10);

//     halphas->GetXaxis()->SetRangeUser(min, max);
//     TF1 *f2 = new TF1("f2", tripleskew, min, max, 15);
//     f2->SetNpx(1000000);
//     f2->SetParameter(0, 75000);
//     f2->SetParLimits(0, min, 75500);
//     f2->SetParLimits(1, 1, 500);
//     f2->SetParLimits(2, -10, -2);
//     f2->SetParLimits(3, 0, 1000000);
//     f2->FixParameter(4, 0);
//     f2->SetParameter(5, 76000);
//     f2->SetParLimits(5, min, max);
//     f2->FixParameter(6, 0);
//     f2->FixParameter(7, 0);
//     f2->SetParLimits(8, 0, 1000000);
//     f2->FixParameter(9, 0);
//     f2->SetParameter(10, 76500);
//     f2->SetParLimits(10, min, max);
//     f2->FixParameter(11, 0);
//     f2->FixParameter(12, 0);
//     f2->SetParLimits(13, 0, 1000000);
//     f2->FixParameter(14, 0);

    // cout << "Fitting" << endl;
    // TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    // halphas->Draw("HIST");
    // f2->Draw("SAME");
    // c2->Write();

//     TF1 *f21 = new TF1("f21", skewedgauss, min, max, 5);
//     f21->SetNpx(1000000);
//     f21->SetParameter(0, f2->GetParameter(0));
//     f21->SetParameter(1, f2->GetParameter(1));
//     f21->SetParameter(2, f2->GetParameter(2));
//     f21->SetParameter(3, f2->GetParameter(3));
//     f21->SetParameter(4, f2->GetParameter(4));

//     TF1 *f22 = new TF1("f22", skewedgauss, min, max, 5);
//     f22->SetNpx(1000000);
//     f22->SetParameter(0, f2->GetParameter(5));
//     f22->SetParameter(1, f2->GetParameter(1));
//     f22->SetParameter(2, f2->GetParameter(2));
//     f22->SetParameter(3, f2->GetParameter(8));
//     f22->SetParameter(4, f2->GetParameter(9));

//     TF1 *f23 = new TF1("f23", skewedgauss, min, max, 5);
//     f23->SetNpx(1000000);
//     f23->SetParameter(0, f2->GetParameter(10));
//     f23->SetParameter(1, f2->GetParameter(1));
//     f23->SetParameter(2, f2->GetParameter(2));
//     f23->SetParameter(3, f2->GetParameter(13));
//     f23->SetParameter(4, f2->GetParameter(14));
    // TFitResultPtr r2 = halphas->Fit("f2", "RS", "", min, max);
    // TMatrixD c20 = r2->GetCovarianceMatrix();
//     TMatrixD c21 = c20.GetSub(0, 4, 0, 4);
//     TMatrixD c22 = c20.GetSub(5, 9, 5, 9);
//     TMatrixD c23 = c20.GetSub(10, 14, 10, 14);

//     double I1 = f21->Integral(min, max);
//     double I2 = f22->Integral(min, max);
//     double I3 = f23->Integral(min, max);

//     double dI1 = f21->IntegralError(min, max, f21->GetParameters(), c21.GetMatrixArray());
//     double dI2 = f22->IntegralError(min, max, f22->GetParameters(), c22.GetMatrixArray());
//     double dI3 = f23->IntegralError(min, max, f23->GetParameters(), c23.GetMatrixArray());

//     cout << 84.8 / I3 * I1;
//     cout << " +/- " << sqrt(pow(84.4 * I1 / pow(I3, 2) * dI3, 2) + pow(84.8 / I3 * dI1, 2)) << endl;

//     cout << 84.8 / I3 * I2;
//     cout << " +/- " << sqrt(pow(84.4 * I2. / pow(I3, 2) * dI3, 2) + pow(84.8 / I3 * dI2, 2)) << endl;
