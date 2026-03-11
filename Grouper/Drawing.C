#include "include/Detectors.hh"
#include "include/Utilities.hh"

TGraph *DATA_A;
TGraph *DATA_B;


int Drawing()
{


    FLAG2025 = true;
    InitDetectors("Config_Files/sample.pid");
    InitRuns();

    double guess_a, guess_b;

    TFile *f = MyTFile((DIR_ROOT_DATA_MATCHED + "SiPM_Matching_test_woONOFF.root").c_str(), "READ");

    TCanvas *c = (TCanvas *)f->Get("070/Graph_HighLowSiPM_Fitting_070");
    for (auto key : *c->GetListOfPrimitives())
    {
        TPad *pad = (TPad *)key;
        string name = pad->GetName();
        cout << name << endl;
        if (name.find("070_1") != string::npos)
        {
            for (auto key2 : *pad->GetListOfPrimitives())
            {
                if (string(key2->ClassName()) == "TF1")
                {
                    TF1 * fkey2 = (TF1 *)key2;
                    guess_a = fkey2->GetParameter(0);
                    guess_b = fkey2->GetParameter(1);
                }

                // if type  == TGRaphErrors
                if (string(key2->ClassName()) != "TGraph")
                    continue;
                TGraph * graph = (TGraph *)key2;
                // cRaw->cd();
                // graph->Draw("AP");
                DATA_A = (TGraph*)graph->Clone(); 
            }
            break;
        }
    }
    // TDirectory *dir = (TDirectory *)f->Get("070/SiPM_9");
    // TIter next(dir->GetListOfKeys());
    // TKey *key;
    // for (; (key = (TKey *)next());)
    // {
    //     string name = key->GetClassName();
    //     if (name.find("TGraph") != string::npos)
    //     {
    //         TGraph * graph = (TGraph *)key->ReadObj();
    //         DATA_A = (TGraph*)graph->Clone(); 
    //         break;
    //     }
    // }

    TFile *f1 = MyTFile((DIR_ROOT_DATA_MATCHED + "SiPM_Matching_test_woONOFF.root").c_str(), "READ");

    TCanvas *c1 = (TCanvas *)f1->Get("105/Graph_HighLowSiPM_Fitting_105");
    for (auto key : *c1->GetListOfPrimitives())
    {
        TPad *pad = (TPad *)key;
        string name = pad->GetName();
        cout << name << endl;
        if (name.find("105_1") != string::npos)
        {
            for (auto key2 : *pad->GetListOfPrimitives())
            {
                if (string(key2->ClassName()) == "TF1")
                {
                    TF1 * fkey2 = (TF1 *)key2;
                    guess_a = fkey2->GetParameter(0);
                    guess_b = fkey2->GetParameter(1);
                }

                // if type  == TGRaphErrors
                if (string(key2->ClassName()) != "TGraph")
                    continue;
                TGraph * graph = (TGraph *)key2;
                // cRaw->cd();
                // graph->Draw("AP");
                DATA_B = (TGraph*)graph->Clone(); 
            }
            break;
        }
    }

    // loop on keys in the folder 105/SiPM_9/ in f1
    // dir = (TDirectory *)f1->Get("105/SiPM_9");
    // TIter next1(dir->GetListOfKeys());
    // for (; (key = (TKey *)next1());)
    // {
    //     string name = key->GetClassName();
    //     if (name.find("TGraph") != string::npos)
    //     {
    //         TGraph * graph = (TGraph *)key->ReadObj();
    //         DATA_B = (TGraph*)graph->Clone(); 
    //         break;
    //     }
    // }



    TCanvas *cCompare = new TCanvas ("cCompare", "cCompare", 800, 600);
    DATA_A->SetMarkerColor(kRed);
    DATA_A->SetMarkerStyle(20);
    DATA_A->SetMarkerSize(0.2);
    DATA_B->SetMarkerColor(kBlue);
    DATA_B->SetMarkerStyle(20);
    DATA_B->SetMarkerSize(0.2);
    TMultiGraph *mg_compare = new TMultiGraph();
    mg_compare->Add(DATA_A, "AP");
    mg_compare->Add(DATA_B, "AP");
    mg_compare->Draw("AP");
    cCompare->Draw();





    return 0;
}