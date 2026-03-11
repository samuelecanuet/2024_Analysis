#include "Grouper/include/Detectors.hh"
#include "Grouper/include/Utilities.hh"

double Chi2(const TGraphErrors* g1, const TGraphErrors* g2) {
    if (!g1 || !g2) {
        Error("One of the graphs is null.");
        return -1;
    }

    int n1 = g1->GetN();
    int n2 = g2->GetN();

    if (n1 != n2) {
        Error("Graphs have different number of points.");
        return -1;
    }

    double chi2 = 0.0;
    for (int i = 0; i < n1; ++i) {
        double x1, y1, x2, y2;
        g1->GetPoint(i, x1, y1);
        g2->GetPoint(i, x2, y2);

        if (std::abs(x1 - x2) > 1e-6) {
            Error("x-values at point " + std::to_string(i) + " differ: " + std::to_string(x1) + " vs " + std::to_string(x2));   
            continue;
        }

        // if (x1 > 50)
        //     continue;

        double err1 = g1->GetErrorY(i);
        double err2 = g2->GetErrorY(i);
        double err_tot2 = err1 * err1 + err2 * err2;

        if (err_tot2 <= 0) {
            Error("Total error squared is zero or negative at point " + std::to_string(i));
            continue;
        }

        double delta = y1 - y2;
        double contrib = (delta * delta) / err_tot2;
        chi2 += contrib;
    }

    // if nan warning
    if (std::isnan(chi2)) {
        Warning("Chi2 calculation resulted in NaN.");
        return 1e6;
    }

    return chi2;
}

int EshiftCompareSimExpRelativePosition()
{
    FLAG2025 = true;
    InitDetectors("Grouper/Config_Files/sample.pid");


    double x0, y0, z0, angle0;

    if (YEAR == 2025)
    {
        x0 = -0.0536;
        y0 = 0.4103;
        z0 = 0.65;
        angle0 = 0.0;
    }
    else if (YEAR == 2024)
    {
        x0 = -0.1;
        y0 = 0.1;
        z0 = -0.5;
        angle0 = 0.0;
    }
    else
    {
        Error("YEAR not recognized");
    }

    // exp
    TFile *Exp = MyTFile((DIR_ROOT_DATA_ANALYSED + "32Ar_" + to_string(YEAR) + "_analysed_1.root").c_str(), "READ");    
    TCanvas *cEshiftCompare_Corrected = (TCanvas *)Exp->Get("Eshift/Eshift_M3/Eshift_3");
    TGraphErrors *G_Exp = nullptr;
    for (auto key : *cEshiftCompare_Corrected->GetListOfPrimitives())
    {
        if (string(key->ClassName()).find("TGraphErrors") != string::npos)
        {
            G_Exp = (TGraphErrors *)key;
        }
    }

    vector<string> filename;
    string path = "/run/media/local1/DATANEX/Samuel-G4/RelativePosition/";
    DIR* dir = opendir(path.c_str());
    if (dir == nullptr) {
        Error("Unable to open directory: " + path);
        return 1;
    }

    //analysed files
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        string file_name = entry->d_name;

        if (file_name == "." || file_name == "..") continue;

        if (file_name.find("result.root") != string::npos && file_name.find(to_string(YEAR)) != string::npos){
            TFile *file = new TFile((path + file_name).c_str(), "READ");
            if (file != nullptr)
            {
                filename.push_back(path+file_name);
                file->Close();
            }
        }
    }
    closedir(dir);

        map<string, TGraphErrors *>
            SimGraphs;
        map<string, double>
            chi2;

    map<double, TH2D*> Chi2DmapAngle;
    map<double, TH2D*> Chi2Dmapz;
    map<double, TGraph2D*> TGraphDmap;

    TCanvas *cEshiftCompare = new TCanvas("cEshiftCompare", "cEshiftCompare", 1200, 800);
    TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);

    int i = 0;
    for (auto f : filename)
    {
        // string name = Form("32Ar_%d_result_%ddeg.root", YEAR, i);
        TFile *Sim = MyTFile(f.c_str(), "READ");

        // sim
        SimGraphs[f] = nullptr;
        TCanvas *c_Sim = (TCanvas *)Sim->Get("Eshift");
        if (c_Sim == nullptr)
            continue;
        
        for (auto key : *c_Sim->GetListOfPrimitives())
        {
            if (string(key->ClassName()).find("TGraphErrors") != string::npos)
            {
                SimGraphs[f] = (TGraphErrors *)key;
            }
        }
       

        chi2[f] = Chi2(G_Exp, SimGraphs[f]);


        // regex to extract position from filename
        std::regex rgx("catcherz([-+]?[0-9]*\\.?[0-9]+)_catcher([-+]?[0-9]*\\.?[0-9]+)deg_x([-+]?[0-9]*\\.?[0-9]+)_y([-+]?[0-9]*\\.?[0-9]+)");
        std::smatch match;
        if (std::regex_search(f, match, rgx) && match.size() == 5) 
        {
            string z_pos = to_string(stod(match[1]) - z0);
            string angle = to_string(stod(match[2]) - angle0);
            string x_pos = to_string(stod(match[3]) - x0);
            string y_pos = to_string(stod(match[4]) - y0);

            if (stod(z_pos) == 0.0)
            {
                if (Chi2DmapAngle.find(stod(angle)) == Chi2DmapAngle.end())
                {
                    Chi2DmapAngle[stod(angle)] = new TH2D(Form("Chi2D_%s", angle.c_str()), Form("Chi2D Angle %s deg;X Position;Y Position;Chi2", angle.c_str()), 26, -6.5, 6.5, 53, -6.5, 6.5);
                }
                Chi2DmapAngle[stod(angle)]->Fill(stod(x_pos), stod(y_pos), chi2[f]);
            }
            if (stod(angle) == 0.0)
            {
                if (Chi2Dmapz.find(stod(z_pos)) == Chi2Dmapz.end())
                {
                    Chi2Dmapz[stod(z_pos)] = new TH2D(Form("Chi2D_z%s", z_pos.c_str()), Form("Chi2D z %s mm;X Position;Y Position;Chi2", z_pos.c_str()), 26, -6.5, 6.5, 53, -6.5, 6.5);
                }
                Chi2Dmapz[stod(z_pos)]->Fill(stod(x_pos), stod(y_pos), chi2[f]);
            }
            Info(Form("Angle: %.2f, X: %.2f, Y: %.2f, Z: %.2f, Chi2: %.2f", stod(angle), stod(x_pos), stod(y_pos), stod(z_pos), chi2[f]));
        }
        i++;
    }


    // printing 5 best chi2
    vector<pair<string, double>> chi2_vec(chi2.begin(), chi2.end());
    sort(chi2_vec.begin(), chi2_vec.end(), [](const pair<string, double>& a, const pair<string, double>& b) {
        return a.second < b.second;
    });
    Info("Top 5 Chi2 Results:");
    for (int j = 0; j < 5 && j < chi2_vec.size(); j++) {

        string f = chi2_vec[j].first;
        //regex to extract position from filename
        std::regex rgx("catcherz([-+]?[0-9]*\\.?[0-9]+)_catcher([-+]?[0-9]*\\.?[0-9]+)deg_x([-+]?[0-9]*\\.?[0-9]+)_y([-+]?[0-9]*\\.?[0-9]+)");
        std::smatch match;
        if (std::regex_search(f, match, rgx) && match.size() == 5) 
        {
            string z_pos = to_string(stod(match[1]) - z0);
            string angle = to_string(stod(match[2]) - angle0);
            string x_pos = to_string(stod(match[3]) - x0);
            string y_pos = to_string(stod(match[4]) - y0);

            Info(Form("%d) Angle: %.2f deg, X: %.2f mm, Y: %.2f mm, Z: %.2f mm, Chi2: %.2f", j+1, stod(angle), stod(x_pos), stod(y_pos), stod(z_pos), chi2_vec[j].second));
        }
    }

    i = 0;
    for (auto pair : chi2_vec) 
    {
        cEshiftCompare->cd();
        string f = pair.first;
        SimGraphs[f]->SetName(Form("G_Sim_%d", i));
        SimGraphs[f]->SetTitle(Form("Eshift Simulated %d", i));
        SimGraphs[f]->SetMarkerColor(i+2);
        SimGraphs[f]->SetLineColor(i+2);
        SimGraphs[f]->SetMarkerStyle(21);
        if (i == 0)
            SimGraphs[f]->Draw("AP");
        else
            SimGraphs[f]->Draw("P SAME");

        leg->AddEntry(SimGraphs[f], Form("Simulated %.0f", chi2[f]), "lp");

        i++;
        if (i >= 1) break;
    }


    cEshiftCompare->cd();
    G_Exp->SetMarkerColor(kBlack);
    G_Exp->SetLineColor(kBlack);
    G_Exp->SetMarkerStyle(20);
    G_Exp->SetMarkerSize(2);
    G_Exp->SetTitle("Eshift Comparison Simulated vs Experimental;Code;Eshift (keV)");
    G_Exp->Draw("P SAME");
    leg->AddEntry(G_Exp, "Experimental", "lp");
    leg->Draw("SAME");
    cEshiftCompare->Draw();

    // Chi2D ANGLE
    TCanvas *cChi2D = new TCanvas("cChi2DAngle_z=0.0mm", "cChi2DAngle_z=0.0mm", 1200, 800);
    int size = Chi2DmapAngle.size();
    double max = 0;
    double min = 1e9;
    for (auto pair : Chi2DmapAngle)
    {
        double local_max = pair.second->GetMaximum();
        if (local_max > max)
            max = local_max;
        // get min non zero
        double local_min = pair.second->GetMinimum(0.0001);
        if (local_min < min)
            min = local_min;
    }
    cChi2D->Divide(ceil(sqrt(size)), ceil(sqrt(size))); 
    for (auto pair : Chi2DmapAngle)
    {
        cChi2D->cd(++i);
        pair.second->SetMaximum(1000);
        pair.second->SetMinimum(min);
        pair.second->Draw("COLZ");
    }
    cChi2D->Draw();

    // Chi2D Z
    TCanvas *cChi2Dz = new TCanvas("cChi2Dz_Angle=0.0deg", "cChi2Dz_Angle=0.0deg", 1200, 800);
    size = Chi2Dmapz.size();
    max = 0;
    min = 1e9;  
    for (auto pair : Chi2Dmapz)
    {
        double local_max = pair.second->GetMaximum();
        if (local_max > max)
            max = local_max;
        // get min non zero
        double local_min = pair.second->GetMinimum(0.0001);
        if (local_min < min)
            min = local_min;
    }
    cChi2Dz->Divide(ceil(sqrt(size)), ceil(sqrt(size)));
    i = 0;
    for (auto pair : Chi2Dmapz)
    {
        cChi2Dz->cd(++i);
        pair.second->SetMaximum(1000);
        pair.second->SetMinimum(min);
        // displaying only unit ticks and labels
        pair.second->GetXaxis()->SetNdivisions(111);
        pair.second->GetYaxis()->SetNdivisions(11);
        pair.second->Draw("COLZ");
    }
    cChi2Dz->Draw();
    
    return 0;
}