

void Plotting_function()
{
    TFile *f = new TFile("Functions.root", "RECREATE");

    TCanvas *c = new TCanvas("","",800,600);
    TF1 *f0 = new TF1("f0","(x<0) ? 0 : 1",-5,5);
    f0->SetNpx(1000);
    f0->SetLineColor(kRed);
    f0->Draw();
    TLegend *legend0 = new TLegend(0.1,0.7,0.48,0.9);
    legend0->AddEntry(f0,"Step function","l");
    legend0->Draw();
    c->Write();

    TCanvas *c1 = new TCanvas("","",800,600);
    TF1 *f1 = new TF1("f1","TMath::Erf(x)",-5,5);
    f1->SetLineColor(kRed);
    f1->Draw();
    TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->AddEntry(f1,"erf(#frac{x}{#sigma})","l");
    legend->SetTextSize(0.1);
    legend->Draw();
    c1->Write();

    TCanvas *c2 = new TCanvas("","",800,600);
    TF1 *f2 = new TF1("f2","TMath::Erf((x+2)/0.5)-TMath::Erf((x-2)/0.5)",-5,5);
    f2->SetLineColor(kRed);
    f2->Draw();
    TLegend *legend2 = new TLegend(0.1,0.7,0.48,0.9);
    legend2->AddEntry(f2,"erf(#frac{x+l}{#sigma})-erf(#frac{x-l}{#sigma})","l");
    legend2->Draw();
    c2->Write();

    TCanvas *c3 = new TCanvas("","",800,600);
    TF2 *f3 = new TF2("f3","(TMath::Erf((x+2)/0.5)-TMath::Erf((x-2)/0.5)) * (TMath::Erf((y+2)/0.5)-TMath::Erf((y-2)/0.5))",-5,5,-5,5);
    f3->SetLineColor(kRed);
    f3->Draw("SURF2");
    TLegend *legend3 = new TLegend(0.1,0.7,0.48,0.9);
    legend3->AddEntry(f3,"(erf(#frac{x+l}{#sigma})-erf(#frac{x-l}{#sigma})*(erf(#frac{y+l}{#sigma})-erf(#frac{y-l}{#sigma}))","l");
    legend3->Draw();
    c3->Write();

    TCanvas *c4 = new TCanvas("", "", 800, 600);

    std::string total_function_expression;

    for (int i = 0 ; i < 5 ; i++)
    {
        for (int j = 0 ; j < 5 ; j++)
        {
            double x_mean = 10 * i - 25 ;
            double y_mean = 10 * j - 25;

            std::string function_expression = Form("(TMath::Erf((x-%f+2)/0.99)-TMath::Erf((x-%f-2)/0.99)) * (TMath::Erf((y-%f+2)/0.99)-TMath::Erf((y-%f-2)/0.99))", x_mean, x_mean, y_mean, y_mean);
            total_function_expression += function_expression + " + ";
        }
    }

    // Remove the last " + "
    total_function_expression = total_function_expression.substr(0, total_function_expression.size() - 3);

    TF2 *f_total = new TF2("f_total", total_function_expression.c_str(), -40, 40, -40, 40);
    f_total->SetNpx(600);
    f_total->SetNpy(600);
    f_total->SetLineColor(kRed);
    f_total->Draw("SURF2");
    
    TLegend *legend4 = new TLegend(0.1, 0.7, 0.48, 0.9);
    legend4->AddEntry(f_total, "Complete Grid", "l");
    legend4->Draw();
    c4->Write();

    f->Close();

    // TLegend *legend4 = new TLegend(0.1, 0.7, 0.48, 0.9);
    // legend4->AddEntry(f4, "Complete Grid", "l");
}