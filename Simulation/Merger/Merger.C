
void Merger()
{
  gSystem->AddIncludePath("-I/softs/clhep/2.4.6.2/include/");
  gSystem->AddIncludePath("-I/softs/clhep/2.4.6.2/include/CLHEP/Vector/");
  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<vector<double> > >", "vector");
  gInterpreter->GenerateDictionary("vector<CLHEP::Hep3Vector>", "vector;ThreeVector.h");
  gInterpreter->GenerateDictionary("vector<vector<CLHEP::Hep3Vector>>", "vector;ThreeVector.h");
  gInterpreter->GenerateDictionary("vector<vector<vector<CLHEP::Hep3Vector>>>", "vector;ThreeVector.h");

  // TString directory = "/mnt/hgfs/shared-2/2024_DATA/SIMULATED/DATA/";

  // TString nuclei = "32ArRMATRIX_";
  // TString param = "_CS0_CSP0_CV1_CVP1";
  // TString output_file = nuclei + "_" + param + ".root";

  // TChain chain("Tree");
  // for (int i = 1; i <= 5; i++)
  // { // replace 10 with the actual number of files
  //   TString file = nuclei + param + Form("_%d.root", i);
  //   if (!gSystem->AccessPathName(file))
  //   { // check if file exists
  //     TFile f(file);
  //     if (!f.IsZombie())
  //     {
  //       chain.Add(file);
  //     }
  //   }
  // }

  // chain.Merge(output_file);
}
