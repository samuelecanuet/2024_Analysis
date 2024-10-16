#include "../../Grouper/include/Detectors.hh"
#include "CLHEP/Vector/ThreeVector.h"

using namespace CLHEP;

int Verbosee = 0;

TFile *ANALYSIS_File;

TTreeReader *Reader;
TTreeReaderArray<int> *Tree_PDG;
TTreeReaderArray<double> *Tree_E0;
TTreeReaderArray<double> *Tree_x;
TTreeReaderArray<double> *Tree_y;
TTreeReaderArray<double> *Tree_z;
TTreeReaderArray<double> *Tree_Catcher_Central_Energy_Deposit;
TTreeReaderArray<double> *Tree_Catcher_Side_Energy_Deposit;
TTreeReaderArray<double> *Tree_PlasticScintillator_Energy_Deposit;
TTreeReaderArray<Hep3Vector> *Tree_PlasticScintillator_Hit_Position;
TTreeReaderArray<double> *Tree_PlasticScintillator_Hit_Angle;
TTreeReaderArray<vector<double>> *Tree_Silicon_Detector_Energy_Deposit;
TTreeReaderArray<vector<Hep3Vector>> *Tree_Silicon_Detector_Hit_Position;
TTreeReaderArray<vector<double>> *Tree_Silicon_Detector_Hit_Angle;
TTreeReaderArray<vector<int>> *Tree_Silicon_Detector_Code;
TTreeReaderArray<vector<double>> *Tree_Silicon_Detector_DL_Energy_Deposit;
TTreeReaderArray<vector<int>> *Tree_Silicon_Detector_InterStrip_Code;
vector<vector<vector<double>>> *Tree_Silicon_Detector_InterStrip_Energy_Deposit;
vector<vector<vector<Hep3Vector>>> *Tree_Silicon_Detector_InterStrip_Hit_Position;

TTree * OutTree[SIGNAL_MAX];
double e;
TTree * PlasticIASTree;
int sili_code;  
double sili_e;
double SiPM_e;

double size_interstrip = 0.07;
const int Particle_Size = 50;
TH1D *H_E0[Particle_Size];
TH1D *H_x[Particle_Size];
TH1D *H_y[Particle_Size];
TH1D *H_z[Particle_Size];
TH1D *H_Catcher_Central_Energy_Deposit[Particle_Size];
TH1D *H_Catcher_Side_Energy_Deposit[Particle_Size];
TH1D *H_PlasticScintillator_Energy_Deposit[Particle_Size];
TH1D *H_PlasticScintillator_Hit_Position_x[Particle_Size];
TH1D *H_PlasticScintillator_Hit_Position_y[Particle_Size];
TH1D *H_PlasticScintillator_Hit_Position_z[Particle_Size];
TH2D *H_PlasticScintillator_Hit_Position_xy[Particle_Size];
TH1D *H_PlasticScintillator_Hit_Angle[Particle_Size];
TH1D *H_Silicon_Detector_Energy_Deposit_Det[Particle_Size][SIGNAL_MAX];
TH1D *H_Silicon_Detector_Hit_Position_x[Particle_Size];
TH1D *H_Silicon_Detector_Hit_Position_y[Particle_Size];
TH1D *H_Silicon_Detector_Hit_Position_z[Particle_Size];
TH2D *H_Silicon_Detector_Hit_Position_xy[Particle_Size];
TH2D *H_Silicon_Detector_Hit_Position_xz[Particle_Size];
TH2D *H_Silicon_Detector_Hit_Position_yz[Particle_Size];
TH3D *H_Silicon_Detector_Hit_Position_xyz[Particle_Size];
TH1D *H_Silicon_Detector_Hit_Angle[Particle_Size];
TH1D *H_Silicon_Detector_Hit_Angle_Det[Particle_Size][SIGNAL_MAX];
TH2D *H_Silicon_Detector_Hit_Anglez[Particle_Size];
TH2D *H_Silicon_Detector_Hit_Anglexy[Particle_Size];
TH1D *H_Silicon_Detector_Code[Particle_Size];
TH1D *H_Silicon_Detector_DL_Energy_Deposit_Det[Particle_Size][SIGNAL_MAX];
TH2D *H_Silicon_Detector_InterStrip_xy[Particle_Size];
TH1D *H_Silicon_Detector_InterStrip_z[Particle_Size];
TH1D *H_Silicon_Detector_InterStrip_Energy_Deposit_Det[Particle_Size][100 * SIGNAL_MAX];

TDirectory *dir_Particle[Particle_Size];
TDirectory *dir_Initial[Particle_Size];
TDirectory *dir_Catcher[Particle_Size];
TDirectory *dir_PlasticScintillator[Particle_Size];
TDirectory *dir_Silicon_Detector[Particle_Size];
TDirectory *dir_Silicon_Detector_Det[Particle_Size][SIGNAL_MAX];
TDirectory *dir_Silicon_Detector_DL[Particle_Size];
TDirectory *dir_Silicon_Detector_InterStrip[Particle_Size];

double Plastic_Full_energy = 0;
double Full_energy[SIGNAL_MAX];
double InterStrip[SIGNAL_MAX];
vector<string> Particle_Used_String = {"32Ar", "32Cl", "31S", "31P", "33Ar", "33Cl", "32S", "18N",
                                       "239Pu", "148Gd", "244Cm", "241Am",
                                       "p", "n", "enu", "enubar", "e+", "e-", "gamma", "alpha", "4He"};

map<int, int> PDGtoIndex;

string SearchName(int PDG)
{
    for (const auto &pair : NametoCode_map)
    {
        if (pair.second == PDG)
        {
            return pair.first;
        }
    }

    Error("No name found for PDG: " + to_string(PDG));
}

void Init()
{
    gSystem->AddIncludePath("-I/softs/clhep/2.4.6.2/include/");
    gSystem->AddIncludePath("-I/softs/clhep/2.4.6.2/include/CLHEP/Vector/");
    gSystem->AddIncludePath("Merger/");
    Info("Initializing");
    for (int i = 0; i < Particle_Used_String.size(); i++)
    {
        PDGtoIndex[NametoCode_map.at(Particle_Used_String.at(i))] = i;
    }
}

void InitHistograms()
{
    Info("Initializing Histograms");
    for (int i = 0; i < Particle_Used_String.size(); i++)
    {
        string name = Particle_Used_String[i];

        H_E0[i] = new TH1D(("E0_" + name).c_str(), ("E0_" + name).c_str(), 15000, 0, 15000);
        H_E0[i]->GetXaxis()->SetTitle("Energy [keV]");
        H_E0[i]->GetYaxis()->SetTitle("Counts");
        H_E0[i]->GetXaxis()->CenterTitle();
        H_E0[i]->GetYaxis()->CenterTitle();

        H_x[i] = new TH1D(("x_" + name).c_str(), ("x_" + name).c_str(), 1000, -1000, 1000);
        H_x[i]->GetXaxis()->SetTitle("x [um]");
        H_x[i]->GetYaxis()->SetTitle("Counts");
        H_x[i]->GetXaxis()->CenterTitle();
        H_x[i]->GetYaxis()->CenterTitle();

        H_y[i] = new TH1D(("y_" + name).c_str(), ("y_" + name).c_str(), 1000, -1000, 1000);
        H_y[i]->GetXaxis()->SetTitle("y [um]");
        H_y[i]->GetYaxis()->SetTitle("Counts");
        H_y[i]->GetXaxis()->CenterTitle();
        H_y[i]->GetYaxis()->CenterTitle();

        H_z[i] = new TH1D(("z_" + name).c_str(), ("z_" + name).c_str(), 1000, -1000, 1000);
        H_z[i]->GetXaxis()->SetTitle("z [nm]");
        H_z[i]->GetYaxis()->SetTitle("Counts");
        H_z[i]->GetXaxis()->CenterTitle();
        H_z[i]->GetYaxis()->CenterTitle();

        H_Catcher_Central_Energy_Deposit[i] = new TH1D(("Catcher_Central_Energy_Deposit_" + name).c_str(), ("Catcher_Central_Energy_Deposit_" + name).c_str(), 15000, 0, 15000);
        H_Catcher_Central_Energy_Deposit[i]->GetXaxis()->SetTitle("Energy [keV]");
        H_Catcher_Central_Energy_Deposit[i]->GetYaxis()->SetTitle("Counts");
        H_Catcher_Central_Energy_Deposit[i]->GetXaxis()->CenterTitle();
        H_Catcher_Central_Energy_Deposit[i]->GetYaxis()->CenterTitle();

        H_Catcher_Side_Energy_Deposit[i] = new TH1D(("Catcher_Side_Energy_Deposit_" + name).c_str(), ("Catcher_Side_Energy_Deposit_" + name).c_str(), 15000, 0, 15000);
        H_Catcher_Side_Energy_Deposit[i]->GetXaxis()->SetTitle("Energy [keV]");
        H_Catcher_Side_Energy_Deposit[i]->GetYaxis()->SetTitle("Counts");
        H_Catcher_Side_Energy_Deposit[i]->GetXaxis()->CenterTitle();
        H_Catcher_Side_Energy_Deposit[i]->GetYaxis()->CenterTitle();

        H_PlasticScintillator_Energy_Deposit[i] = new TH1D(("PlasticScintillator_Energy_Deposit_" + name).c_str(), ("PlasticScintillator_Energy_Deposit_" + name).c_str(), 15000, 0, 15000);
        H_PlasticScintillator_Energy_Deposit[i]->GetXaxis()->SetTitle("Energy [keV]");
        H_PlasticScintillator_Energy_Deposit[i]->GetYaxis()->SetTitle("Counts");
        H_PlasticScintillator_Energy_Deposit[i]->GetXaxis()->CenterTitle();
        H_PlasticScintillator_Energy_Deposit[i]->GetYaxis()->CenterTitle();

        H_PlasticScintillator_Hit_Position_x[i] = new TH1D(("PlasticScintillator_Hit_Position_x_" + name).c_str(), ("PlasticScintillator_Hit_Position_x_" + name).c_str(), 1000, -50, 50);
        H_PlasticScintillator_Hit_Position_x[i]->GetXaxis()->SetTitle("X [mm]");
        H_PlasticScintillator_Hit_Position_x[i]->GetYaxis()->SetTitle("Counts");
        H_PlasticScintillator_Hit_Position_x[i]->GetXaxis()->CenterTitle();
        H_PlasticScintillator_Hit_Position_x[i]->GetYaxis()->CenterTitle();

        H_PlasticScintillator_Hit_Position_y[i] = new TH1D(("PlasticScintillator_Hit_Position_y_" + name).c_str(), ("PlasticScintillator_Hit_Position_y_" + name).c_str(), 1000, -50, 50);
        H_PlasticScintillator_Hit_Position_y[i]->GetXaxis()->SetTitle("Y [mm]");
        H_PlasticScintillator_Hit_Position_y[i]->GetYaxis()->SetTitle("Counts");
        H_PlasticScintillator_Hit_Position_y[i]->GetXaxis()->CenterTitle();
        H_PlasticScintillator_Hit_Position_y[i]->GetYaxis()->CenterTitle();

        H_PlasticScintillator_Hit_Position_z[i] = new TH1D(("PlasticScintillator_Hit_Position_z_" + name).c_str(), ("PlasticScintillator_Hit_Position_z_" + name).c_str(), 600, 0, 60);
        H_PlasticScintillator_Hit_Position_z[i]->GetXaxis()->SetTitle("Z [mm]");
        H_PlasticScintillator_Hit_Position_z[i]->GetYaxis()->SetTitle("Counts");
        H_PlasticScintillator_Hit_Position_z[i]->GetXaxis()->CenterTitle();
        H_PlasticScintillator_Hit_Position_z[i]->GetYaxis()->CenterTitle();

        H_PlasticScintillator_Hit_Position_xy[i] = new TH2D(("PlasticScintillator_Hit_Position_xy_" + name).c_str(), ("PlasticScintillator_Hit_Position_xy_" + name).c_str(), 400, -20, 20, 400, -20, 20);
        H_PlasticScintillator_Hit_Position_xy[i]->GetXaxis()->SetTitle("X [mm]");
        H_PlasticScintillator_Hit_Position_xy[i]->GetYaxis()->SetTitle("Y [mm]");
        H_PlasticScintillator_Hit_Position_xy[i]->GetXaxis()->CenterTitle();
        H_PlasticScintillator_Hit_Position_xy[i]->GetYaxis()->CenterTitle();

        H_PlasticScintillator_Hit_Angle[i] = new TH1D(("PlasticScintillator_Hit_Angle_" + name).c_str(), ("PlasticScintillator_Hit_Angle_" + name).c_str(), 900, 0, 90);
        H_PlasticScintillator_Hit_Angle[i]->GetXaxis()->SetTitle("Angle [deg]");
        H_PlasticScintillator_Hit_Angle[i]->GetYaxis()->SetTitle("Counts");
        H_PlasticScintillator_Hit_Angle[i]->GetXaxis()->CenterTitle();
        H_PlasticScintillator_Hit_Angle[i]->GetYaxis()->CenterTitle();

        H_Silicon_Detector_Hit_Position_x[i] = new TH1D(("Silicon_Detector_Hit_Position_x_" + name).c_str(), ("Silicon_Detector_Hit_Position_x_" + name).c_str(), 400, -20, 20);
        H_Silicon_Detector_Hit_Position_x[i]->GetXaxis()->SetTitle("X [mm]");
        H_Silicon_Detector_Hit_Position_x[i]->GetYaxis()->SetTitle("Counts");
        H_Silicon_Detector_Hit_Position_x[i]->GetXaxis()->CenterTitle();
        H_Silicon_Detector_Hit_Position_x[i]->GetYaxis()->CenterTitle();

        H_Silicon_Detector_Hit_Position_y[i] = new TH1D(("Silicon_Detector_Hit_Position_y_" + name).c_str(), ("Silicon_Detector_Hit_Position_y_" + name).c_str(), 400, -20, 20);
        H_Silicon_Detector_Hit_Position_y[i]->GetXaxis()->SetTitle("Y [mm]");
        H_Silicon_Detector_Hit_Position_y[i]->GetYaxis()->SetTitle("Counts");
        H_Silicon_Detector_Hit_Position_y[i]->GetXaxis()->CenterTitle();
        H_Silicon_Detector_Hit_Position_y[i]->GetYaxis()->CenterTitle();

        H_Silicon_Detector_Hit_Position_z[i] = new TH1D(("Silicon_Detector_Hit_Position_z_" + name).c_str(), ("Silicon_Detector_Hit_Position_z_" + name).c_str(), 1200, -60, 60);
        H_Silicon_Detector_Hit_Position_z[i]->GetXaxis()->SetTitle("Z [mm]");
        H_Silicon_Detector_Hit_Position_z[i]->GetYaxis()->SetTitle("Counts");
        H_Silicon_Detector_Hit_Position_z[i]->GetXaxis()->CenterTitle();
        H_Silicon_Detector_Hit_Position_z[i]->GetYaxis()->CenterTitle();

        H_Silicon_Detector_Hit_Position_xy[i] = new TH2D(("Silicon_Detector_Hit_Position_xy_" + name).c_str(), ("Silicon_Detector_Hit_Position_xy_" + name).c_str(), 1000, -50, 50, 1000, -50, 50);
        H_Silicon_Detector_Hit_Position_xy[i]->GetXaxis()->SetTitle("X [mm]");
        H_Silicon_Detector_Hit_Position_xy[i]->GetYaxis()->SetTitle("Y [mm]");
        H_Silicon_Detector_Hit_Position_xy[i]->GetXaxis()->CenterTitle();
        H_Silicon_Detector_Hit_Position_xy[i]->GetYaxis()->CenterTitle();

        H_Silicon_Detector_Hit_Position_xz[i] = new TH2D(("Silicon_Detector_Hit_Position_xz_" + name).c_str(), ("Silicon_Detector_Hit_Position_xz_" + name).c_str(), 1000, -50, 50, 1200, -60, 60);
        H_Silicon_Detector_Hit_Position_xz[i]->GetXaxis()->SetTitle("X [mm]");
        H_Silicon_Detector_Hit_Position_xz[i]->GetYaxis()->SetTitle("Z [mm]");
        H_Silicon_Detector_Hit_Position_xz[i]->GetXaxis()->CenterTitle();
        H_Silicon_Detector_Hit_Position_xz[i]->GetYaxis()->CenterTitle();

        H_Silicon_Detector_Hit_Position_yz[i] = new TH2D(("Silicon_Detector_Hit_Position_yz_" + name).c_str(), ("Silicon_Detector_Hit_Position_yz_" + name).c_str(), 1000, -50, 50, 1200, -60, 60);
        H_Silicon_Detector_Hit_Position_yz[i]->GetXaxis()->SetTitle("Y [mm]");
        H_Silicon_Detector_Hit_Position_yz[i]->GetYaxis()->SetTitle("Z [mm]");
        H_Silicon_Detector_Hit_Position_yz[i]->GetXaxis()->CenterTitle();
        H_Silicon_Detector_Hit_Position_yz[i]->GetYaxis()->CenterTitle();

        // H_Silicon_Detector_Hit_Position_xyz[i] = new TH3D(("Silicon_Detector_Hit_Position_xyz_" + name).c_str(), ("Silicon_Detector_Hit_Position_xyz_" + name).c_str(), 500, 0, 50, 500, 0, 50, 600, 0, 60);

        H_Silicon_Detector_Hit_Angle[i] = new TH1D(("Silicon_Detector_Hit_Angle_" + name).c_str(), ("Silicon_Detector_Hit_Angle_" + name).c_str(), 900, 0, 90);
        H_Silicon_Detector_Hit_Angle[i]->GetXaxis()->SetTitle("Angle [deg]");
        H_Silicon_Detector_Hit_Angle[i]->GetYaxis()->SetTitle("Counts");
        H_Silicon_Detector_Hit_Angle[i]->GetXaxis()->CenterTitle();
        H_Silicon_Detector_Hit_Angle[i]->GetYaxis()->CenterTitle();

        H_Silicon_Detector_Hit_Anglez[i] = new TH2D(("Silicon_Detector_Hit_Anglez_" + name).c_str(), ("Silicon_Detector_Hit_Anglez_" + name).c_str(), 900, 0, 90, 1200, -60, 60);
        H_Silicon_Detector_Hit_Anglez[i]->GetXaxis()->SetTitle("Angle [deg]");
        H_Silicon_Detector_Hit_Anglez[i]->GetYaxis()->SetTitle("Z [mm]");
        H_Silicon_Detector_Hit_Anglez[i]->GetXaxis()->CenterTitle();
        H_Silicon_Detector_Hit_Anglez[i]->GetYaxis()->CenterTitle();

        H_Silicon_Detector_Hit_Anglexy[i] = new TH2D(("Silicon_Detector_Hit_Anglexy_" + name).c_str(), ("Silicon_Detector_Hit_Anglexy_" + name).c_str(), 900, 0, 90, 1000, -50, 50);
        H_Silicon_Detector_Hit_Anglexy[i]->GetXaxis()->SetTitle("Angle [deg]");
        H_Silicon_Detector_Hit_Anglexy[i]->GetYaxis()->SetTitle("X [mm]");
        H_Silicon_Detector_Hit_Anglexy[i]->GetXaxis()->CenterTitle();
        H_Silicon_Detector_Hit_Anglexy[i]->GetYaxis()->CenterTitle();

        H_Silicon_Detector_Code[i] = new TH1D(("Silicon_Detector_Code_" + name).c_str(), ("Silicon_Detector_Code_" + name).c_str(), 100, 0, 100);
        H_Silicon_Detector_Code[i]->GetXaxis()->SetTitle("Code");
        H_Silicon_Detector_Code[i]->GetYaxis()->SetTitle("Counts");
        H_Silicon_Detector_Code[i]->GetXaxis()->CenterTitle();
        H_Silicon_Detector_Code[i]->GetYaxis()->CenterTitle();

        H_Silicon_Detector_InterStrip_xy[i] = new TH2D(("Silicon_Detector_InterStrip_xy_" + name).c_str(), ("Silicon_Detector_InterStrip_xy_" + name).c_str(), 8000, -40, 40, 1000, -0.05, 0.05);
        H_Silicon_Detector_InterStrip_xy[i]->GetXaxis()->SetTitle("X [mm]");
        H_Silicon_Detector_InterStrip_xy[i]->GetYaxis()->SetTitle("Y [mm]");
        H_Silicon_Detector_InterStrip_xy[i]->GetXaxis()->CenterTitle();
        H_Silicon_Detector_InterStrip_xy[i]->GetYaxis()->CenterTitle();

        H_Silicon_Detector_InterStrip_z[i] = new TH1D(("Silicon_Detector_InterStrip_z_" + name).c_str(), ("Silicon_Detector_InterStrip_z_" + name).c_str(), 4000, -0.2, 0.2);
        H_Silicon_Detector_InterStrip_z[i]->GetXaxis()->SetTitle("Z [um]");
        H_Silicon_Detector_InterStrip_z[i]->GetYaxis()->SetTitle("Counts");
        H_Silicon_Detector_InterStrip_z[i]->GetXaxis()->CenterTitle();
        H_Silicon_Detector_InterStrip_z[i]->GetYaxis()->CenterTitle();

        for (int det = 0; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorSiliStrip(det))
            {
                H_Silicon_Detector_Energy_Deposit_Det[i][det] = new TH1D(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_" + name).c_str(), ("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_" + name).c_str(), 10000, 0, 10000);
                H_Silicon_Detector_Energy_Deposit_Det[i][det]->GetXaxis()->SetTitle("Energy [keV]");
                H_Silicon_Detector_Energy_Deposit_Det[i][det]->GetYaxis()->SetTitle("Counts");
                H_Silicon_Detector_Energy_Deposit_Det[i][det]->GetXaxis()->CenterTitle();
                H_Silicon_Detector_Energy_Deposit_Det[i][det]->GetYaxis()->CenterTitle();

                H_Silicon_Detector_Hit_Angle_Det[i][det] = new TH1D(("Silicon_Detector_Hit_Angle_" + detectorName[det] + "_" + name).c_str(), ("Silicon_Detector_Hit_Angle_" + detectorName[det] + "_" + name).c_str(), 900, 0, 90);
                H_Silicon_Detector_Hit_Angle_Det[i][det]->GetXaxis()->SetTitle("Angle [deg]");
                H_Silicon_Detector_Hit_Angle_Det[i][det]->GetYaxis()->SetTitle("Counts");
                H_Silicon_Detector_Hit_Angle_Det[i][det]->GetXaxis()->CenterTitle();
                H_Silicon_Detector_Hit_Angle_Det[i][det]->GetYaxis()->CenterTitle();

                H_Silicon_Detector_DL_Energy_Deposit_Det[i][det] = new TH1D(("Silicon_Detector_DL_Energy_Deposit_" + detectorName[det] + "_" + name).c_str(), ("Silicon_Detector_DL_Energy_Deposit_" + detectorName[det] + "_" + name).c_str(), 10000, 0, 10000);
                H_Silicon_Detector_DL_Energy_Deposit_Det[i][det]->GetXaxis()->SetTitle("Energy [keV]");
                H_Silicon_Detector_DL_Energy_Deposit_Det[i][det]->GetYaxis()->SetTitle("Counts");
                H_Silicon_Detector_DL_Energy_Deposit_Det[i][det]->GetXaxis()->CenterTitle();
                H_Silicon_Detector_DL_Energy_Deposit_Det[i][det]->GetYaxis()->CenterTitle();

                if (GetDetectorChannel(det) != 5)
                {
                    int interstrip = GetDetector(det) * 1000 + ((GetDetectorChannel(det) * 2 + 1)) * 100 / 2;
                    H_Silicon_Detector_InterStrip_Energy_Deposit_Det[i][interstrip] = new TH1D(("Silicon_Detector_InterStrip_Energy_Deposit_" + detectorName[det] + to_string(GetDetectorChannel(det) + 1) + "_" + name).c_str(), ("Silicon_Detector_InterStrip_Energy_Deposit_" + detectorName[det] + to_string(GetDetectorChannel(det) + 1) + "_" + name).c_str(), 10000, 0, 10000);
                    H_Silicon_Detector_InterStrip_Energy_Deposit_Det[i][interstrip]->GetXaxis()->SetTitle("Energy [keV]");
                    H_Silicon_Detector_InterStrip_Energy_Deposit_Det[i][interstrip]->GetYaxis()->SetTitle("Counts");
                    H_Silicon_Detector_InterStrip_Energy_Deposit_Det[i][interstrip]->GetXaxis()->CenterTitle();
                    H_Silicon_Detector_InterStrip_Energy_Deposit_Det[i][interstrip]->GetYaxis()->CenterTitle();
                }
            }
        }
    }

    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det))
        {
            H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 1][det] = new TH1D(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All").c_str(), ("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All").c_str(), 10000, 0, 10000);
            H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 1][det]->GetXaxis()->SetTitle("Energy [keV]");
            H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 1][det]->GetYaxis()->SetTitle("Counts");
            H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 1][det]->GetXaxis()->CenterTitle();
            H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 1][det]->GetYaxis()->CenterTitle();

            H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 2][det] = new TH1D(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All_InterStrip").c_str(), ("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All_InterStrip").c_str(), 10000, 0, 10000);
            H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 2][det]->GetXaxis()->SetTitle("Energy [keV]");
            H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 2][det]->GetYaxis()->SetTitle("Counts");
            H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 2][det]->GetXaxis()->CenterTitle();
            H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 2][det]->GetYaxis()->CenterTitle();
        }
    }
}




void WriteHistograms()
{

    Info("Writing Histograms");
    ANALYSIS_File->cd();
    for (int i = 0; i < Particle_Used_String.size(); i++)
    {
        if (H_E0[i]->GetEntries() == 0)
        {
            continue;
        }

        string name = Particle_Used_String[i];

        dir_Particle[i] = ANALYSIS_File->mkdir(name.c_str());
        dir_Initial[i] = dir_Particle[i]->mkdir("Initial");
        dir_Catcher[i] = dir_Particle[i]->mkdir("Catcher");
        dir_PlasticScintillator[i] = dir_Particle[i]->mkdir("PlasticScintillator");
        dir_Silicon_Detector[i] = dir_Particle[i]->mkdir("Silicon_Detector");
        dir_Silicon_Detector_InterStrip[i] = dir_Particle[i]->mkdir("Silicon_Detector_InterStrip");
        for (int det = 0; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorSiliBack(det))
            {
                dir_Silicon_Detector_Det[i][GetDetector(det)] = dir_Silicon_Detector[i]->mkdir(("D." + to_string(GetDetector(det))).c_str());
            }
        }

        dir_Initial[i]->cd();
        H_E0[i]->Write();
        H_x[i]->Write();
        H_y[i]->Write();
        H_z[i]->Write();

        dir_Catcher[i]->cd();
        H_Catcher_Central_Energy_Deposit[i]->Write();
        H_Catcher_Side_Energy_Deposit[i]->Write();

        dir_PlasticScintillator[i]->cd();
        H_PlasticScintillator_Energy_Deposit[i]->Write();
        H_PlasticScintillator_Hit_Position_x[i]->Write();
        H_PlasticScintillator_Hit_Position_y[i]->Write();
        H_PlasticScintillator_Hit_Position_z[i]->Write();
        H_PlasticScintillator_Hit_Position_xy[i]->Write();
        H_PlasticScintillator_Hit_Angle[i]->Write();

        dir_Silicon_Detector[i]->cd();
        H_Silicon_Detector_Hit_Position_x[i]->Write();
        H_Silicon_Detector_Hit_Position_y[i]->Write();
        H_Silicon_Detector_Hit_Position_z[i]->Write();
        H_Silicon_Detector_Hit_Position_xy[i]->Write();
        H_Silicon_Detector_Hit_Position_xz[i]->Write();
        H_Silicon_Detector_Hit_Position_yz[i]->Write();
        H_Silicon_Detector_Hit_Angle[i]->Write();
        H_Silicon_Detector_Hit_Anglez[i]->Write();
        H_Silicon_Detector_Hit_Anglexy[i]->Draw("COLZ");
        H_Silicon_Detector_Hit_Anglexy[i]->Write();
        H_Silicon_Detector_Code[i]->Write();
        dir_Silicon_Detector_InterStrip[i]->cd();
        H_Silicon_Detector_InterStrip_xy[i]->Write();
        H_Silicon_Detector_InterStrip_z[i]->Write();

        for (int det = 0; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorSiliStrip(det))
            {
                dir_Silicon_Detector_Det[i][GetDetector(det)]->cd();
                H_Silicon_Detector_Energy_Deposit_Det[i][det]->Write();
                H_Silicon_Detector_Hit_Angle_Det[i][det]->Write();               
            }
        }
    }

    ANALYSIS_File->cd();
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det))
        {
            H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 1][det]->Write();

            TCanvas *c = new TCanvas(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All_WithWithoutInterStrip").c_str(), ("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All_WithWithoutInterStrip").c_str(), 800, 600);
            c->cd();
            H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 1][det]->SetLineColor(kBlack);
            H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 1][det]->Draw("HIST");
            H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 2][det]->SetLineColor(kRed);
            H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 2][det]->Draw("HIST SAME");
            c->Write();
            H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 2][det]->Write();


            TCanvas *c2 = new TCanvas(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_ParticlevsAll").c_str(), ("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_pvsAll").c_str(), 800, 600);
            c2->cd();
            TH1D* sub = (TH1D*)H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 1][det]->Clone();
            sub->Add(H_Silicon_Detector_Energy_Deposit_Det[Particle_Size - 2][det], -1);
            sub->SetLineColor(kBlack);
            sub->Draw("HIST");
            
            for (int i = 0; i < Particle_Used_String.size(); i++)
            {
                if (H_E0[i]->GetEntries() == 0)
                {
                    continue;
                }
                H_Silicon_Detector_Energy_Deposit_Det[i][det]->SetLineColor(kBlue + i*2);
                H_Silicon_Detector_Energy_Deposit_Det[i][det]->Draw("HIST SAME");
            }
            c2->Write();
        }
    }
}