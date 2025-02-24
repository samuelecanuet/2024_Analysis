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
TTreeReaderArray<double> *Tree_px;
TTreeReaderArray<double> *Tree_py;
TTreeReaderArray<double> *Tree_pz;
TTreeReaderArray<double> *Tree_T0;
TTreeReaderArray<double> *Tree_Catcher_Central_Energy_Deposit;
TTreeReaderArray<double> *Tree_Catcher_Side_Energy_Deposit;
TTreeReaderArray<double> *Tree_PlasticScintillator_Energy_Deposit;
TTreeReaderArray<Hep3Vector> *Tree_PlasticScintillator_Hit_Position;
TTreeReaderArray<double> *Tree_PlasticScintillator_Hit_Angle;
TTreeReaderArray<double> *Tree_PlasticScintillator_Hit_Time;
TTreeReaderArray<vector<double>> *Tree_Silicon_Detector_Energy_Deposit;
TTreeReaderArray<vector<Hep3Vector>> *Tree_Silicon_Detector_Hit_Position;
TTreeReaderArray<vector<double>> *Tree_Silicon_Detector_Hit_Angle;
TTreeReaderArray<vector<double>> *Tree_Silicon_Detector_Hit_Time;
TTreeReaderArray<vector<int>> *Tree_Silicon_Detector_Code;
TTreeReaderArray<vector<double>> *Tree_Silicon_Detector_DL_Energy_Deposit;
TTreeReaderArray<vector<int>> *Tree_Silicon_Detector_InterStrip_Code;
vector<vector<vector<double>>> *Tree_Silicon_Detector_InterStrip_Energy_Deposit;
vector<vector<vector<Hep3Vector>>> *Tree_Silicon_Detector_InterStrip_Hit_Position;

TTree *OutTree[SIGNAL_MAX];
double e;
TTree *PlasticIASTree;
int sili_code;
double sili_e;
double SiPM_e;

double size_interstrip = 0.07;

double time_e_IAS;

map<int, TH1D *> H_E0;
map<int, TH1D *> H_x;
map<int, TH1D *> H_y;
map<int, TH1D *> H_z;
map<int, TH1D *> H_px;
map<int, TH1D *> H_py;
map<int, TH1D *> H_pz;
map<int, TH1D *> H_T0;
map<int, TH1D *> H_Catcher_Central_Energy_Deposit;
map<int, TH1D *> H_Catcher_Side_Energy_Deposit;
map<int, TH1D *> H_PlasticScintillator_Energy_Deposit;
map<int, TH1D *> H_PlasticScintillator_Hit_Position_x;
map<int, TH1D *> H_PlasticScintillator_Hit_Position_y;
map<int, TH1D *> H_PlasticScintillator_Hit_Position_z;
map<int, TH2D *> H_PlasticScintillator_Hit_Position_xy;
map<int, TH1D *> H_PlasticScintillator_Hit_Angle;
map<int, TH1D *> H_PlasticScintillator_Hit_Time;
map<int, TH1D *[SIGNAL_MAX]> H_Silicon_Detector_Energy_Deposit_Det;
map<int, TH1D *[SIGNAL_MAX]> H_Silicon_Detector_Energy_Deposit_Det_without_interstrip;
map<int, TH2D *[SIGNAL_MAX]> H_Silicon_Detector_Energy_Deposit_Det_Rear;
map<int, TH1D *[SIGNAL_MAX]> H_Silicon_Detector_Energy_Deposit_SINGLE;
map<int, TH1D *[SIGNAL_MAX]> H_Silicon_Detector_Energy_Deposit_COINC;
map<int, TH1D *[SIGNAL_MAX]> H_Silicon_Detector_Energy_Deposit_NOCOINC;
map<int, TH1D *> H_Silicon_Detector_Hit_Position_x;
map<int, TH1D *> H_Silicon_Detector_Hit_Position_y;
map<int, TH1D *> H_Silicon_Detector_Hit_Position_z;
map<int, TH2D *> H_Silicon_Detector_Hit_Position_xy;
map<int, TH2D *> H_Silicon_Detector_Hit_Position_xz;
map<int, TH2D *> H_Silicon_Detector_Hit_Position_yz;
map<int, TH3D *> H_Silicon_Detector_Hit_Position_xyz;
map<int, TH1D *> H_Silicon_Detector_Hit_Angle;
map<int, TH1D *> H_Silicon_Detector_Hit_Time;
map<int, TH1D *[SIGNAL_MAX]> H_Silicon_Detector_Hit_Angle_Det;
map<int, TH2D *> H_Silicon_Detector_Hit_Anglez;
map<int, TH2D *> H_Silicon_Detector_Hit_Anglexy;
map<int, TH1D *> H_Silicon_Detector_Code;
map<int, TH1D *[SIGNAL_MAX]> H_Silicon_Detector_DL_Energy_Deposit_Det;
map<int, TH2D *> H_Silicon_Detector_InterStrip_xy;
map<int, TH1D *> H_Silicon_Detector_InterStrip_z;
map<int, TH1D * [100 * SIGNAL_MAX]> H_Silicon_Detector_InterStrip_Energy_Deposit_Det;
TH1D *H = new TH1D("H", "H", 10000, 0, 10000);

map<int, TDirectory *> dir_Particle;
map<int, TDirectory *> dir_Initial;
map<int, TDirectory *> dir_Catcher;
map<int, TDirectory *> dir_PlasticScintillator;
map<int, TDirectory *> dir_Silicon_Detector;
map<int, array<TDirectory *, SIGNAL_MAX>> dir_Silicon_Detector_Det;
map<int, TDirectory *> dir_Silicon_Detector_DL;
map<int, TDirectory *> dir_Silicon_Detector_InterStrip;

// TH2D* H[SIGNAL_MAX];
double Plastic_Full_energy = 0;
double Full_energy[SIGNAL_MAX];
double Full_energy_without_interstrip[SIGNAL_MAX];
double RearEnergy[SIGNAL_MAX];
double proton_energy[SIGNAL_MAX];
double positron_energy[SIGNAL_MAX];
double InterStrip[SIGNAL_MAX];
vector<string> Particle_Used_String = {"32Ar", "32Cl", "31S", "31P", "33Ar", "33Cl", "32S", "18N",
                                       "239Pu", "148Gd", "244Cm", "241Am",
                                       "p", "n", "enu", "enubar", "e+", "e-", "gamma", "alpha", "4He"};

map<int, int> PDGtoIndex;

double offset[SIGNAL_MAX];
string SearchName(int PDG)
{
    if (PDG == 0)
    {
        return "All";
    }

    if (PDG > 1000000000)
    {
        PDG = (PDG/10)*10;
    }
    for (const auto &pair : NametoCode_map)
    {
        if (pair.second == PDG)
        {
            return pair.first;
        }
    }

    Error("No name found for PDG: " + to_string(PDG));
    return "";
}

double Detector_Resolution[SIGNAL_MAX];

void InitElectronicResolution()
{
    ifstream file("../Grouper/Config_Files/res_electronic.txt");
    if (!file.is_open())
    {
        Error("Impossible to open res_electronic.txt");
    }
    Info("Reading res_electronic.txt");

    string line;
    int code;
    double resolution;
    double err_resolution;
    int counter = 0;
    while (getline(file, line))
    {
        counter++;
        stringstream ss(line);
        ss >> code >> resolution >> err_resolution;

        Detector_Resolution[code] = resolution / 1000.;
    }

    file.close();
}

TF1 *Calibration_Function[SIGNAL_MAX];
void InitCalib()
{
    TFile *CALIBRATED_File = MyTFile(DIR_ROOT_DATA_CALIBRATED + "Calibrated.root", "READ");
    for (int i = 0; i < SIGNAL_MAX; i++)
        {
            if (IsDetectorSiliStrip(i))
            {
                Calibration_Function[i] = (TF1 *)CALIBRATED_File->Get(("Calibration_" + detectorName[i]).c_str());

                if (Calibration_Function[i] == NULL)
                {
                    Calibration_Function[i] = new TF1(("Calibration_" + detectorName[i]).c_str(), "x", 0, 10000);
                    Warning("No calibration found for " + detectorName[i]);
                }
            }
        }
}

void Init()
{
    gSystem->AddIncludePath("-I/softs/clhep/2.4.6.2/include/");
    gSystem->AddIncludePath("-I/softs/clhep/2.4.6.2/include/CLHEP/Vector/");
    gSystem->AddIncludePath("Merger/");
    Info("Initialisation");
    for (int i = 0; i < Particle_Used_String.size(); i++)
    {
        PDGtoIndex[NametoCode_map.at(Particle_Used_String.at(i))] = i;
    }
}

void InitHistograms(int PDG_code)
{

    string name = SearchName(PDG_code);

    H_E0[PDG_code] = new TH1D(("E0_" + name).c_str(), ("E0_" + name).c_str(), 15000, 0, 15000);
    H_E0[PDG_code]->GetXaxis()->SetTitle("Energy [keV]");
    H_E0[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_E0[PDG_code]->GetXaxis()->CenterTitle();
    H_E0[PDG_code]->GetYaxis()->CenterTitle();

    H_x[PDG_code] = new TH1D(("x_" + name).c_str(), ("x_" + name).c_str(), 1000, -1000, 1000);
    H_x[PDG_code]->GetXaxis()->SetTitle("x [um]");
    H_x[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_x[PDG_code]->GetXaxis()->CenterTitle();
    H_x[PDG_code]->GetYaxis()->CenterTitle();

    H_y[PDG_code] = new TH1D(("y_" + name).c_str(), ("y_" + name).c_str(), 1000, -1000, 1000);
    H_y[PDG_code]->GetXaxis()->SetTitle("y [um]");
    H_y[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_y[PDG_code]->GetXaxis()->CenterTitle();
    H_y[PDG_code]->GetYaxis()->CenterTitle();

    H_z[PDG_code] = new TH1D(("z_" + name).c_str(), ("z_" + name).c_str(), 1000, -1000, 1000);
    H_z[PDG_code]->GetXaxis()->SetTitle("z [nm]");
    H_z[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_z[PDG_code]->GetXaxis()->CenterTitle();
    H_z[PDG_code]->GetYaxis()->CenterTitle();

    H_px[PDG_code] = new TH1D(("px_" + name).c_str(), ("px_" + name).c_str(), 1000, -1, 1);
    H_px[PDG_code]->GetXaxis()->SetTitle("px");
    H_px[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_px[PDG_code]->GetXaxis()->CenterTitle();
    H_px[PDG_code]->GetYaxis()->CenterTitle();

    H_py[PDG_code] = new TH1D(("py_" + name).c_str(), ("py_" + name).c_str(), 1000, -1, 1);
    H_py[PDG_code]->GetXaxis()->SetTitle("py");
    H_py[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_py[PDG_code]->GetXaxis()->CenterTitle();
    H_py[PDG_code]->GetYaxis()->CenterTitle();

    H_pz[PDG_code] = new TH1D(("pz_" + name).c_str(), ("pz_" + name).c_str(), 1000, -1, 1);
    H_pz[PDG_code]->GetXaxis()->SetTitle("pz");
    H_pz[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_pz[PDG_code]->GetXaxis()->CenterTitle();
    H_pz[PDG_code]->GetYaxis()->CenterTitle();

    H_T0[PDG_code] = new TH1D(("T0_" + name).c_str(), ("T0_" + name).c_str(), 1000, 0, 1000);
    H_T0[PDG_code]->GetXaxis()->SetTitle("#T_0 [ns]");
    H_T0[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_T0[PDG_code]->GetXaxis()->CenterTitle();
    H_T0[PDG_code]->GetYaxis()->CenterTitle();

    H_Catcher_Central_Energy_Deposit[PDG_code] = new TH1D(("Catcher_Central_Energy_Deposit_" + name).c_str(), ("Catcher_Central_Energy_Deposit_" + name).c_str(), 15000, 0, 15000);
    H_Catcher_Central_Energy_Deposit[PDG_code]->GetXaxis()->SetTitle("Energy [keV]");
    H_Catcher_Central_Energy_Deposit[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_Catcher_Central_Energy_Deposit[PDG_code]->GetXaxis()->CenterTitle();
    H_Catcher_Central_Energy_Deposit[PDG_code]->GetYaxis()->CenterTitle();

    H_Catcher_Side_Energy_Deposit[PDG_code] = new TH1D(("Catcher_Side_Energy_Deposit_" + name).c_str(), ("Catcher_Side_Energy_Deposit_" + name).c_str(), 15000, 0, 15000);
    H_Catcher_Side_Energy_Deposit[PDG_code]->GetXaxis()->SetTitle("Energy [keV]");
    H_Catcher_Side_Energy_Deposit[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_Catcher_Side_Energy_Deposit[PDG_code]->GetXaxis()->CenterTitle();
    H_Catcher_Side_Energy_Deposit[PDG_code]->GetYaxis()->CenterTitle();

    H_PlasticScintillator_Energy_Deposit[PDG_code] = new TH1D(("PlasticScintillator_Energy_Deposit_" + name).c_str(), ("PlasticScintillator_Energy_Deposit_" + name).c_str(), 15000, 0, 15000);
    H_PlasticScintillator_Energy_Deposit[PDG_code]->GetXaxis()->SetTitle("Energy [keV]");
    H_PlasticScintillator_Energy_Deposit[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_PlasticScintillator_Energy_Deposit[PDG_code]->GetXaxis()->CenterTitle();
    H_PlasticScintillator_Energy_Deposit[PDG_code]->GetYaxis()->CenterTitle();

    H_PlasticScintillator_Hit_Position_x[PDG_code] = new TH1D(("PlasticScintillator_Hit_Position_x_" + name).c_str(), ("PlasticScintillator_Hit_Position_x_" + name).c_str(), 1000, -50, 50);
    H_PlasticScintillator_Hit_Position_x[PDG_code]->GetXaxis()->SetTitle("X [mm]");
    H_PlasticScintillator_Hit_Position_x[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_PlasticScintillator_Hit_Position_x[PDG_code]->GetXaxis()->CenterTitle();
    H_PlasticScintillator_Hit_Position_x[PDG_code]->GetYaxis()->CenterTitle();

    H_PlasticScintillator_Hit_Position_y[PDG_code] = new TH1D(("PlasticScintillator_Hit_Position_y_" + name).c_str(), ("PlasticScintillator_Hit_Position_y_" + name).c_str(), 1000, -50, 50);
    H_PlasticScintillator_Hit_Position_y[PDG_code]->GetXaxis()->SetTitle("Y [mm]");
    H_PlasticScintillator_Hit_Position_y[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_PlasticScintillator_Hit_Position_y[PDG_code]->GetXaxis()->CenterTitle();
    H_PlasticScintillator_Hit_Position_y[PDG_code]->GetYaxis()->CenterTitle();

    H_PlasticScintillator_Hit_Position_z[PDG_code] = new TH1D(("PlasticScintillator_Hit_Position_z_" + name).c_str(), ("PlasticScintillator_Hit_Position_z_" + name).c_str(), 600, 0, 60);
    H_PlasticScintillator_Hit_Position_z[PDG_code]->GetXaxis()->SetTitle("Z [mm]");
    H_PlasticScintillator_Hit_Position_z[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_PlasticScintillator_Hit_Position_z[PDG_code]->GetXaxis()->CenterTitle();
    H_PlasticScintillator_Hit_Position_z[PDG_code]->GetYaxis()->CenterTitle();

    H_PlasticScintillator_Hit_Position_xy[PDG_code] = new TH2D(("PlasticScintillator_Hit_Position_xy_" + name).c_str(), ("PlasticScintillator_Hit_Position_xy_" + name).c_str(), 400, -20, 20, 400, -20, 20);
    H_PlasticScintillator_Hit_Position_xy[PDG_code]->GetXaxis()->SetTitle("X [mm]");
    H_PlasticScintillator_Hit_Position_xy[PDG_code]->GetYaxis()->SetTitle("Y [mm]");
    H_PlasticScintillator_Hit_Position_xy[PDG_code]->GetXaxis()->CenterTitle();
    H_PlasticScintillator_Hit_Position_xy[PDG_code]->GetYaxis()->CenterTitle();

    H_PlasticScintillator_Hit_Angle[PDG_code] = new TH1D(("PlasticScintillator_Hit_Angle_" + name).c_str(), ("PlasticScintillator_Hit_Angle_" + name).c_str(), 900, 0, 90);
    H_PlasticScintillator_Hit_Angle[PDG_code]->GetXaxis()->SetTitle("Angle [deg]");
    H_PlasticScintillator_Hit_Angle[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_PlasticScintillator_Hit_Angle[PDG_code]->GetXaxis()->CenterTitle();
    H_PlasticScintillator_Hit_Angle[PDG_code]->GetYaxis()->CenterTitle();

    H_PlasticScintillator_Hit_Time[PDG_code] = new TH1D(("PlasticScintillator_Time_" + name).c_str(), ("PlasticScintillator_Time_" + name).c_str(), 1000, 0, 1000);
    H_PlasticScintillator_Hit_Time[PDG_code]->GetXaxis()->SetTitle("Time [ns]");
    H_PlasticScintillator_Hit_Time[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_PlasticScintillator_Hit_Time[PDG_code]->GetXaxis()->CenterTitle();
    H_PlasticScintillator_Hit_Time[PDG_code]->GetYaxis()->CenterTitle();

    H_Silicon_Detector_Hit_Position_x[PDG_code] = new TH1D(("Silicon_Detector_Hit_Position_x_" + name).c_str(), ("Silicon_Detector_Hit_Position_x_" + name).c_str(), 400, -20, 20);
    H_Silicon_Detector_Hit_Position_x[PDG_code]->GetXaxis()->SetTitle("X [mm]");
    H_Silicon_Detector_Hit_Position_x[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_Silicon_Detector_Hit_Position_x[PDG_code]->GetXaxis()->CenterTitle();
    H_Silicon_Detector_Hit_Position_x[PDG_code]->GetYaxis()->CenterTitle();

    H_Silicon_Detector_Hit_Position_y[PDG_code] = new TH1D(("Silicon_Detector_Hit_Position_y_" + name).c_str(), ("Silicon_Detector_Hit_Position_y_" + name).c_str(), 400, -20, 20);
    H_Silicon_Detector_Hit_Position_y[PDG_code]->GetXaxis()->SetTitle("Y [mm]");
    H_Silicon_Detector_Hit_Position_y[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_Silicon_Detector_Hit_Position_y[PDG_code]->GetXaxis()->CenterTitle();
    H_Silicon_Detector_Hit_Position_y[PDG_code]->GetYaxis()->CenterTitle();

    H_Silicon_Detector_Hit_Position_z[PDG_code] = new TH1D(("Silicon_Detector_Hit_Position_z_" + name).c_str(), ("Silicon_Detector_Hit_Position_z_" + name).c_str(), 1200, -60, 60);
    H_Silicon_Detector_Hit_Position_z[PDG_code]->GetXaxis()->SetTitle("Z [mm]");
    H_Silicon_Detector_Hit_Position_z[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_Silicon_Detector_Hit_Position_z[PDG_code]->GetXaxis()->CenterTitle();
    H_Silicon_Detector_Hit_Position_z[PDG_code]->GetYaxis()->CenterTitle();

    H_Silicon_Detector_Hit_Position_xy[PDG_code] = new TH2D(("Silicon_Detector_Hit_Position_xy_" + name).c_str(), ("Silicon_Detector_Hit_Position_xy_" + name).c_str(), 1000, -50, 50, 1000, -50, 50);
    H_Silicon_Detector_Hit_Position_xy[PDG_code]->GetXaxis()->SetTitle("X [mm]");
    H_Silicon_Detector_Hit_Position_xy[PDG_code]->GetYaxis()->SetTitle("Y [mm]");
    H_Silicon_Detector_Hit_Position_xy[PDG_code]->GetXaxis()->CenterTitle();
    H_Silicon_Detector_Hit_Position_xy[PDG_code]->GetYaxis()->CenterTitle();

    H_Silicon_Detector_Hit_Position_xz[PDG_code] = new TH2D(("Silicon_Detector_Hit_Position_xz_" + name).c_str(), ("Silicon_Detector_Hit_Position_xz_" + name).c_str(), 1000, -50, 50, 1200, -60, 60);
    H_Silicon_Detector_Hit_Position_xz[PDG_code]->GetXaxis()->SetTitle("X [mm]");
    H_Silicon_Detector_Hit_Position_xz[PDG_code]->GetYaxis()->SetTitle("Z [mm]");
    H_Silicon_Detector_Hit_Position_xz[PDG_code]->GetXaxis()->CenterTitle();
    H_Silicon_Detector_Hit_Position_xz[PDG_code]->GetYaxis()->CenterTitle();

    H_Silicon_Detector_Hit_Position_yz[PDG_code] = new TH2D(("Silicon_Detector_Hit_Position_yz_" + name).c_str(), ("Silicon_Detector_Hit_Position_yz_" + name).c_str(), 1000, -50, 50, 1200, -60, 60);
    H_Silicon_Detector_Hit_Position_yz[PDG_code]->GetXaxis()->SetTitle("Y [mm]");
    H_Silicon_Detector_Hit_Position_yz[PDG_code]->GetYaxis()->SetTitle("Z [mm]");
    H_Silicon_Detector_Hit_Position_yz[PDG_code]->GetXaxis()->CenterTitle();
    H_Silicon_Detector_Hit_Position_yz[PDG_code]->GetYaxis()->CenterTitle();

    // H_Silicon_Detector_Hit_Position_xyz[PDG_code] = new TH3D(("Silicon_Detector_Hit_Position_xyz_" + name).c_str(), ("Silicon_Detector_Hit_Position_xyz_" + name).c_str(), 500, 0, 50, 500, 0, 50, 600, 0, 60);

    H_Silicon_Detector_Hit_Angle[PDG_code] = new TH1D(("Silicon_Detector_Hit_Angle_" + name).c_str(), ("Silicon_Detector_Hit_Angle_" + name).c_str(), 900, 0, 90);
    H_Silicon_Detector_Hit_Angle[PDG_code]->GetXaxis()->SetTitle("Angle [deg]");
    H_Silicon_Detector_Hit_Angle[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_Silicon_Detector_Hit_Angle[PDG_code]->GetXaxis()->CenterTitle();
    H_Silicon_Detector_Hit_Angle[PDG_code]->GetYaxis()->CenterTitle();

    H_Silicon_Detector_Hit_Anglez[PDG_code] = new TH2D(("Silicon_Detector_Hit_Anglez_" + name).c_str(), ("Silicon_Detector_Hit_Anglez_" + name).c_str(), 900, 0, 90, 1200, -60, 60);
    H_Silicon_Detector_Hit_Anglez[PDG_code]->GetXaxis()->SetTitle("Angle [deg]");
    H_Silicon_Detector_Hit_Anglez[PDG_code]->GetYaxis()->SetTitle("Z [mm]");
    H_Silicon_Detector_Hit_Anglez[PDG_code]->GetXaxis()->CenterTitle();
    H_Silicon_Detector_Hit_Anglez[PDG_code]->GetYaxis()->CenterTitle();

    H_Silicon_Detector_Hit_Anglexy[PDG_code] = new TH2D(("Silicon_Detector_Hit_Anglexy_" + name).c_str(), ("Silicon_Detector_Hit_Anglexy_" + name).c_str(), 900, 0, 90, 1000, -50, 50);
    H_Silicon_Detector_Hit_Anglexy[PDG_code]->GetXaxis()->SetTitle("Angle [deg]");
    H_Silicon_Detector_Hit_Anglexy[PDG_code]->GetYaxis()->SetTitle("X [mm]");
    H_Silicon_Detector_Hit_Anglexy[PDG_code]->GetXaxis()->CenterTitle();
    H_Silicon_Detector_Hit_Anglexy[PDG_code]->GetYaxis()->CenterTitle();

    H_Silicon_Detector_Hit_Time[PDG_code] = new TH1D(("Silicon_Detector_Time_" + name).c_str(), ("Silicon_Detector_Time_" + name).c_str(), 1000, 0, 1000);
    H_Silicon_Detector_Hit_Time[PDG_code]->GetXaxis()->SetTitle("Time [ns]");
    H_Silicon_Detector_Hit_Time[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_Silicon_Detector_Hit_Time[PDG_code]->GetXaxis()->CenterTitle();
    H_Silicon_Detector_Hit_Time[PDG_code]->GetYaxis()->CenterTitle();

    H_Silicon_Detector_Code[PDG_code] = new TH1D(("Silicon_Detector_Code_" + name).c_str(), ("Silicon_Detector_Code_" + name).c_str(), 100, 0, 100);
    H_Silicon_Detector_Code[PDG_code]->GetXaxis()->SetTitle("Code");
    H_Silicon_Detector_Code[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_Silicon_Detector_Code[PDG_code]->GetXaxis()->CenterTitle();
    H_Silicon_Detector_Code[PDG_code]->GetYaxis()->CenterTitle();

    H_Silicon_Detector_InterStrip_xy[PDG_code] = new TH2D(("Silicon_Detector_InterStrip_xy_" + name).c_str(), ("Silicon_Detector_InterStrip_xy_" + name).c_str(), 8000, -40, 40, 1000, -0.05, 0.05);
    H_Silicon_Detector_InterStrip_xy[PDG_code]->GetXaxis()->SetTitle("X [mm]");
    H_Silicon_Detector_InterStrip_xy[PDG_code]->GetYaxis()->SetTitle("Y [mm]");
    H_Silicon_Detector_InterStrip_xy[PDG_code]->GetXaxis()->CenterTitle();
    H_Silicon_Detector_InterStrip_xy[PDG_code]->GetYaxis()->CenterTitle();

    H_Silicon_Detector_InterStrip_z[PDG_code] = new TH1D(("Silicon_Detector_InterStrip_z_" + name).c_str(), ("Silicon_Detector_InterStrip_z_" + name).c_str(), 4000, -0.2, 0.2);
    H_Silicon_Detector_InterStrip_z[PDG_code]->GetXaxis()->SetTitle("Z [um]");
    H_Silicon_Detector_InterStrip_z[PDG_code]->GetYaxis()->SetTitle("Counts");
    H_Silicon_Detector_InterStrip_z[PDG_code]->GetXaxis()->CenterTitle();
    H_Silicon_Detector_InterStrip_z[PDG_code]->GetYaxis()->CenterTitle();

    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det))
        {
            H_Silicon_Detector_Energy_Deposit_Det[PDG_code][det] = new TH1D(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_" + name).c_str(), ("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_" + name).c_str(), eSiliN_cal, 0, 10000);
            H_Silicon_Detector_Energy_Deposit_Det[PDG_code][det]->GetXaxis()->SetTitle("Energy [keV]");
            H_Silicon_Detector_Energy_Deposit_Det[PDG_code][det]->GetYaxis()->SetTitle("Counts");
            H_Silicon_Detector_Energy_Deposit_Det[PDG_code][det]->GetXaxis()->CenterTitle();
            H_Silicon_Detector_Energy_Deposit_Det[PDG_code][det]->GetYaxis()->CenterTitle();

            if (PDG_code == 0)
            {
                H_Silicon_Detector_Energy_Deposit_Det_without_interstrip[PDG_code][det] = new TH1D(("Silicon_Detector_Energy_Deposit_WithoutInterStrip_" + detectorName[det] + "_" + name).c_str(), ("Silicon_Detector_Energy_Deposit_WithoutInterStrip_" + detectorName[det] + "_" + name).c_str(), eSiliN_cal, 0, 10000);
                H_Silicon_Detector_Energy_Deposit_Det_without_interstrip[PDG_code][det]->GetXaxis()->SetTitle("Energy [keV]");
                H_Silicon_Detector_Energy_Deposit_Det_without_interstrip[PDG_code][det]->GetYaxis()->SetTitle("Counts");
                H_Silicon_Detector_Energy_Deposit_Det_without_interstrip[PDG_code][det]->GetXaxis()->CenterTitle();
                H_Silicon_Detector_Energy_Deposit_Det_without_interstrip[PDG_code][det]->GetYaxis()->CenterTitle();

                H_Silicon_Detector_Energy_Deposit_Det_Rear[PDG_code][det] = new TH2D(("Silicon_Detector_Energy_Deposit_Rear_" + detectorName[det] + "_" + name).c_str(), ("Silicon_Detector_Energy_Deposit_Rear_" + detectorName[det] + "_" + name).c_str(), 1000, 0, 10000, 1000, 0, 10000);
                H_Silicon_Detector_Energy_Deposit_Det_Rear[PDG_code][det]->GetXaxis()->SetTitle("Energy [keV]");
                H_Silicon_Detector_Energy_Deposit_Det_Rear[PDG_code][det]->GetYaxis()->SetTitle("Counts");
                H_Silicon_Detector_Energy_Deposit_Det_Rear[PDG_code][det]->GetXaxis()->CenterTitle();
                H_Silicon_Detector_Energy_Deposit_Det_Rear[PDG_code][det]->GetYaxis()->CenterTitle();
                
            }

            H_Silicon_Detector_Hit_Angle_Det[PDG_code][det] = new TH1D(("Silicon_Detector_Hit_Angle_" + detectorName[det] + "_" + name).c_str(), ("Silicon_Detector_Hit_Angle_" + detectorName[det] + "_" + name).c_str(), 900, 0, 90);
            H_Silicon_Detector_Hit_Angle_Det[PDG_code][det]->GetXaxis()->SetTitle("Angle [deg]");
            H_Silicon_Detector_Hit_Angle_Det[PDG_code][det]->GetYaxis()->SetTitle("Counts");
            H_Silicon_Detector_Hit_Angle_Det[PDG_code][det]->GetXaxis()->CenterTitle();
            H_Silicon_Detector_Hit_Angle_Det[PDG_code][det]->GetYaxis()->CenterTitle();

            H_Silicon_Detector_DL_Energy_Deposit_Det[PDG_code][det] = new TH1D(("Silicon_Detector_DL_Energy_Deposit_" + detectorName[det] + "_" + name).c_str(), ("Silicon_Detector_DL_Energy_Deposit_" + detectorName[det] + "_" + name).c_str(), eSiliN_cal, 0, 10000);
            H_Silicon_Detector_DL_Energy_Deposit_Det[PDG_code][det]->GetXaxis()->SetTitle("Energy [keV]");
            H_Silicon_Detector_DL_Energy_Deposit_Det[PDG_code][det]->GetYaxis()->SetTitle("Counts");
            H_Silicon_Detector_DL_Energy_Deposit_Det[PDG_code][det]->GetXaxis()->CenterTitle();
            H_Silicon_Detector_DL_Energy_Deposit_Det[PDG_code][det]->GetYaxis()->CenterTitle();

            H_Silicon_Detector_Energy_Deposit_SINGLE[PDG_code][det] = new TH1D(("Silicon_Detector_Energy_Deposit_SINGLE_" + detectorName[det] + "_" + name).c_str(), ("Silicon_Detector_Energy_Deposit_SINGLE_" + detectorName[det] + "_" + name).c_str(), eSiliN_cal, 0, 10000);
            H_Silicon_Detector_Energy_Deposit_SINGLE[PDG_code][det]->GetXaxis()->SetTitle("Energy [keV]");
            H_Silicon_Detector_Energy_Deposit_SINGLE[PDG_code][det]->GetYaxis()->SetTitle("Counts");
            H_Silicon_Detector_Energy_Deposit_SINGLE[PDG_code][det]->GetXaxis()->CenterTitle();
            H_Silicon_Detector_Energy_Deposit_SINGLE[PDG_code][det]->GetYaxis()->CenterTitle();

            H_Silicon_Detector_Energy_Deposit_NOCOINC[PDG_code][det] = new TH1D(("Silicon_Detector_Energy_Deposit_NOCOINC_" + detectorName[det] + "_" + name).c_str(), ("Silicon_Detector_Energy_Deposit_NOCOINC_" + detectorName[det] + "_" + name).c_str(), eSiliN_cal, 0, 10000);
            H_Silicon_Detector_Energy_Deposit_NOCOINC[PDG_code][det]->GetXaxis()->SetTitle("Energy [keV]");
            H_Silicon_Detector_Energy_Deposit_NOCOINC[PDG_code][det]->GetYaxis()->SetTitle("Counts");
            H_Silicon_Detector_Energy_Deposit_NOCOINC[PDG_code][det]->GetXaxis()->CenterTitle();
            H_Silicon_Detector_Energy_Deposit_NOCOINC[PDG_code][det]->GetYaxis()->CenterTitle();

            H_Silicon_Detector_Energy_Deposit_COINC[PDG_code][det] = new TH1D(("Silicon_Detector_Energy_Deposit_COINC_" + detectorName[det] + "_" + name).c_str(), ("Silicon_Detector_Energy_Deposit_COINC_" + detectorName[det] + "_" + name).c_str(), eSiliN_cal, 0, 10000);
            H_Silicon_Detector_Energy_Deposit_COINC[PDG_code][det]->GetXaxis()->SetTitle("Energy [keV]");
            H_Silicon_Detector_Energy_Deposit_COINC[PDG_code][det]->GetYaxis()->SetTitle("Counts");
            H_Silicon_Detector_Energy_Deposit_COINC[PDG_code][det]->GetXaxis()->CenterTitle();
            H_Silicon_Detector_Energy_Deposit_COINC[PDG_code][det]->GetYaxis()->CenterTitle();

            if (GetDetectorChannel(det) != 5)
            {
                int interstrip = GetDetector(det) * 1000 + ((GetDetectorChannel(det) * 2 + 1)) * 100 / 2;
                H_Silicon_Detector_InterStrip_Energy_Deposit_Det[PDG_code][interstrip] = new TH1D(("Silicon_Detector_InterStrip_Energy_Deposit_" + detectorName[det] + to_string(GetDetectorChannel(det) + 1) + "_" + name).c_str(), ("Silicon_Detector_InterStrip_Energy_Deposit_" + detectorName[det] + to_string(GetDetectorChannel(det) + 1) + "_" + name).c_str(), eSiliN_cal, 0, 10000);
                H_Silicon_Detector_InterStrip_Energy_Deposit_Det[PDG_code][interstrip]->GetXaxis()->SetTitle("Energy [keV]");
                H_Silicon_Detector_InterStrip_Energy_Deposit_Det[PDG_code][interstrip]->GetYaxis()->SetTitle("Counts");
                H_Silicon_Detector_InterStrip_Energy_Deposit_Det[PDG_code][interstrip]->GetXaxis()->CenterTitle();
                H_Silicon_Detector_InterStrip_Energy_Deposit_Det[PDG_code][interstrip]->GetYaxis()->CenterTitle();
            }
        }
    }
}

void WriteHistograms()
{

    Info("Writing Histograms");
    ANALYSIS_File->cd();
    for (auto &pair : H_E0)
    {
        int PDG = pair.first;
        if (H_E0[PDG]->GetEntries() == 0)
        {
            continue;
        }

        string name = SearchName(PDG);

        dir_Particle[PDG] = ANALYSIS_File->mkdir(name.c_str());
        dir_Initial[PDG] = dir_Particle[PDG]->mkdir("Initial");
        dir_Catcher[PDG] = dir_Particle[PDG]->mkdir("Catcher");
        dir_PlasticScintillator[PDG] = dir_Particle[PDG]->mkdir("PlasticScintillator");
        dir_Silicon_Detector[PDG] = dir_Particle[PDG]->mkdir("Silicon_Detector");
        dir_Silicon_Detector_InterStrip[PDG] = dir_Particle[PDG]->mkdir("Silicon_Detector_InterStrip");
        for (int det = 0; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorSiliBack(det))
            {
                dir_Silicon_Detector_Det[PDG][GetDetector(det)] = dir_Silicon_Detector[PDG]->mkdir(("D." + to_string(GetDetector(det))).c_str());
            }
        }

        dir_Initial[PDG]->cd();
        H_E0[PDG]->Write();
        H_x[PDG]->Write();
        H_y[PDG]->Write();
        H_z[PDG]->Write();
        H_px[PDG]->Write();
        H_py[PDG]->Write();
        H_pz[PDG]->Write();

        dir_Catcher[PDG]->cd();
        H_Catcher_Central_Energy_Deposit[PDG]->Write();
        H_Catcher_Side_Energy_Deposit[PDG]->Write();

        dir_PlasticScintillator[PDG]->cd();
        H_PlasticScintillator_Energy_Deposit[PDG]->Write();
        H_PlasticScintillator_Hit_Position_x[PDG]->Write();
        H_PlasticScintillator_Hit_Position_y[PDG]->Write();
        H_PlasticScintillator_Hit_Position_z[PDG]->Write();
        H_PlasticScintillator_Hit_Position_xy[PDG]->Write();
        H_PlasticScintillator_Hit_Angle[PDG]->Write();
        H_PlasticScintillator_Hit_Time[PDG]->Write();

        dir_Silicon_Detector[PDG]->cd();
        H_Silicon_Detector_Hit_Position_x[PDG]->Write();
        H_Silicon_Detector_Hit_Position_y[PDG]->Write();
        H_Silicon_Detector_Hit_Position_z[PDG]->Write();
        H_Silicon_Detector_Hit_Position_xy[PDG]->Write();
        H_Silicon_Detector_Hit_Position_xz[PDG]->Write();
        H_Silicon_Detector_Hit_Position_yz[PDG]->Write();
        H_Silicon_Detector_Hit_Angle[PDG]->Write();
        H_Silicon_Detector_Hit_Anglez[PDG]->Write();
        H_Silicon_Detector_Hit_Anglexy[PDG]->Draw("COLZ");
        H_Silicon_Detector_Hit_Anglexy[PDG]->Write();
        H_Silicon_Detector_Hit_Time[PDG]->Write();
        H_Silicon_Detector_Code[PDG]->Write();
        dir_Silicon_Detector_InterStrip[PDG]->cd();
        H_Silicon_Detector_InterStrip_xy[PDG]->Write();
        H_Silicon_Detector_InterStrip_z[PDG]->Write();

        for (int det = 0; det < SIGNAL_MAX; det++)
        {
            if (IsDetectorSiliStrip(det))
            {
                dir_Silicon_Detector_Det[PDG][GetDetector(det)]->cd();
                H_Silicon_Detector_Energy_Deposit_Det[PDG][det]->Write();
                H_Silicon_Detector_Hit_Angle_Det[PDG][det]->Write();

                if (GetDetectorChannel(det) != 5)
                {
                    dir_Silicon_Detector_InterStrip[PDG]->cd();
                    int interstrip = GetDetector(det) * 1000 + ((GetDetectorChannel(det) * 2 + 1)) * 100 / 2;
                    H_Silicon_Detector_InterStrip_Energy_Deposit_Det[PDG][interstrip]->Write();
                }
            }
        }
    }

    ANALYSIS_File->cd();
    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det))
        {
            H_Silicon_Detector_Energy_Deposit_Det[0][det]->Write();
            H_Silicon_Detector_Energy_Deposit_Det_Rear[0][det]->Write();

            TCanvas *c = new TCanvas(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All_WithWithoutInterStrip").c_str(), ("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All_WithWithoutInterStrip").c_str(), 800, 600);
            c->cd();
            H_Silicon_Detector_Energy_Deposit_Det[0][det]->SetLineColor(kBlack);
            H_Silicon_Detector_Energy_Deposit_Det[0][det]->Draw("HIST");
            H_Silicon_Detector_Energy_Deposit_Det_without_interstrip[0][det]->SetLineColor(kRed);
            H_Silicon_Detector_Energy_Deposit_Det_without_interstrip[0][det]->Draw("HIST SAME");
            c->Write();
            H_Silicon_Detector_Energy_Deposit_Det_without_interstrip[0][det]->Write();

            TCanvas *c2 = new TCanvas(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_pvsAll").c_str(), ("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_pvsAll").c_str(), 800, 600);
            c2->cd();
            TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
            int counter = 1;
            H_Silicon_Detector_Energy_Deposit_Det[0][det]->SetLineColor(counter);
            H_Silicon_Detector_Energy_Deposit_Det[0][det]->Draw("HIST");
            leg->AddEntry(H_Silicon_Detector_Energy_Deposit_Det[0][det], "All", "l");
            for (auto pair : H_E0)
            {
                if (pair.first == 0)
                    continue;

                string name = SearchName(pair.first);
                counter++;
                int PDG = pair.first;
                H_Silicon_Detector_Energy_Deposit_Det[PDG][det]->SetLineColor(counter);
                H_Silicon_Detector_Energy_Deposit_Det[PDG][det]->Draw("HIST SAME");
                leg->AddEntry(H_Silicon_Detector_Energy_Deposit_Det[PDG][det], name.c_str(), "l");
            }
            leg->Draw("SAME");
            c2->Write();

            TCanvas *c3 = new TCanvas(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_Coinc").c_str(), ("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_Coinc").c_str(), 800, 600);
            c3->cd();
            H_Silicon_Detector_Energy_Deposit_SINGLE[0][det]->SetLineColor(kBlack);
            H_Silicon_Detector_Energy_Deposit_SINGLE[0][det]->Draw("HIST");
            H_Silicon_Detector_Energy_Deposit_NOCOINC[0][det]->SetLineColor(kBlue);
            H_Silicon_Detector_Energy_Deposit_NOCOINC[0][det]->Draw("HIST SAME");
            H_Silicon_Detector_Energy_Deposit_COINC[0][det]->SetLineColor(kRed);
            H_Silicon_Detector_Energy_Deposit_COINC[0][det]->Draw("HIST SAME");

            H_Silicon_Detector_Energy_Deposit_SINGLE[0][det]->GetXaxis()->SetRangeUser(3200, 3400);
            H_Silicon_Detector_Energy_Deposit_COINC[0][det]->GetXaxis()->SetRangeUser(3200, 3400);
            
            cout << "Detector: " << detectorName[det] << endl;
            cout << "--- Ratio : " << H_Silicon_Detector_Energy_Deposit_COINC[0][det]->Integral() / H_Silicon_Detector_Energy_Deposit_SINGLE[0][det]->Integral() * 100 << endl;

            TLatex *text = new TLatex();
            text->SetNDC();
            text->SetTextSize(0.03);
            text->DrawLatex(0.1, 0.92, "Single");
            text->SetTextColor(kBlue);

            c3->Write();
            
        
        }
    }

    TH1D* H_StripH_Single[SIGNAL_MAX];
    TH1D* H_StripH_Coinc[SIGNAL_MAX];
    TH1D* H_StripH_NOCoinc[SIGNAL_MAX];

    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det))
        {
            if (GetDetector(det) == 1 || GetDetector(det) == 5)
            {
                if (GetDetectorChannel(det) == 1)
                {
                    H_StripH_Single[det] = (TH1D*)H_Silicon_Detector_Energy_Deposit_SINGLE[0][det]->Clone();
                    H_StripH_Coinc[det] = (TH1D*)H_Silicon_Detector_Energy_Deposit_COINC[0][det]->Clone();
                    H_StripH_NOCoinc[det] = (TH1D*)H_Silicon_Detector_Energy_Deposit_NOCOINC[0][det]->Clone();
                }
            }
            else
            {
                if (GetDetector(det) <= 4)
                {
                    H_StripH_Single[11]->Add(H_Silicon_Detector_Energy_Deposit_SINGLE[0][det]);
                    H_StripH_Coinc[11]->Add(H_Silicon_Detector_Energy_Deposit_COINC[0][det]);
                    H_StripH_NOCoinc[11]->Add(H_Silicon_Detector_Energy_Deposit_NOCOINC[0][det]);
                }
                if (GetDetector(det) >= 5)
                {
                    H_StripH_Single[51]->Add(H_Silicon_Detector_Energy_Deposit_COINC[0][det]);
                    H_StripH_Coinc[51]->Add(H_Silicon_Detector_Energy_Deposit_SINGLE[0][det]);
                    H_StripH_NOCoinc[51]->Add(H_Silicon_Detector_Energy_Deposit_NOCOINC[0][det]);
                }
            }
        }
    }

    for (int det = 0; det < SIGNAL_MAX; det++)
    {
        if (IsDetectorSiliStrip(det))
        {
            if (GetDetector(det) == 1 || GetDetector(det) == 5)
            {
                if (GetDetectorChannel(det) == 1)
                {
                    TCanvas *c = new TCanvas(("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_SUMMED").c_str(), ("Silicon_Detector_Energy_Deposit_" + detectorName[det] + "_All_WithWithoutInterStrip").c_str(), 800, 600);
                    H_StripH_Single[det]->Draw("HIST");
                    H_StripH_Coinc[det]->SetLineColor(kRed);
                    H_StripH_Coinc[det]->Write("HIST SAME");
                    H_StripH_NOCoinc[det]->SetLineColor(kBlue);
                    H_StripH_NOCoinc[det]->Write("HIST SAME");
                    c->Write();
                }
            }
        }
    }


    TCanvas *c3 = new TCanvas("PlasticScintillator_Energy_Deposit_All", "PlasticScintillator_Energy_Deposit_All", 800, 600);
    c3->cd();
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    int counter = 1;
    for (auto pair : H_E0)
    {
        if (pair.first == 0)
            continue;

        string name = SearchName(pair.first);
        counter++;
        int PDG = pair.first;
        // H_PlasticScintillator_Energy_Deposit[0]->Add(H_PlasticScintillator_Energy_Deposit[PDG], 1);
    }
    counter = 1 ;
    H_PlasticScintillator_Energy_Deposit[0]->SetLineColor(counter);
    H_PlasticScintillator_Energy_Deposit[0]->Draw("HIST");
    leg->AddEntry(H_PlasticScintillator_Energy_Deposit[0], "All", "l");
    for (auto pair : H_E0)
    {
        if (pair.first == 0)
            continue;

        string name = SearchName(pair.first);
        counter++;
        int PDG = pair.first;
        H_PlasticScintillator_Energy_Deposit[PDG]->SetLineColor(counter);
        H_PlasticScintillator_Energy_Deposit[PDG]->Draw("HIST SAME");
        leg->AddEntry(H_PlasticScintillator_Energy_Deposit[PDG], name.c_str(), "l");
    }
    leg->Draw("SAME");
    c3->Write();

    H->Write();

}