#include "Grouper/include/Detectors.hh"
#include "Grouper/include/Utilities.hh"

int syscalib()
{
    FLAG2025 = true;
    InitDetectors("Grouper/Config_Files/sample.pid");

    TFile *Exp = MyTFile((DIR_ROOT_DATA_CALIBRATED + "Calibrated_" + to_string(YEAR) + "_new.root").c_str(), "READ");
    
    for (int det = 1; det < SIGNAL_MAX; det ++)
    {
        if (IsDetectorSiliStrip(det))
        {
            TF1 *f = (TF1*)Exp->Get(("Calibration_" + detectorName[det]).c_str());

            cout << "Detector " << detectorName[det] << ":" << endl;
            cout << f->GetParameter(0) << "+/-" << f->GetParError(0) << endl;
            cout << f->GetParameter(1) << "+/-" << f->GetParError(1) << endl;
            cout << f->GetParameter(2) << "+/-" << f->GetParError(2) << endl;

            TF1 *bij = InvertFunction(f);
            double CH = bij->Eval(3350) - bij->Eval(3354);

            cout << " 4 keV in channels: " << CH << endl;
            cout << "error on 4 keV : " << CH * f->GetParError(1)/sqrt(40) << endl;
        }
    }


    return 0;
}