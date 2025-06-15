#include "Detectors.hh"

int strip_ref = 1;
int detector_ref = 1;

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

        if (x1 == strip_ref)
            continue;

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

    return chi2;
}

double Chi2(const TMultiGraph* mg1, const TMultiGraph* mg2) {
    if (!mg1 || !mg2) {
        Error("One of the graphs is null.");
        return -1;
    }

    // transformting TMultiGraph to TGraphErrors point y point

    int counter1 = 0;
    TGraphErrors *g1 = new TGraphErrors();
    for (int i = 0; i < mg1->GetListOfGraphs()->GetSize(); ++i)
    {
        TGraphErrors *temp = dynamic_cast<TGraphErrors *>(mg1->GetListOfGraphs()->At(i));
        for (int j = 0; j < temp->GetN(); ++j)
        {
            g1->AddPoint(temp->GetX()[j], temp->GetY()[j]);
            g1->SetPointError(g1->GetN() - 1, temp->GetErrorX(j), temp->GetErrorY(j));

            counter1++;
        }
    }

    int counter2 = 0;
    TGraphErrors *g2 = new TGraphErrors();
    for (int i = 0; i < mg2->GetListOfGraphs()->GetSize(); ++i)
    {
        TGraphErrors *temp = dynamic_cast<TGraphErrors *>(mg2->GetListOfGraphs()->At(i));
        for (int j = 0; j < temp->GetN(); ++j)
        {
            g2->AddPoint(temp->GetX()[j], temp->GetY()[j]);
            g2->SetPointError(g2->GetN() - 1, temp->GetErrorX(j), temp->GetErrorY(j));
            counter2++;
        }
    }

    cout << "Counter1: " << counter1 << ", Counter2: " << counter2 << endl;

    return Chi2(g1, g2);
}


struct Position
{
    double x;
    double y;
    double z;
    double theta;

    bool operator<(const Position& other) const {
        if (x*x + y*y + z*z <= other.x*other.x + other.y*other.y + other.z*other.z)
            return true;
        else
            return false;
    }

    bool operator==(const Position& other) const {
        return (x == other.x && y == other.y && z == other.z && theta == other.theta);
    }

    bool operator!=(const Position& other) const {
        return !(*this == other);
    }

    friend ostream& operator<<(ostream& os, const Position& pos) {
        os << "x: " << pos.x << ", y: " << pos.y << ", z: " << pos.z << ", theta: " << pos.theta;
        return os;
    }
};

bool Accepting_Position(Position pos)
{
    // if (pos.y > 0.0)
    //     return false;
    // if (abs(pos.y - pos.z) > 1.)
    //     return false;
    // if (pos.y != 0.5 || pos.z != 0.0)
    //     return false;
    // if (pos.theta != 0.0)
    //     return false;
    return true;
}



