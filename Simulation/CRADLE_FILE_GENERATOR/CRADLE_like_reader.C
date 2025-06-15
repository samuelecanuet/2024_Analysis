#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

using namespace std;

int CRADLE_like_reader()
{

    int Z = 17;
    int A = 33;

    std::ostringstream filename;
    filename << "33Ar";
    filename << "/z" << Z << ".a" << A << "_ENSDF";
    std::ifstream radDataFile((filename.str()).c_str());

    string line;
    double excitationEnergy = 0.;
    double lifetime;

    int counter=0;
    while (getline(radDataFile, line))
    {
      // cout<<line<<endl;

      if (!line.compare(0, 1, "#"))
      {
        // Comment line
        continue;
      }
      else if (!line.compare(0, 1, "P"))
      {
        // Parent line
        std::istringstream iss(line);
        string p;
        string flag;

        iss >> p >> excitationEnergy >> flag >> lifetime;
        counter++;
        continue;
      }
      // cout << "Lifetime: " << lifetime << endl;
      string mode;
      double daughterExcitationEnergy = 0;
      double intensity = 0;
      double Q = 0;
      string modifier;
      string flag;

      std::istringstream iss(line);
      iss >> mode >> daughterExcitationEnergy >> flag >> intensity >> Q >> modifier;     
    
    }
    cout << counter << endl;

    return 0;

}