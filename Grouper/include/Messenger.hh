#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream> 
#include <dirent.h>
#include <gsl/gsl_statistics.h>
#include "TFile.h"
#include "../../../lib/SignalDict/Signal.h"

#include "/home/local1/Documents/lib/GTools1.0/include/GString.hh"

using namespace std;

string RED = "\033[1;31m";
string GREEN = "\033[1;32m";
string YELLOW = "\033[1;33m";
string BLUE = "\033[1;34m";
string MAGENTA = "\033[1;35m";
string CYAN = "\033[1;36m";
string WHITE = "\033[1;37m";
string RESET = "\033[0m";

void Error(const string& message) {
    cout << RED << "<ERROR> " << message << RESET << endl;
    exit(0);
}

void Error(const char *message) {
    cout << RED << "<ERROR> " << message << RESET << endl;
    exit(0);
}

void Warning(const string& message) {
    cout << YELLOW << "<WARNING> " << message << RESET << endl;
}

void Warning(const char *message) {
    cout << YELLOW << "<WARNING> " << message << RESET << endl;
}

void Info(const string& message) {
    cout << BLUE << " <INFO>  " << message << RESET << endl;
}

void Info(const char *message) {
    cout << BLUE << " <INFO>  " << message << RESET << endl;
}

void Success(const string& message) {
    cout << GREEN << "<SUCCESS> " << message << RESET << endl;
}

void Success(const char *message) {
    cout << GREEN << "<SUCCESS> " << message << RESET << endl;
}

void Verbose(const string& message, int verbose, int this_verbose) {
    if (verbose >= this_verbose) cout << MAGENTA << "<VERBOSE " << verbose << ">" << message << RESET << endl;
}

void Verbose(Signal signal, int verbose, int this_verbose) {
    if (verbose >= this_verbose) cout << MAGENTA << "<VERBOSE " << verbose << ">" << signal << RESET << endl;
}

void Start(const string& message) {
    cout << CYAN << "<START> " << message << RESET << endl;
}

void Start(const char *message) {
    cout << CYAN << "<START> " << message << RESET << endl;
}

void ProgressBar(int cEntry, int TotalEntries, clock_t start, clock_t Current, string Prefix = "")
{
  if (cEntry % 10000 == 0 && cEntry > 2 * 1000)
  {
    Current = clock();
    const Char_t *Color;
    Double_t Frac = 1.0 * cEntry / TotalEntries;
    Double_t Timeclock = ((double)(Current - start) / CLOCKS_PER_SEC);
    Double_t TimeLeft = Timeclock * (1 / Frac - 1.);
    Color = "\e[1;31m";

    cout << Form(("\r"+Prefix+" Entries : ").c_str())
         << TotalEntries
         << " --- "
         << Form("%4.2f", 100. * cEntry / TotalEntries) << " %"
         << " --- "
         << " Time Left : " << Form("%2d min ", (int)TimeLeft / 60)
         << Form("%02d sec", (int)TimeLeft % 60)
         << flush;
  }

  if (cEntry == TotalEntries-1)
  {
    Current = clock();
    const Char_t *Color;
    Double_t Frac = 1.0 * cEntry / TotalEntries;
    Double_t Timeclock = ((double)(Current - start) / CLOCKS_PER_SEC);
    Double_t TimeLeft = Timeclock * (1 / Frac - 1.);
    Color = "\e[1;31m";
    cout << Form(("\r"+Prefix+" Entries : ").c_str())
         << TotalEntries
         << " --- "
         << Form("%4.2f", 100. * cEntry / TotalEntries) << " %"
         << " --- "
         << " Time Left : " << Form("%2d min ", (int)TimeLeft / 60)
         << Form("%02d sec", (int)TimeLeft % 60)
         << flush;
    cout << endl;
  }
}

void ProgressCounter(int cEntry, int TotalEntries, string Prefix = "")
{
    cout << BLUE
         << Form(("\r"+Prefix+" : ").c_str())
         << cEntry
         << " / "
         << TotalEntries
         << flush;

  if (cEntry == TotalEntries-1)
  {
    cout << BLUE
         << Form(("\r"+Prefix+" : ").c_str())
         << "Completed "
         << RESET
         << endl;
  }
}

string formatValueWithError(double value, double error, string type = "latex") {
    // Calculate the number of decimal places in the error
    int errorDecimals = error > 0 ? -floor(log10(error)) : 0;

    // Calculate the multiplier for rounding the error up
    double multiplier = pow(10, errorDecimals);

    // Round the error up to the next higher number with the same number of significant digits
    double roundedError = ceil(error * multiplier) / multiplier;

    // Calculate the number of decimal places in the rounded error
    int valueDecimals = roundedError > 1.0 ? 0 : errorDecimals;

    // Format the value and rounded error
    stringstream ssValue, ssError;
    ssValue.precision(valueDecimals);
    ssValue << fixed << value;
    ssError.precision(errorDecimals);
    ssError << fixed << roundedError;

    stringstream ss;
    if (type == "latex") {
        ss << ssValue.str() << " #pm " << ssError.str();
    } else {
        ss << ssValue.str() << " +/- " << ssError.str();
    }

    return ss.str();
}

