#ifndef MESSENGER_HH
#define MESSENGER_HH

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream> 
#include <sstream>
#include <iomanip>
#include <dirent.h>
#include <gsl/gsl_statistics.h>
#include "TFile.h"
#include "../../../lib/SignalDict/Signal.h"

using namespace std;

string RED = "\033[1;31m";
string GREEN = "\033[1;32m";
string YELLOW = "\033[1;33m";
string BLUE = "\033[1;34m";
string MAGENTA = "\033[1;35m";
string CYAN = "\033[1;36m";
string WHITE = "\033[1;37m";
string RESET = "\033[0m";

int GENERAL_INDENT = 12;

void Error(const string& message) {
    cout << RED << left << setw(GENERAL_INDENT) << " <ERROR>" << message << RESET << endl;
    exit(0);
}

void Error(const char *message) {
    cout << RED << left << setw(GENERAL_INDENT) << " <ERROR>" << message << RESET << endl;
    exit(0);
}

void Warning(const string& message, int indent=0) {
    cout << YELLOW << left << setw(GENERAL_INDENT) << " <WARNING>";
    if (indent > 0)
    {
        for (int i = 0; i < indent; i++) cout << " ";
        cout << "├─";
    }
    cout << message << RESET << endl;

}

void Warning(const char *message, int indent = 0) {
    cout << YELLOW << left << setw(GENERAL_INDENT) << " <WARNING>";
    if (indent > 0)
    {
        for (int i = 0; i < indent; i++) cout << " ";
        cout << "├─";
    }
    cout << message << RESET << endl;
}

void Info(const string& message, int indent=0) {
    cout << BLUE << left << setw(GENERAL_INDENT) << " <INFO>  ";
    if (indent > 0)
    {
        for (int i = 0; i < indent; i++) cout << " ";
        cout << "├─";
    }
    cout << message << RESET << endl;
}

void Info(const char *message, int indent=0) {
    Info(string(message), indent);
}

void Success(const string& message, int indent=0) {
    cout << GREEN << left << setw(GENERAL_INDENT) << " <SUCCESS>";
    if (indent > 0)
    {
        for (int i = 0; i < indent; i++) cout << " ";
        cout << "├─";
    }
    cout << message << RESET << endl;
}

void Success(const char *message, int indent=0) {
    Success(string(message), indent);
}

void Verbose(const string& message, int verbose, int this_verbose) {
    if (verbose >= this_verbose) cout << MAGENTA << "<VERBOSE " << verbose << ">" << message << RESET << endl;
}

void Verbose(Signal signal, int verbose, int this_verbose) {
    if (verbose >= this_verbose) cout << MAGENTA << "<VERBOSE " << verbose << ">" << signal << RESET << endl;
}

void Start(const string& message, int indent=0) {
    cout << CYAN << left << setw(GENERAL_INDENT) << " <START>  ";
    if (indent > 0)
    {
        for (int i = 0; i < indent; i++) cout << " ";
        cout << "├─";
    }
    cout << message << RESET << endl;
}

void Start(const char *message, int indent = 0) {
    Start(string(message), indent);
}

void ProgressBar(int cEntry, int TotalEntries, clock_t start, clock_t Current, string Prefix = "", int Step=10000)
{
  if (cEntry % Step == 0 && cEntry > 2 * 1000)
  {
    Current = clock();
    const Char_t *Color;
    Double_t Frac = 1.0 * cEntry / TotalEntries;
    Double_t Timeclock = ((double)(Current - start) / CLOCKS_PER_SEC);
    Double_t TimeLeft = Timeclock * (1. / Frac - 1.);
    Color = "\e[1;31m";

    cout << ("\r"+Prefix+" Entries : ").c_str()
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
    Double_t TimeLeft = Timeclock * (1. / Frac - 1.);
    Color = "\e[1;31m";
    cout << ("\r"+Prefix+" Entries : ").c_str()
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

void ProgressCounter(int cEntry, int TotalEntries, string Prefix = "", int stepping = 1, int indent = 0)
{
    if (cEntry % stepping == 0 || cEntry == TotalEntries - 1) // Refresh only at stepping intervals or the last entry
    {
        cout << BLUE << setw(GENERAL_INDENT); 
        if (indent > 0)
        {
            for (int i = 0; i < indent; i++)
                cout << " ";
            cout << "├─";
        }
        cout << Form(("\r" + Prefix + " : ").c_str())
        << cEntry
        << " / "
        << TotalEntries
        << flush;
    }

    if (cEntry == TotalEntries - 1)
    {
        cout << BLUE
             << Form(("\r" + Prefix + " : ").c_str())
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
        ss << ssValue.str() << " ± " << ssError.str();
    }

    return ss.str();
}

#endif

