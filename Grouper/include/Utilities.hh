#ifndef UTILITIES_HH
#define UTILITIES_HH

#include "Messenger.hh"
#include "Detectors.hh"


// Open ROOT File
TFile* MyTFile(string filename, string option)
{
    TFile *file = new TFile(filename.c_str(), option.c_str());
    filename = filename.substr(filename.find_last_of("/\\") + 1);
    if (file->IsZombie())
    {
        Error("Not Found: " + filename);
    }
    else if (option == "READ")
    {
        Success("Read: " + filename);
    }
    else if (option == "RECREATE")
    {
        Success("Recreated: " + filename);
    }
    else if (option == "UPDATE")
    {
        Success("Updated: " + filename);
    }
    else if (option == "NEW")
    {
        Success("Created: " + filename);
    }
    else
    {
        Error("Unknown option: " + filename);
    }
    return file;
}

// Open ROOT Tree
TTree* MyTTree(TFile* file, string tree_name)
{
    TTree *tree = (TTree *)file->Get(tree_name.c_str());
    if (tree == NULL)
    {
        Error("Not Found: " + tree_name + " in " + file->GetName());   
    }
    return tree;
}

#endif

