//======================================================================
/*! \file FWisardReader.cpp
 *
 *  Source file for WISArD experiment analysis class
 */
//======================================================================

#include "FWisardReader.hh"

//======================================================================
//  CONSTRUCTOR /DESTRUCTOR
//======================================================================

/*! NuBall / FASTER data reader constructor.
 */
FWisardReader::FWisardReader() : FReader()
{
  // signals definitions
  for (size_t i = 0; i < FDATA_MAX; ++i)
    coderDetector[i] = -1;
  for (size_t i = 0; i < SIGNAL_MAX; ++i)
  {
    detectorCoder[i] = -1;
    detectorInfo[i] = -1;
  }

  for (size_t i = 0; i < BETA_NUM; ++i)
  {
    for (size_t j = 0; j < BETA_SIZE; ++j)
      detBeta[i][j] = -1;
  }
  for (size_t i = 0; i < SILI_NUM; ++i)
  {
    for (size_t j = 0; j < SILI_SIZE; ++j)
      detSili[i][j] = -1;
  }

  detectorNum = 0;
  detectorFileName = "Data/detectors.dat";

  //------------------------------------------------------------
  // outputs
  optStat = false;

  histoSave = false;
  histoDir = DEFAULT_HISTO_DIR;
}

//----------------------------------------------------------------------

/*! NuBall / FASTER data reader destructor.
 */
FWisardReader::~FWisardReader()
{
}

//======================================================================
//
//======================================================================

/*! Print the list of options for the class on the terminal.*/
void FWisardReader::PrintOptionInfo()
{
  FReader::PrintOptionInfo();
  printf("* FWisardReader analysis class options\n");
  printf("    -c<file>       set the configuration file (not used)\n");
  printf("    -s             print analysis statistics for each run\n");
}

//----------------------------------------------------------------------
/*! Function processing the command-line options.
 *  The valid options are removed from the list for processing in derived
 *  classes.
 *  The function returns the list of unprocessed options.
 *
 *  If redefined in derived class, it should be built as follows in the
 *  derived class:
 *  - check options in derived class;
 *  - remove recognised options from the list;
 *  - if there are remaining options, call the parent class ProcessOptions
 *    function with the remaining list.
 *
 *  \param  warn    warn on unrecognized options
 */
GStringList FWisardReader::ProcessOptions(bool warn)
{
  GStringList newOptions;

  // process base class options
  FReader::ProcessOptions(false);

  // process remaining options
  for (auto &p : cmdOptions)
  {
    bool ok = false; // identified option
    // string  argstr = cmdOptions.At(i);
    string argstr = *p;

    //----------------------------------------------------------
    // configuration file
    if (argstr.substr(0, 2) == "-c")
    {
      // *** not yet implemented ***
      ok = true;
      string cfg = GString(argstr.substr(2)).NoEndSpace(" \n\t\r\"'");
      GLogMessage("<WIS> Configuration file: " + cfg, 1);
    }

    // configuration file
    if (argstr.substr(0, 2) == "-s")
    {
      ok = true;
      optStat = true;
    }

    //----------------------------------------------------------
    if (!ok) // set in remaining options
    {
      newOptions.push_back(GString(argstr));

      if (warn)
        GLogWarning("Undefined argument " + argstr);
    }
  }

  return (newOptions);
}

//======================================================================
//  DATA PROCESSING FUNCTIONS
//======================================================================

//----------------------------------------------------------------------
/*! User function initializing the analysis
 */
int FWisardReader::Init()
{
  WisMessage("<WIS> Initialization function");
  int error = 0;

  // read the detectors information
  error = InitDetectors(detectorFileName);

  // create histograms
  if (error == 0)
    error = InitHistograms();

  if (error == 0)
    error =0;

  // initializes the coders
  FReader::Init();

  cout << endl;

  return (error);
}

//----------------------------------------------------------------------
/*! Function reading the detectors definition to set the coder channels,
 *  detector numbers, detector names...
 *  \param  fname   definition file name
 */
int FWisardReader::InitDetectors(const string &fname)
{
  int error = 0;

  FILE *fp = fopen(fname.c_str(), "r");
  if (fp != NULL)
  {
    detectorNum = 0;
    while (!feof(fp))
    {
      char fline[512];
      fgets(fline, 511, fp);

      if (!feof(fp))
      {
        string defline = GString(fline).NoEndSpace();
        if (defline.length() > 0)
        {
          if (defline.at(0) != '#') // not a comment
          {
            istringstream iss(defline);
            int coder;
            string name;
            iss >> coder >> name;

            if (!iss.fail())
            {
              detectorCoder[detectorNum] = coder;
              detectorName[detectorNum] = name;

              // set detector groups
              if (name.substr(0, 4) == "Beta")
              {
                detectorInfo[detectorNum] = 0;
                int d = -1;
                sscanf(&name.c_str()[6], "%d", &d);
                if ((d > 0) && (d <= 9))
                {
                  if (name.substr(4, 2) == "Hi")
                  {
                    detBeta[0][d - 1] = detectorNum;
                    detectorInfo[detectorNum] += 10 + d;
                    // cerr << "** BetaHi " << d << " = " << detectorInfo[detectorNum] << endl;
                  }
                  else if (name.substr(4, 2) == "Lo")
                  {
                    detBeta[1][d - 1] = detectorNum;
                    detectorInfo[detectorNum] += 20 + d;
                    // cerr << "** BetaLo " << d << " = " << detectorInfo[detectorNum] << endl;
                  }
                  else
                    GLogWarning("Bad beta detector: " + defline);
                }
                else
                  GLogWarning("Bad beta detector number: " + defline);
              }
              else if (name.substr(0, 1) == "D")
              {
                detectorInfo[detectorNum] = 100;
                // name.at(2) = ' ';
                if (name.at(3) == 'R')
                  name.at(3) = '6';
                int s = -1;
                int d = -1;
                sscanf(&name.c_str()[1], "%d.%d", &s, &d);
                if ((s > 0) && (s <= 8) && (d > 0) && (d <= 6))
                {
                  detSili[s - 1][d - 1] = detectorNum;
                  detectorInfo[detectorNum] += 10 * s + d;
                }
                else
                {
                  GLogWarning("Bad silicon detector number: " + defline);
                }
              }
              else
                GLogWarning("Error in definition line: " + defline);

              detectorNum++;
            }
            else
              GLogWarning("Error in definition line: " + defline);
          }
        }
      }
    } // read loop

    fclose(fp);
    GLogMessage("<WIS> Number of defined detectors: " + GGetString((int)detectorNum));

    // update coder -> detector conversion table
    for (size_t i = 0; i < detectorNum; ++i)
    {
      // cerr << detectorName[i] << " -> " << GGetString(detectorInfo[i],3) << endl;
      if (detectorCoder[i] >= 0)
        
        {
          coderDetector[detectorCoder[i]] = i;

        }
    }
  }
  else
    GLogWarning("Error opening definition file: " + fname);

  return (error);
}

int FWisardReader::InitTrees()
{
  int error = 0;

  treeSave = true;

  rf = new TFile((histoDir + "/" + baseFileName + ".root").c_str(), "RECREATE");

  Tree = new TTree("Tree", "Tree");
  Tree->Branch("Label", &Label, "Label/I");
  Tree->Branch("Channel", &Channel, "Channel/D");
  Tree->Branch("Time", &Time, "Time/D");
  Tree->Branch("Energy", &Energy, "Energy/D");

  return (error);
}

int FWisardReader::WriteRunTime()
{
  int error = 0;

  string start_time;
  string stop_time;
  string line;

  ifstream file(fasterFile->GetFileName().substr(0, fasterFile->GetFileName().length() - 10) + ".setup");
  if (file.is_open())
  {
    GLogMessage("<SAM> Setup file found");
    while (std::getline(file, line))
    {
      if (line.find("- Start date") != std::string::npos)
      {
        start_time = line.substr(17, 50);
      }
      if (line.find("- Stop date") != std::string::npos)
      {
        stop_time = line.substr(16, 50);
      }
      if (line.find("- Stop  date") != std::string::npos)
      {
        stop_time = line.substr(16, 50);
      }
    }
    file.close();

    rf->cd();
    TNamed("Start_time", start_time.c_str()).Write();
    TNamed("Stop_time", stop_time.c_str()).Write();
    cout << "Start time: " << start_time << endl;
    cout << "Stop time: " << stop_time << endl;
  }

  else
  {
    GLogWarning("Error opening .setup file: " + fasterFile->GetFileName().substr(0, fasterFile->GetFileName().length() - 10) + ".setup");
  }

  return (error);
}
//----------------------------------------------------------------------
/*! Function defining histograms
 */
int FWisardReader::InitHistograms()
{
  int error = 0;
  return (error);
}
//----------------------------------------------------------------------
/*! Function called at end of analysis.*/
int FWisardReader::End()
{
  int error = 0;
  GLogMessage("");
  GLogMessage("<WIS> ---- Analysis End ----", 2);

  return (error);
}

//======================================================================

//----------------------------------------------------------------------

/*! Function that analyzes the group data currently read from acquisition.
 *  It is redefined from FReader base class.
 *  \param  data    data from acquisition
 *  \param  level   group depth
 */
int FWisardReader::ProcessGroup(faster_data_p &data, u_int level)
{
  tGroup = faster_data_hr_clock_ns(data);
  mGroup = 0;
  mGroupBetaHi = 0;
  mGroupBetaLo = 0;
  mGroupSiStrip = 0;
  mGroupSiBack = 0;

  event++;
  multi = 0;

  memset(dtGroup, 0, SIGNAL_MAX * sizeof(double));

  // call the base class function
  int error = FReader::ProcessGroup(data, level);

  return (error);
}

//----------------------------------------------------------------------

/*! Function overloaded from base class to fill histograms, if defined.
 *  \param  label   coder acquisition identifier
 *  \param  coder   coder with signal
 *  \param  level   group depth level
 */
int FWisardReader::CoderProcessed(u_int label, FCoder *coder, u_int level)
{
  int error = 0;

  if (coder != NULL)
  {
    if (level == 0)
    {
      // get the last data for the coder (in case of multiplicity > 1)
      FCoderData *cdata = coder->GetData(coder->GetMultiplicity() - 1);
      int det = label;
      long double time_ns = cdata->GetTime();
      long double trun_ns = time_ns - runT0;
      double trun_s = trun_ns * 0.000000001L;
      double value = cdata->GetMeasure(); // raw data
      
      Channel = value + ((double) rand() / (RAND_MAX)) - 0.5;
      Time = trun_ns;
      Label = det;
      Tree->Fill();
    }
    else
    {
      error = 3;
      GLogWarning("FWisardReader::CoderProcessed(...): Bad group level " + GGetString(level));
    }
  }
  else
  {
    error = 4;
    GLogWarning("FWisardReader::CoderProcessed(...): NULL coder data pointer");
  }

  return (error);
}

//======================================================================
//
//======================================================================

void FWisardReader::ReadConfigFile(const string &file)
{
  FILE *fp = fopen(file.c_str(), "r");
  if (fp != NULL)
  {
    while (!feof(fp))
    {
      char fline[512];
      fgets(fline, 511, fp);

      if (!feof(fp))
      {
        string cmdline = GString(fline).NoEndSpace();
        if (cmdline.length() > 0)
        {
          if (cmdline.at(0) != '#')
          {
            ConfigCommand(cmdline);
          }
        }
      }
    } // read loop

    fclose(fp);
  }
  else
    GLogWarning("Error opening configuration file: " + file);
}

//----------------------------------------------------------------------
int FWisardReader::ConfigCommand(const string &cmdline)
{
  int error = 0;

  GString line(cmdline);

  std::cout<< line << std::endl;

  size_t i0, i1;
  string cmd = (line.GetWord(i0, i1, 0, " ,;:.!?'[]{}()=+-*&%$#\n\t\r")).ToUpper();
  string args = line.substr(i1 + 1);
  istringstream isargs(args);

  bool ok = false;
  bool argErr = false;

  // verbosity
  if (cmd == "VERB")
  {
    ok = true;
    int verb = 0;
    isargs >> verb;
    if (!isargs.fail())
      GSetVerboseLevel(verb);
    else
      argErr = true;
  }

  // event information
  if (cmd == "INFORATE")
  {
    ok = true;
    int rate = 0;
    isargs >> rate;
    if (!isargs.fail())
      SetInfoRate(rate);
    else
      argErr = true;
  }

  // output histogram directory
  if (cmd == "OUTPUT_DIR")
  {
    ok = true;
    histoDir = GString(args).GetWord(0, " ,;:!?'[]{}()=*&%$#\n\t\r");
  }

  //------------------------------------------------------------
  if (argErr)
  {
    GLogWarning("Bad arguments for command " + cmd);
  }
  //------------------------------------------------------------
  else if (!ok)
  {
    GLogWarning("Undefined configuration command " + cmd + " ignored");
  }

  return (error);
}

//======================================================================
//
//======================================================================

//======================================================================
/*! Save histograms in a ROOT file.
 *  \param  fname   ROOT file name
 */
int FWisardReader::SaveHisto(const string &fname)
{
  int error = 0;

  if (rf->IsOpen())
  {
    GLogMessage("<WIS> Output histograms file: " + fname);

    Tree->Write();
    //------------------------------------------------------------
    rf->Close();
  }
  else
  {
    error = 1;
    GLogWarning("Error creating histograms file: " + fname);
  }
  delete rf;

  return (error);
}

//======================================================================
//  USER HOOKS FROM BASE CLASSES
//======================================================================

//----------------------------------------------------------------------
/*! Function called at start of a new run.
 *  It clears the monitor histograms
 */
int FWisardReader::UserRunStart()
{
  WisMessage("<WIS> ---- Analysis Run Start ----", 1);

  int error = 0;

  string fileName = fasterFile->GetFileName();

  // in case of a set of run files, extract the run number
  runNumber = 9999;
  FRunSerie *fastRunSerie = dynamic_cast<FRunSerie *>(fasterFile);
  if (fastRunSerie != NULL)
  {
    runNumber = fastRunSerie->GetCurrentRun();
    if (runNumber == 0)
    {
      size_t p = fileName.find("Run_");
      if (p != string::npos)
      {
        string runStr = fileName.substr(p + 4, 4);
        runNumber = GString(runStr).ReadInt();
      }
    }
    GLogMessage("<WIS> Current run: " + GGetString(runNumber));
  }
  else
    GLogWarning("Undefined run number for run files group - set to " + GGetString(runNumber));

  // extract the base file name (for hosto saving)
  string baseName = gSystem->BaseName(fileName.c_str());
  size_t p = baseName.find(".fast");
  if (p != string::npos)
    baseName = baseName.substr(0, p - 5);

  baseFileName = baseName;
  GLogMessage("<WIS> Base file name: " + baseFileName + "\n");
  InitTrees();
  WriteRunTime();


  // cerr << "t(FILE *) = " << typeid(FILE *).name() << "  S=" << sizeof(FILE) << endl;
  // cerr << "t(gzFile) = " << typeid(gzFile).name() << "  S=" << sizeof(gzFile_s) << endl;
  return (error);
}

//----------------------------------------------------------------------
/*! Function called at start of a new single file.*/
int FWisardReader::UserFileStart()
{
  WisMessage("<WIS> ---- Analysis File Start ----", 4);

  return (0);
}

/*! Function called at end of a single file.*/
int FWisardReader::UserFileStop()
{
  WisMessage("<WIS> ---- Analysis File Stop ----", 4);

  return (0);
}

//----------------------------------------------------------------------
/*! Function called at end of run processing.*/
int FWisardReader::UserRunStop()
{
  WisMessage("<WIS> ---- Analysis Run Stop ----", 1);

  int error = 0;

  if (optStat)
  {
    GLogMessage("");
    GLogMessage("<WIS> Statistics for run " + GGetString(runNumber));
    GLogMessage("<F>    coder           data        scalers    mult.max      missed");
    for (size_t ic = 0; ic <= maxCoders; ++ic)
    {
      FCoder *coder = GetCoder(ic);
      if (coder != NULL)
      {
        ULong64_t nc = coder->GetDataCounts();
        ULong64_t ns = coder->GetScalerCounts();
        u_int mm = coder->GetMaxMultiplicity();
        size_t md = coder->GetScalerRunMismatch();
        if ((nc > 0) || (ns > 0))
          GLogMessage("<F>    " + GGetString(ic, 5) + " " + GGetString(nc, 14) + " " + GGetString(ns, 14) + "    " + GGetString(mm, 8) + "    " + GGetString(md, 8));
      }
    }

    GLogMessage("<F>  grp level   counts");
    for (u_int i = 0; i < FGROUP_DEPTH_MAX; ++i)
    {
      ULong64_t ng = GetGroupCounts(i);
      if (ng > 0)
        GLogMessage("<F>  " + GGetString(i, 3) + " " + GGetString(ng, 14));
    }
    /*
       GLogMessage ( "<WIS> - Groups" );
       for ( u_int i = 0; i < FGROUP_DEPTH_MAX; ++i )
       {
         if (groupCount[i] > 0)
           GLogMessage ( "<WIS>        level " + GGetString(i) + " : "
                    + GGetString(groupCount[i],12) );
       }

       GLogMessage ( "<WIS> - Signals" );
       for ( u_int i = 0; i < detectorNum; ++i )
       {
         if ( detectorCoder[i] > 0 )
         {
           int     coder = detectorCoder[i];
           string  name  = detectorName[i];
           while (name.length() < 16) name = " " + name;

           GLogMessage ( "<WIS> " + name + " :  "
                       + GGetString(coderCount[coder],12) + " data  /"
                       + GGetString(scalerCount[coder],8) + " scaler" );
         }
       }
   */
  }

  // save histograms
  
  
    gSystem->mkdir(histoDir.c_str(), kTRUE);
    string hfile = histoDir + "/" + baseFileName + ".root";

    SaveHisto(hfile);
  

  return (error);
}


