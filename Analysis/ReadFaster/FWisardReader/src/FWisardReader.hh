//======================================================================
/*! \file FWisardReader.hh
 *
 *  Include file for WISArD experiment analysis class
 */
//======================================================================

#ifndef FWISARD_READER_HH
#define FWISARD_READER_HH

#include "FReader.hh"
#include "FColor.hh"
//#include <vector>
//#include <deque>

//----------------------------------------------------------------------

#define FDATA_MAX   190      ///< Maximum coder label
#define SIGNAL_MAX  190      ///< Maximum number of signals
#define BETA_NUM     2      ///< Number of beta (SiPM) detector groups
#define SILI_NUM     8      ///< Number of silicon detector groups
#define BETA_SIZE    9      ///< Number of channels in beta (SiPM) detector groups
#define SILI_SIZE    6      ///< Number of channels in silicon detector groups

#define BETA_HI      1
#define BETA_LO      2

#define DEFAULT_HISTO_DIR   "../../../../../../../mnt/hgfs/shared-2/"

//======================================================================
/*! \class FWisardReader */
//======================================================================
class FWisardReader : public FReader
{
  protected:

  //------------------------------------------------------------
  /** @name Class variables */
  //@{
    //----------------------------------------------------------
    int     runNumber;

    string  detectorFileName;               ///< Detectors definition file name
    size_t  detectorNum;                    ///< Number of defined detectors channels
    string  detectorName [FDATA_MAX];
    int     detectorCoder[FDATA_MAX];      ///< Coder label of all signals
    int     detectorInfo [FDATA_MAX];      ///< Detector information of all signals (type identifier)

    int     detBeta[BETA_NUM][BETA_SIZE];   ///< Detector number of beta signals
    int     detSili[SILI_NUM][SILI_SIZE];   ///< Detector number of silicon signals

    int     coderDetector[FDATA_MAX];       ///< Detector number of all coders

    double detectorCalib [SIGNAL_MAX];
    bool Calibration = false;

    string  baseFileName;                   ///< Base name for histogram saving (from run file name)

    //////////////////
    //int SiStrip_table[SILI_NUM][SILI_SIZE];
    TFile *rf;
    //////////////////

    //----------------------------------------------------------
    //    histograms

    // groups analysis
    long double tGroup;                     ///< Time of group
    double      dtGroup[SIGNAL_MAX];
    int         mGroup;                     ///< Multiplicity
    int         mGroupBetaHi;               ///< Multiplicity
    int         mGroupBetaLo;               ///< Multiplicity
    int         mGroupSiStrip;              ///< Multiplicity
    int         mGroupSiBack;               ///< Multiplicity


    double  dtBetaHiMin;      ///< Coincidence window
    double  dtBetaHiMax;      ///< Coincidence window
    double  dtBetaLoMin;      ///< Coincidence window
    double  dtBetaLoMax;      ///< Coincidence window
    double  dtSiliMin;        ///< Coincidence window
    double  dtSiliMax;        ///< Coincidence window

    // energy
    int     nSili;
    double  eSiliMin;
    double  eSiliMax;
    int     nBetaHi;
    double  qBetaHiMin;
    double  qBetaHiMax;
    int     nBetaLo;
    double  qBetaLoMin;
    double  qBetaLoMax;


    /////////
    TTree * Tree;
    double Channel;
    double Energy;
    double Time;
    int Label;
    bool    treeSave;    ///< Some trees to save
    int event = 0;
    int multi;
      
    
    /////////

    int     nTrun;
    double  tRunMin;
    double  tRunMax;

    //----------------------------------------------------------
    // outputs information
    bool    optStat;      ///< Display statistics

    bool    histoSave;    ///< Some histograms to save
    string  histoDir;     ///< Histogram output files directory
  //@}

  public:

  //------------------------------------------------------------
  /** @name Constructor, destructor */
  //@{
             FWisardReader ( );
    virtual ~FWisardReader ( );
  //@}

  //----------------------------------------------------------
  /** @name Analysis setting functions */
  //@{
    static  void        PrintOptionInfo ( );
    virtual GStringList ProcessOptions ( bool warn = true );
  //@}

  //------------------------------------------------------------
  /** @name Data processing functions */
  //@{
    virtual int   InitDetectors  ( const string & fname );
    virtual int   InitHistograms ( );

    ////////////
    virtual int   InitTrees ( );
    virtual int WriteRunTime();
    ////////////

    virtual int   Init     ( );
    virtual int   End      ( );

      // redefinition of the base class function for group analysis
      // initialization
    virtual int   ProcessGroup   ( faster_data_p & data, u_int level = 0 );

      // user function (does nothing in the base class)
    virtual int   CoderProcessed ( u_int label, FCoder * coder, u_int level = 0 );

    virtual void  ReadConfigFile ( const string & file );
    virtual int   ConfigCommand ( const string & cmdline );
  //@}

  //------------------------------------------------------------
  /** @name Outputs related functions */
  //@{
    virtual int   SaveHisto ( const string & fname );
  //@}

  //----------------------------------------------------------
  /** @name user hooks from base classes */
  //@{
    virtual int   UserRunStart  ( );
    virtual int   UserFileStart ( );
    virtual int   UserFileStop  ( );
    virtual int   UserRunStop   ( );
  //@}

  //------------------------------------------------------------
  /** @name ROOT related function */
  //@{
    //  For ROOT encapsulation
    ClassDef(FWisardReader,0);
  //@}
};

//----------------------------------------------------------------------
//  Inline functions
#include "icc/FWisardReader.icc"


//======================================================================
#endif
