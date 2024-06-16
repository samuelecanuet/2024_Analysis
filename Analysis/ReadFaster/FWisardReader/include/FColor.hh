//======================================================================
/*! \file FColor.hh
 *
 *  Include file for terminal colors definitions
 */
//======================================================================

#ifndef FIS_COLOR_HH
#define FIS_COLOR_HH

#include "FCommon.hh"
//----------------------------------------------------------------------

extern string wisColTxt;
extern string wisColEnd;

inline void WisMessage ( const string & text, u_int verb=0 )
  { GLogMessage ( wisColTxt + text + wisColEnd, verb ); }

//======================================================================
#endif
