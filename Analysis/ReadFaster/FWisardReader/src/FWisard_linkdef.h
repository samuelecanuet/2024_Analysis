//======================================================================

#ifndef FWISARD_LINKDEF_H
#define FWISARD_LINKDEF_H

//======================================================================

#ifdef __MAKECINT__

//----------------------------------------------------------------------
// functions from WISArD analysis layer

//#pragma link C++ namespace FWisard;

#pragma link C++ class FWisardReader;

//----------------------------------------------------------------------
// functions from FisColor.hh

#pragma link C++ global wisColTxt;
#pragma link C++ global wisColEnd;
#pragma link C++ function WisMessage;

//----------------------------------------------------------------------

#endif

//======================================================================
#endif
