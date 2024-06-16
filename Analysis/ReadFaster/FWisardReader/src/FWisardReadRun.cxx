// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME FWisardReadRun
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "Signal.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_Signal(void *p = nullptr);
   static void *newArray_Signal(Long_t size, void *p);
   static void delete_Signal(void *p);
   static void deleteArray_Signal(void *p);
   static void destruct_Signal(void *p);
   static void streamer_Signal(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Signal*)
   {
      ::Signal *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Signal >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("Signal", ::Signal::Class_Version(), "Signal.h", 7,
                  typeid(::Signal), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Signal::Dictionary, isa_proxy, 16,
                  sizeof(::Signal) );
      instance.SetNew(&new_Signal);
      instance.SetNewArray(&newArray_Signal);
      instance.SetDelete(&delete_Signal);
      instance.SetDeleteArray(&deleteArray_Signal);
      instance.SetDestructor(&destruct_Signal);
      instance.SetStreamerFunc(&streamer_Signal);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Signal*)
   {
      return GenerateInitInstanceLocal((::Signal*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Signal*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Signal::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *Signal::Class_Name()
{
   return "Signal";
}

//______________________________________________________________________________
const char *Signal::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Signal*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int Signal::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Signal*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Signal::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Signal*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Signal::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Signal*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void Signal::Streamer(TBuffer &R__b)
{
   // Stream an object of class Signal.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      R__b >> Label;
      R__b >> Time;
      R__b >> Channel;
      R__b >> Multiplicity;
      R__b >> isValid;
      R__b.CheckByteCount(R__s, R__c, Signal::IsA());
   } else {
      R__c = R__b.WriteVersion(Signal::IsA(), kTRUE);
      R__b << Label;
      R__b << Time;
      R__b << Channel;
      R__b << Multiplicity;
      R__b << isValid;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Signal(void *p) {
      return  p ? new(p) ::Signal : new ::Signal;
   }
   static void *newArray_Signal(Long_t nElements, void *p) {
      return p ? new(p) ::Signal[nElements] : new ::Signal[nElements];
   }
   // Wrapper around operator delete
   static void delete_Signal(void *p) {
      delete ((::Signal*)p);
   }
   static void deleteArray_Signal(void *p) {
      delete [] ((::Signal*)p);
   }
   static void destruct_Signal(void *p) {
      typedef ::Signal current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_Signal(TBuffer &buf, void *obj) {
      ((::Signal*)obj)->::Signal::Streamer(buf);
   }
} // end of namespace ROOT for class ::Signal

namespace {
  void TriggerDictionaryInitialization_FWisardReadRun_Impl() {
    static const char* headers[] = {
"Signal.h",
nullptr
    };
    static const char* includePaths[] = {
"/home/local1/Documents/lib/FasterProcess2.0/include",
"/usr/local/root/6.26.10/include",
"/home/local1/Documents/lib/GTools1.0/include",
"/usr/local/XercesC/3.2.4/include",
"/home/local1/Documents/lib/fasterac-2.11/include",
"/softs/root/6.26.10/include/",
"/home/local1/Documents/2024_Analysis/Analysis/ReadFaster/FWisardReader/src/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "FWisardReadRun dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate(R"ATTRDUMP(Signal)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$Signal.h")))  Signal;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "FWisardReadRun dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "Signal.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"Signal", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("FWisardReadRun",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_FWisardReadRun_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_FWisardReadRun_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_FWisardReadRun() {
  TriggerDictionaryInitialization_FWisardReadRun_Impl();
}
