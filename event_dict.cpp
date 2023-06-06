// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME event_dict
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
#include "Event.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_NewEvent(void *p = nullptr);
   static void *newArray_NewEvent(Long_t size, void *p);
   static void delete_NewEvent(void *p);
   static void deleteArray_NewEvent(void *p);
   static void destruct_NewEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NewEvent*)
   {
      ::NewEvent *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::NewEvent >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("NewEvent", ::NewEvent::Class_Version(), "Event.h", 15,
                  typeid(::NewEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::NewEvent::Dictionary, isa_proxy, 4,
                  sizeof(::NewEvent) );
      instance.SetNew(&new_NewEvent);
      instance.SetNewArray(&newArray_NewEvent);
      instance.SetDelete(&delete_NewEvent);
      instance.SetDeleteArray(&deleteArray_NewEvent);
      instance.SetDestructor(&destruct_NewEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NewEvent*)
   {
      return GenerateInitInstanceLocal((::NewEvent*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::NewEvent*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr NewEvent::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *NewEvent::Class_Name()
{
   return "NewEvent";
}

//______________________________________________________________________________
const char *NewEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NewEvent*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int NewEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NewEvent*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *NewEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NewEvent*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *NewEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NewEvent*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void NewEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class NewEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(NewEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(NewEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_NewEvent(void *p) {
      return  p ? new(p) ::NewEvent : new ::NewEvent;
   }
   static void *newArray_NewEvent(Long_t nElements, void *p) {
      return p ? new(p) ::NewEvent[nElements] : new ::NewEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_NewEvent(void *p) {
      delete ((::NewEvent*)p);
   }
   static void deleteArray_NewEvent(void *p) {
      delete [] ((::NewEvent*)p);
   }
   static void destruct_NewEvent(void *p) {
      typedef ::NewEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::NewEvent

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = nullptr);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 216,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));

      ::ROOT::AddClassAlternate("vector<int>","std::vector<int, std::allocator<int> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)nullptr)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEfloatgR_Dictionary();
   static void vectorlEfloatgR_TClassManip(TClass*);
   static void *new_vectorlEfloatgR(void *p = nullptr);
   static void *newArray_vectorlEfloatgR(Long_t size, void *p);
   static void delete_vectorlEfloatgR(void *p);
   static void deleteArray_vectorlEfloatgR(void *p);
   static void destruct_vectorlEfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<float>*)
   {
      vector<float> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<float>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<float>", -2, "vector", 216,
                  typeid(vector<float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));

      ::ROOT::AddClassAlternate("vector<float>","std::vector<float, std::allocator<float> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<float>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<float>*)nullptr)->GetClass();
      vectorlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float> : new vector<float>;
   }
   static void *newArray_vectorlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float>[nElements] : new vector<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEfloatgR(void *p) {
      delete ((vector<float>*)p);
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] ((vector<float>*)p);
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<float>

namespace ROOT {
   static TClass *vectorlENewEventgR_Dictionary();
   static void vectorlENewEventgR_TClassManip(TClass*);
   static void *new_vectorlENewEventgR(void *p = nullptr);
   static void *newArray_vectorlENewEventgR(Long_t size, void *p);
   static void delete_vectorlENewEventgR(void *p);
   static void deleteArray_vectorlENewEventgR(void *p);
   static void destruct_vectorlENewEventgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<NewEvent>*)
   {
      vector<NewEvent> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<NewEvent>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<NewEvent>", -2, "vector", 216,
                  typeid(vector<NewEvent>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlENewEventgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<NewEvent>) );
      instance.SetNew(&new_vectorlENewEventgR);
      instance.SetNewArray(&newArray_vectorlENewEventgR);
      instance.SetDelete(&delete_vectorlENewEventgR);
      instance.SetDeleteArray(&deleteArray_vectorlENewEventgR);
      instance.SetDestructor(&destruct_vectorlENewEventgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<NewEvent> >()));

      ::ROOT::AddClassAlternate("vector<NewEvent>","std::vector<NewEvent, std::allocator<NewEvent> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<NewEvent>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlENewEventgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<NewEvent>*)nullptr)->GetClass();
      vectorlENewEventgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlENewEventgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlENewEventgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<NewEvent> : new vector<NewEvent>;
   }
   static void *newArray_vectorlENewEventgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<NewEvent>[nElements] : new vector<NewEvent>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlENewEventgR(void *p) {
      delete ((vector<NewEvent>*)p);
   }
   static void deleteArray_vectorlENewEventgR(void *p) {
      delete [] ((vector<NewEvent>*)p);
   }
   static void destruct_vectorlENewEventgR(void *p) {
      typedef vector<NewEvent> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<NewEvent>

namespace {
  void TriggerDictionaryInitialization_event_dict_Impl() {
    static const char* headers[] = {
"Event.h",
nullptr
    };
    static const char* includePaths[] = {
"/home/noah/root/include/",
"/home/noah/mdt/reconstruction/mdtreco/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "event_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$Event.h")))  NewEvent;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "event_dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "Event.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"NewEvent", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("event_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_event_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_event_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_event_dict() {
  TriggerDictionaryInitialization_event_dict_Impl();
}
