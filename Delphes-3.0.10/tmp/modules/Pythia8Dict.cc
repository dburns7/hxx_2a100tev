#define private public
#define protected public

/** \class
 *
 *  Lists classes to be included in cint dicitonary
 *
 *  $Date: 2013-07-04 13:19:15 +0200 (Thu, 04 Jul 2013) $
 *  $Revision: 1177 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/PileUpMergerPythia8.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PileUpMergerPythia8+;

#endif
//
// File generated by rootcint at Tue May  6 15:22:56 2014

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME tmpdImodulesdIPythia8Dict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "Pythia8Dict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void PileUpMergerPythia8_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_PileUpMergerPythia8(void *p = 0);
   static void *newArray_PileUpMergerPythia8(Long_t size, void *p);
   static void delete_PileUpMergerPythia8(void *p);
   static void deleteArray_PileUpMergerPythia8(void *p);
   static void destruct_PileUpMergerPythia8(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PileUpMergerPythia8*)
   {
      ::PileUpMergerPythia8 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PileUpMergerPythia8 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PileUpMergerPythia8", ::PileUpMergerPythia8::Class_Version(), "./modules/PileUpMergerPythia8.h", 27,
                  typeid(::PileUpMergerPythia8), DefineBehavior(ptr, ptr),
                  &::PileUpMergerPythia8::Dictionary, isa_proxy, 4,
                  sizeof(::PileUpMergerPythia8) );
      instance.SetNew(&new_PileUpMergerPythia8);
      instance.SetNewArray(&newArray_PileUpMergerPythia8);
      instance.SetDelete(&delete_PileUpMergerPythia8);
      instance.SetDeleteArray(&deleteArray_PileUpMergerPythia8);
      instance.SetDestructor(&destruct_PileUpMergerPythia8);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PileUpMergerPythia8*)
   {
      return GenerateInitInstanceLocal((::PileUpMergerPythia8*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::PileUpMergerPythia8*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *PileUpMergerPythia8::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *PileUpMergerPythia8::Class_Name()
{
   return "PileUpMergerPythia8";
}

//______________________________________________________________________________
const char *PileUpMergerPythia8::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PileUpMergerPythia8*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PileUpMergerPythia8::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PileUpMergerPythia8*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void PileUpMergerPythia8::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PileUpMergerPythia8*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *PileUpMergerPythia8::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PileUpMergerPythia8*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void PileUpMergerPythia8::Streamer(TBuffer &R__b)
{
   // Stream an object of class PileUpMergerPythia8.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PileUpMergerPythia8::Class(),this);
   } else {
      R__b.WriteClassBuffer(PileUpMergerPythia8::Class(),this);
   }
}

//______________________________________________________________________________
void PileUpMergerPythia8::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class PileUpMergerPythia8.
      TClass *R__cl = ::PileUpMergerPythia8::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fMeanPileUp", &fMeanPileUp);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fZVertexSpread", &fZVertexSpread);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fPTMin", &fPTMin);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*fPythia", &fPythia);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*fItInputArray", &fItInputArray);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*fInputArray", &fInputArray);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*fOutputArray", &fOutputArray);
      DelphesModule::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PileUpMergerPythia8(void *p) {
      return  p ? new(p) ::PileUpMergerPythia8 : new ::PileUpMergerPythia8;
   }
   static void *newArray_PileUpMergerPythia8(Long_t nElements, void *p) {
      return p ? new(p) ::PileUpMergerPythia8[nElements] : new ::PileUpMergerPythia8[nElements];
   }
   // Wrapper around operator delete
   static void delete_PileUpMergerPythia8(void *p) {
      delete ((::PileUpMergerPythia8*)p);
   }
   static void deleteArray_PileUpMergerPythia8(void *p) {
      delete [] ((::PileUpMergerPythia8*)p);
   }
   static void destruct_PileUpMergerPythia8(void *p) {
      typedef ::PileUpMergerPythia8 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PileUpMergerPythia8

/********************************************************
* tmp/modules/Pythia8Dict.cc
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtablePythia8Dict();

extern "C" void G__set_cpp_environmentPythia8Dict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__cpp_reset_tagtablePythia8Dict();
}
#include <new>
extern "C" int G__cpp_dllrevPythia8Dict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* PileUpMergerPythia8 */
static int G__Pythia8Dict_444_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   PileUpMergerPythia8* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new PileUpMergerPythia8[n];
     } else {
       p = new((void*) gvp) PileUpMergerPythia8[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new PileUpMergerPythia8;
     } else {
       p = new((void*) gvp) PileUpMergerPythia8;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__Pythia8DictLN_PileUpMergerPythia8));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Pythia8Dict_444_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) PileUpMergerPythia8::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Pythia8Dict_444_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) PileUpMergerPythia8::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Pythia8Dict_444_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) PileUpMergerPythia8::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Pythia8Dict_444_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      PileUpMergerPythia8::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Pythia8Dict_444_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((PileUpMergerPythia8*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Pythia8Dict_444_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) PileUpMergerPythia8::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Pythia8Dict_444_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) PileUpMergerPythia8::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Pythia8Dict_444_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) PileUpMergerPythia8::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__Pythia8Dict_444_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) PileUpMergerPythia8::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__Pythia8Dict_444_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   PileUpMergerPythia8* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new PileUpMergerPythia8(*(PileUpMergerPythia8*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__Pythia8DictLN_PileUpMergerPythia8));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef PileUpMergerPythia8 G__TPileUpMergerPythia8;
static int G__Pythia8Dict_444_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (PileUpMergerPythia8*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((PileUpMergerPythia8*) (soff+(sizeof(PileUpMergerPythia8)*i)))->~G__TPileUpMergerPythia8();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (PileUpMergerPythia8*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((PileUpMergerPythia8*) (soff))->~G__TPileUpMergerPythia8();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* PileUpMergerPythia8 */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncPythia8Dict {
 public:
  G__Sizep2memfuncPythia8Dict(): p(&G__Sizep2memfuncPythia8Dict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncPythia8Dict::*p)();
};

size_t G__get_sizep2memfuncPythia8Dict()
{
  G__Sizep2memfuncPythia8Dict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritancePythia8Dict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__Pythia8DictLN_PileUpMergerPythia8))) {
     PileUpMergerPythia8 *G__Lderived;
     G__Lderived=(PileUpMergerPythia8*)0x1000;
     {
       DelphesModule *G__Lpbase=(DelphesModule*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__Pythia8DictLN_PileUpMergerPythia8),G__get_linked_tagnum(&G__Pythia8DictLN_DelphesModule),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
     {
       ExRootTask *G__Lpbase=(ExRootTask*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__Pythia8DictLN_PileUpMergerPythia8),G__get_linked_tagnum(&G__Pythia8DictLN_ExRootTask),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
     {
       TTask *G__Lpbase=(TTask*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__Pythia8DictLN_PileUpMergerPythia8),G__get_linked_tagnum(&G__Pythia8DictLN_TTask),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
     {
       TNamed *G__Lpbase=(TNamed*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__Pythia8DictLN_PileUpMergerPythia8),G__get_linked_tagnum(&G__Pythia8DictLN_TNamed),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__Pythia8DictLN_PileUpMergerPythia8),G__get_linked_tagnum(&G__Pythia8DictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetablePythia8Dict() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__Pythia8DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__Pythia8DictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__Pythia8DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__Pythia8DictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__Pythia8DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__Pythia8DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__Pythia8DictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__Pythia8DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__Pythia8DictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__Pythia8DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<std::bidirectional_iterator_tag,TObject*,std::ptrdiff_t,const TObject**,const TObject*&>",117,G__get_linked_tagnum(&G__Pythia8DictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*>",117,G__get_linked_tagnum(&G__Pythia8DictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long>",117,G__get_linked_tagnum(&G__Pythia8DictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long,const TObject**>",117,G__get_linked_tagnum(&G__Pythia8DictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("map<TString,TString>",117,G__get_linked_tagnum(&G__Pythia8DictLN_maplETStringcOTStringcOlesslETStringgRcOallocatorlEpairlEconstsPTStringcOTStringgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("map<TString,TString,less<TString> >",117,G__get_linked_tagnum(&G__Pythia8DictLN_maplETStringcOTStringcOlesslETStringgRcOallocatorlEpairlEconstsPTStringcOTStringgRsPgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* PileUpMergerPythia8 */
static void G__setup_memvarPileUpMergerPythia8(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__Pythia8DictLN_PileUpMergerPythia8));
   { PileUpMergerPythia8 *p; p=(PileUpMergerPythia8*)0x1000; if (p) { }
   G__memvar_setup((void*)0,100,0,0,-1,G__defined_typename("Double_t"),-1,4,"fMeanPileUp=",0,(char*)NULL);
   G__memvar_setup((void*)0,100,0,0,-1,G__defined_typename("Double_t"),-1,4,"fZVertexSpread=",0,(char*)NULL);
   G__memvar_setup((void*)0,100,0,0,-1,G__defined_typename("Double_t"),-1,4,"fPTMin=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__Pythia8DictLN_Pythia8cLcLPythia),-1,-1,4,"fPythia=",0,"!");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__Pythia8DictLN_TIterator),-1,-1,4,"fItInputArray=",0,"!");
   G__memvar_setup((void*)0,85,0,1,G__get_linked_tagnum(&G__Pythia8DictLN_TObjArray),-1,-1,4,"fInputArray=",0,"!");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__Pythia8DictLN_TObjArray),-1,-1,4,"fOutputArray=",0,"!");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__Pythia8DictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarPythia8Dict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncPileUpMergerPythia8(void) {
   /* PileUpMergerPythia8 */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__Pythia8DictLN_PileUpMergerPythia8));
   G__memfunc_setup("PileUpMergerPythia8",1880,G__Pythia8Dict_444_0_1, 105, G__get_linked_tagnum(&G__Pythia8DictLN_PileUpMergerPythia8), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Init",404,(G__InterfaceMethod) NULL,121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Process",735,(G__InterfaceMethod) NULL,121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Finish",609,(G__InterfaceMethod) NULL,121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Class",502,G__Pythia8Dict_444_0_5, 85, G__get_linked_tagnum(&G__Pythia8DictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&PileUpMergerPythia8::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__Pythia8Dict_444_0_6, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&PileUpMergerPythia8::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__Pythia8Dict_444_0_7, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&PileUpMergerPythia8::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__Pythia8Dict_444_0_8, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&PileUpMergerPythia8::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__Pythia8DictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - insp", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__Pythia8Dict_444_0_12, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__Pythia8Dict_444_0_13, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&PileUpMergerPythia8::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__Pythia8Dict_444_0_14, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&PileUpMergerPythia8::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__Pythia8Dict_444_0_15, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&PileUpMergerPythia8::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__Pythia8Dict_444_0_16, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&PileUpMergerPythia8::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("PileUpMergerPythia8", 1880, G__Pythia8Dict_444_0_17, (int) ('i'), G__get_linked_tagnum(&G__Pythia8DictLN_PileUpMergerPythia8), -1, 0, 1, 1, 1, 0, "u 'PileUpMergerPythia8' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~PileUpMergerPythia8", 2006, G__Pythia8Dict_444_0_18, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncPythia8Dict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalPythia8Dict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {
}

static void G__cpp_setup_func13() {
}

static void G__cpp_setup_func14() {
}

static void G__cpp_setup_func15() {
}

static void G__cpp_setup_func16() {
}

static void G__cpp_setup_func17() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcPythia8Dict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
  G__cpp_setup_func13();
  G__cpp_setup_func14();
  G__cpp_setup_func15();
  G__cpp_setup_func16();
  G__cpp_setup_func17();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__Pythia8DictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__Pythia8DictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__Pythia8DictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__Pythia8DictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__Pythia8DictLN_TNamed = { "TNamed" , 99 , -1 };
G__linked_taginfo G__Pythia8DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__Pythia8DictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__Pythia8DictLN_TObjArray = { "TObjArray" , 99 , -1 };
G__linked_taginfo G__Pythia8DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__Pythia8DictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__Pythia8DictLN_TIterator = { "TIterator" , 99 , -1 };
G__linked_taginfo G__Pythia8DictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR = { "iterator<bidirectional_iterator_tag,TObject*,long,const TObject**,const TObject*&>" , 115 , -1 };
G__linked_taginfo G__Pythia8DictLN_TTask = { "TTask" , 99 , -1 };
G__linked_taginfo G__Pythia8DictLN_maplETStringcOTStringcOlesslETStringgRcOallocatorlEpairlEconstsPTStringcOTStringgRsPgRsPgR = { "map<TString,TString,less<TString>,allocator<pair<const TString,TString> > >" , 99 , -1 };
G__linked_taginfo G__Pythia8DictLN_ExRootTask = { "ExRootTask" , 99 , -1 };
G__linked_taginfo G__Pythia8DictLN_DelphesModule = { "DelphesModule" , 99 , -1 };
G__linked_taginfo G__Pythia8DictLN_Pythia8 = { "Pythia8" , 110 , -1 };
G__linked_taginfo G__Pythia8DictLN_Pythia8cLcLPythia = { "Pythia8::Pythia" , 99 , -1 };
G__linked_taginfo G__Pythia8DictLN_PileUpMergerPythia8 = { "PileUpMergerPythia8" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtablePythia8Dict() {
  G__Pythia8DictLN_TClass.tagnum = -1 ;
  G__Pythia8DictLN_TBuffer.tagnum = -1 ;
  G__Pythia8DictLN_TMemberInspector.tagnum = -1 ;
  G__Pythia8DictLN_TObject.tagnum = -1 ;
  G__Pythia8DictLN_TNamed.tagnum = -1 ;
  G__Pythia8DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__Pythia8DictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__Pythia8DictLN_TObjArray.tagnum = -1 ;
  G__Pythia8DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__Pythia8DictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__Pythia8DictLN_TIterator.tagnum = -1 ;
  G__Pythia8DictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR.tagnum = -1 ;
  G__Pythia8DictLN_TTask.tagnum = -1 ;
  G__Pythia8DictLN_maplETStringcOTStringcOlesslETStringgRcOallocatorlEpairlEconstsPTStringcOTStringgRsPgRsPgR.tagnum = -1 ;
  G__Pythia8DictLN_ExRootTask.tagnum = -1 ;
  G__Pythia8DictLN_DelphesModule.tagnum = -1 ;
  G__Pythia8DictLN_Pythia8.tagnum = -1 ;
  G__Pythia8DictLN_Pythia8cLcLPythia.tagnum = -1 ;
  G__Pythia8DictLN_PileUpMergerPythia8.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtablePythia8Dict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_TClass);
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_TObject);
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_TNamed);
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_TObjArray);
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_TIterator);
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR);
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_TTask);
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_maplETStringcOTStringcOlesslETStringgRcOallocatorlEpairlEconstsPTStringcOTStringgRsPgRsPgR);
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_ExRootTask);
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_DelphesModule);
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_Pythia8);
   G__get_linked_tagnum_fwd(&G__Pythia8DictLN_Pythia8cLcLPythia);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__Pythia8DictLN_PileUpMergerPythia8),sizeof(PileUpMergerPythia8),-1,324864,(char*)NULL,G__setup_memvarPileUpMergerPythia8,G__setup_memfuncPileUpMergerPythia8);
}
extern "C" void G__cpp_setupPythia8Dict(void) {
  G__check_setup_version(30051515,"G__cpp_setupPythia8Dict()");
  G__set_cpp_environmentPythia8Dict();
  G__cpp_setup_tagtablePythia8Dict();

  G__cpp_setup_inheritancePythia8Dict();

  G__cpp_setup_typetablePythia8Dict();

  G__cpp_setup_memvarPythia8Dict();

  G__cpp_setup_memfuncPythia8Dict();
  G__cpp_setup_globalPythia8Dict();
  G__cpp_setup_funcPythia8Dict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncPythia8Dict();
  return;
}
class G__cpp_setup_initPythia8Dict {
  public:
    G__cpp_setup_initPythia8Dict() { G__add_setup_func("Pythia8Dict",(G__incsetup)(&G__cpp_setupPythia8Dict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initPythia8Dict() { G__remove_setup_func("Pythia8Dict"); }
};
G__cpp_setup_initPythia8Dict G__cpp_setup_initializerPythia8Dict;

