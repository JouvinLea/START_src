// STL
#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <cstdlib>

//root
#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "Math/QuantFuncMathCore.h"
#include "TUUID.h"

// START
#include "ComputeResults.hh"
#include "Band.hh"
#include "Hypothesis.hh"
#include "BandsFactory.hh"
#include "CheckFit_MC.hh"
// Utilities
#define DEBUG 0
#include "debugging.hh"


//Constructor
START::CheckFit_MC::CheckFit_MC(BandsFactory &BandsFact ,std::vector<Band> &BandArray, Hypothesis &hypo):trandom()
{
  //Initialize the N_On of the array to O
  double n_on=0.;
  BandsFact.Change_N_ON(BandArray,n_on);
  

  //An object hypothesis that contains the fittedparameter of the fit obtained with Starfit
  fHypothesis = &hypo;

  //An object Comput Result in order to use the method FunctionExpectedExcess
  fCompRes = new ComputeResults(BandArray,*fHypothesis);
}

START::CheckFit_MC::~CheckFit_MC(){
  if (fCompRes!=0) delete fCompRes;
  fCompRes = 0;
}

void START::CheckFit_MC::N_On_MC(std::vector<Band> &BandArray)
 { 
   double mean_signal, mean_background;
   double On_mean, On;
   for(unsigned int iband(0); iband<BandArray.size(); iband++) {
    
     for(unsigned int ibin(0); ibin<BandArray[iband].ebin.size(); ibin++) {
       
       //mean signal: calculated with the FunctionExpectedExcess from ComputeResult
       mean_signal=fCompRes->FunctionExpectedExcess(&BandArray[iband],&BandArray[iband].ebin[ibin],fHypothesis->GetFittedParameters(), -1, -1);
       // the background of the ON is calculated by using the formula alpha_iband *N_OFF
       mean_background=BandArray[iband].GetAlphaRun()*BandArray[iband].ebin[ibin].GetOff();
       On_mean= mean_signal + mean_background;
       //the number of events in the ON is calculated with a poisonnian distribution.(Poisson(On_mean)= Poisson(mean_signal) + Poisson(mean_background)
       On = trandom.Poisson(On_mean);
       //Fill de On of the object mydata by the previous value of ON
       BandArray[iband].ebin[ibin].SetOn(On);
       
     }
     
   }
 }
