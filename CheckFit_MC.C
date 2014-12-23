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
START::CheckFit_MC::CheckFit_MC(BandsFactory &BandsFact ,std::vector<Band> &BandArray, Hypothesis &hypo):trandom(0)
{
  fBand_realdata=BandArray;
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
   double OFF;
   double sum_background,int_spectre_background ;
   int i_max, i_min;
   double tab_phiO_ON[BandArray.size()];
   double tab_phiO_Off[BandArray.size()];
   i_max=BandArray[0].ebin.size()-1;
   
   while(BandArray[0].ebin[i_min].GetEmin()<=1){
    i_min=i_min + 1;
  }
   std::cout << i_min << std::endl;
   int_spectre_background=(pow(BandArray[0].ebin[i_max].GetEmax(),-1.7)-pow(BandArray[0].ebin[i_min].GetEmin(),-1.7))/(-1.7);
   
   for(unsigned int iband(0); iband<BandArray.size(); iband++) {   
     sum_background=0;
     for(unsigned int ibin(i_min); ibin<BandArray[iband].ebin.size(); ibin++) {
       sum_background=sum_background+fBand_realdata[iband].ebin[ibin].GetOff();
     }
     
     tab_phiO_ON[iband]=fBand_realdata[iband].GetAlphaRun()*sum_background/(int_spectre_background);
     tab_phiO_Off[iband]=sum_background/(int_spectre_background);
   }
   
   for(unsigned int iband(0); iband<BandArray.size(); iband++) {
     //std::cout << BandArray.size() << std::endl;
     if(BandArray[iband].GetKeepBand()==0) continue; // We are interested only in the selected bands
     for(unsigned int ibin(0); ibin<BandArray[iband].ebin.size(); ibin++) {
       if (BandArray[iband].ebin[ibin].GetKeepBin()==0) continue;// skip energy bins below threshold
       //mean signal: calculated with the FunctionExpectedExcess from ComputeResult
       //std::cout << BandArray[iband].ebin.size() << std::endl;
       
       mean_signal=fCompRes->FunctionExpectedExcess(&BandArray[iband],&BandArray[iband].ebin[ibin],fHypothesis->GetFittedParameters(), -1, -1);
       
       // the background of the ON is calculated by using the formula alpha_iband *N_OFF
       if(ibin<i_min){
	 //std::cout <<ibin << "ibin<i_min" << std::endl;
	 mean_background=BandArray[iband].GetAlphaRun()*BandArray[iband].ebin[ibin].GetOff();
       }
       else{
	 //std::cout <<ibin << "ibin>i_min" << std::endl; 
	 OFF=trandom.Poisson(tab_phiO_Off[iband]*((pow(BandArray[0].ebin[ibin].GetEmax(),-1.7)-pow(BandArray[0].ebin[ibin].GetEmin(),-1.7))/(-1.7)));
	 BandArray[iband].ebin[ibin].SetOff(OFF);
	 mean_background=tab_phiO_ON[iband]*((pow(BandArray[0].ebin[ibin].GetEmax(),-1.7)-pow(BandArray[0].ebin[ibin].GetEmin(),-1.7))/(-1.7));
       }
       
       On_mean= mean_signal + mean_background;
       //the number of events in the ON is calculated with a poisonnian distribution.(Poisson(On_mean)= Poisson(mean_signal) + Poisson(mean_background)
       On = trandom.Poisson(On_mean);
       //Fill de On of the object mydata by the previous value of ON
       BandArray[iband].ebin[ibin].SetOn(On);
       
     }
     
   }
 }
