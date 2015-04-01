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

/**
   *This class is used to remplace the ON and OFF Data by a poissonian distribution in each bands (zenith, offset, area and energy).
   *The mean of the ON in each Band is obtained by calculating the mean of the signal and the mean of the background in the ON.
   *The mean of the signal is calculated with the functionExpectedExcess of the computeResult Class by usinf the spectral parameters of the law generating the MCs given as an input
   *the mean of the background in each band is calculated by using the formula: alpha_iband * N_OFF for E<1Tev and a power law of spectral index 2.7 for energy >1 TeV.
   *The OFF data are not changed for E<1 TeV. For E>1TeV, the off data are calculated by drawing a poissonian law with a mean value for each energy bin determined by the integration of the power law with a spectral index of 2.7 between each boundaries of the bin.
   */


//Constructor
START::CheckFit_MC::CheckFit_MC(BandsFactory &BandsFact ,std::vector<Band> &BandArray, Hypothesis &hypo):trandom(0)
{
  fBand_realdata=BandArray;
  //Initialize the N_On of the array to O
  double n_on=0.;
  BandsFact.Change_N_ON(BandArray,n_on);
  

  //An object hypothesis that contains the spectral parameters of the law given as an input to generate the MCs
  fHypothesis = &hypo;

  //An object Comput Result in order to use the method FunctionExpectedExcess
  fCompRes = new ComputeResults(BandArray,*fHypothesis);
}

//Destructor
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
   //From a energy of 1 TeV, the backgroud fallow a power law of spectral index 2.7:phi_O E^(-2.7)
   int index=2.7;
   int index_int=index-1;
     
   i_max=BandArray[0].ebin.size()-1;
   
   //Determine the index of the energy tab from which the energy is equal to 1 TeV in order to calculate the normalisation of the power law fallow by the background from 1 TeV.
   while(BandArray[0].ebin[i_min].GetEmin()<=1){
    i_min=i_min + 1;
  }
   int_spectre_background=(pow(BandArray[0].ebin[i_max].GetEmax(),-index_int)-pow(BandArray[0].ebin[i_min].GetEmin(),-index_int))/(-index_int);
   //Calcul of the spectrum normalisation for the background (phi_0) from the total events between 1 TeV and Emax of the energy bin for each band in the real data
   for(unsigned int iband(0); iband<BandArray.size(); iband++) {   
     sum_background=0;
     for(unsigned int ibin(i_min); ibin<BandArray[iband].ebin.size(); ibin++) {
       sum_background=sum_background+fBand_realdata[iband].ebin[ibin].GetOff();
     }
     
     tab_phiO_ON[iband]=fBand_realdata[iband].GetAlphaRun()*sum_background/(int_spectre_background);
     tab_phiO_Off[iband]=sum_background/(int_spectre_background);
   }
   
   for(unsigned int iband(0); iband<BandArray.size(); iband++) {
     if(BandArray[iband].GetKeepBand()==0) continue; // We are interested only in the selected bands
     for(unsigned int ibin(0); ibin<BandArray[iband].ebin.size(); ibin++) {
       if (BandArray[iband].ebin[ibin].GetKeepBin()==0) continue;// skip energy bins below threshold
       
       //mean signal: calculated with the FunctionExpectedExcess from ComputeResult
       mean_signal=fCompRes->FunctionExpectedExcess(&BandArray[iband],&BandArray[iband].ebin[ibin],fHypothesis->GetFittedParameters(), -1, -1);
       
       // the background of the ON is calculated by using the formula alpha_iband *N_OFF for energy<1 TeV
       if(ibin<i_min){
	 mean_background=BandArray[iband].GetAlphaRun()*BandArray[iband].ebin[ibin].GetOff();
       }
       //for energy> 1TeV, we calculate the number of OFF events expected by the power law by integrating the law between emin and emax of each energy bin and then we draw a poissonian law over the previous mean value.The OFF data are fill by this poissonian value for E>1TeV.
       //for energy> 1TeV, we calculate the number of background events expected in the ON data by the power law by integrating the law between emin and emax of each energy bin. For the ON data, the normalisation of the power law fallow by the background take into account the value of alpha for each band 
       else{
	 OFF=trandom.Poisson(tab_phiO_Off[iband]*((pow(BandArray[0].ebin[ibin].GetEmax(),-index_int)-pow(BandArray[0].ebin[ibin].GetEmin(),-index_int))/(-index_int)));
	 BandArray[iband].ebin[ibin].SetOff(OFF);
	 mean_background=tab_phiO_ON[iband]*((pow(BandArray[0].ebin[ibin].GetEmax(),-index_int)-pow(BandArray[0].ebin[ibin].GetEmin(),-index_int))/(-index_int));
       }
       
       On_mean= mean_signal + mean_background;
       //the number of events in the ON is calculated with a poisonnian distribution.(Poisson(On_mean)= Poisson(mean_signal) + Poisson(mean_background)
       On = trandom.Poisson(On_mean);
       //Fill de On of the object mydata by the previous value of ON
       BandArray[iband].ebin[ibin].SetOn(On);
       
     }
     
   }
 }
