
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
#include "check_fit.hh"
// Utilities
#define DEBUG 0
#include "debugging.hh"


//Constructor
START::CheckFit::CheckFit(BandsFactory &BandsFact ,std::vector<Band> &BandArray, Hypothesis &hypo):trandom()
{
  //Initialize the N_On of the array to O
  double n_on=0.;
  BandsFact.Change_N_ON(BandArray,n_on);
  //fBandArray = BandArray;

  //An object hypothesis that contains the parameter of the fit
  fHypothesis = &hypo;

  //An object Comput Result in order use the method FunctionExpectedExcess
  fCompRes = new ComputeResults(BandArray,*fHypothesis);
}

START::CheckFit::~CheckFit(){
  if (fCompRes!=0) delete fCompRes;
  fCompRes = 0;
}

void START::CheckFit::add_N_On(std::vector<Band> &BandArray)
 { 
   double mean_signal, mean_background;
   double On_mean, On;
   //std:: cout << fHypothesis->GetFittedParameters()[0] <<std::endl;
   //std:: cout << fHypothesis->GetFittedParameters()[1] <<std::endl;
   std:: cout << "Je suis dans add N ON" <<std::endl;
   for(unsigned int iband(0); iband<BandArray.size(); iband++) {
    
     for(unsigned int ibin(0); ibin<BandArray[iband].ebin.size(); ibin++) {
       // functionexpectedexcess attend des pointeurs pour les bandes et l energie, il faut donc lui passer fbandArray par reference
       //std:: cout << "phi0 =" << fHypothesis->GetFittedParameters()[0] << std::endl;
       //std:: cout << "gamma =" <<fHypothesis->GetFittedParameters()[1] << std::endl;
       
       
       mean_signal=fCompRes->FunctionExpectedExcess(&BandArray[iband],&BandArray[iband].ebin[ibin],fHypothesis->GetFittedParameters(), -1, -1);
     
       mean_background=BandArray[iband].GetAlphaRun()*BandArray[iband].ebin[ibin].GetOff();
       std:: cout << "jai passe mean background" <<std::endl;
       ;
       On_mean= mean_signal + mean_background;
       
       On = trandom.Poisson(On_mean);
       
       BandArray[iband].ebin[ibin].SetOn(On);
       
     }
     
   }
   
 }
//Ensuite faire un fit avec ces donnees. J ai pas compris si la fonction de minimizefactory ou on pouvait donner les donnees du fit a l avance nous sert la. Je pense que c etait peut etre que dans la class hypothesis on peut creeer un tableau qui nous sert vec les donnees gamma et phi 0 qui nous sert avant pour calculer l expected excess. 
