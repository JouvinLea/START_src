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
START::CheckFit::CheckFit(BandsFactory &BandsFact ,std::vector<Band> &BandArray, int N_iteration, Hypothesis &hypo)
{
  //Initialize the N_On of the array to O
  double n_on=0.;
  BandsFact.Change_N_ON(BandArray,n_on);
  fBandArray = BandArray;
  // Define the number of point we want for the phi_O and the gamma reconstructed distribution
  f_iteration= N_iteration;

  //An object hypothesis that contains the parameter of the fit
  fHypothesis = &hypo;

  //An object Comput Result in order use the method FunctionExpectedExcess
  fCompRes = new ComputeResults(BandArray,*fHypothesis);
}

START::CheckFit::~CheckFit()
{
}

void START::CheckFit::add_N_On()
 { 
   double mean_signal, mean_background;
   double On_mean, On;

   for(unsigned int iband(0); iband<fBandArray.size(); iband++) {
    
     for(unsigned int ibin(0); ibin<fBandArray[iband].ebin.size(); ibin++) {
       // functionexpectedexcess attend des pointeurs pour les bandes et l energie, il faut donc lui passer fbandArray par reference
       //std:: cout << "phi0 =" << fHypothesis->GetFittedParameters()[0] << std::endl;
       //std:: cout << "gamma =" <<fHypothesis->GetFittedParameters()[1] << std::endl;
       mean_signal=fCompRes->FunctionExpectedExcess(&fBandArray[iband],&fBandArray[iband].ebin[ibin],fHypothesis->GetFittedParameters(), -1, -1);
       //std:: cout << "mean_signal = " << mean_signal << std::endl ;
       mean_background=fBandArray[iband].GetAlphaRun()*fBandArray[iband].ebin[ibin].GetOff();
       //std:: cout << "background OFF = " << fBandArray[iband].ebin[ibin].GetOff() << std::endl ;
       //std:: cout << "mean_background = " << mean_background << std::endl ;
       On_mean= mean_signal + mean_background;
       TRandom trandom;
       On = trandom.Poisson(On_mean);
       //std:: cout << "On = " << On << std::endl;
       fBandArray[iband].ebin[ibin].SetOn(On);
       
     }
     
   }
 }
//Ensuite faire un fit avec ces donnees. J ai pas compris si la fonction de minimizefactory ou on pouvait donner les donnees du fit a l avance nous sert la. Je pense que c etait peut etre que dans la class hypothesis on peut creeer un tableau qui nous sert vec les donnees gamma et phi 0 qui nous sert avant pour calculer l expected excess. 
