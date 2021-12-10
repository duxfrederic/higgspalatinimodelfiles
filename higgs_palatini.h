#ifndef HIGGS_PALATINI_H //Usual macro guard to prevent multiple inclusion
#define HIGGS_PALATINI_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Frédéric Dux  Year: 2021

#include "CosmoInterface/cosmointerface.h"

//Include cosmointerface to have access to all of the library.

namespace TempLat
{
    /////////
    // Model name and number of fields
    /////////

    // In the following class, we define the defining parameters of your model:
    // number of fields of each species and the type of tinteractions.

    struct ModelPars : public TempLat::DefaultModelPars {
        static constexpr size_t NScalars = 1;
        // In our phi4 example, we only want 2 scalar fields.
        static constexpr size_t NPotTerms = 1;
        // Our potential naturaly splits into two terms: the inflaton potential
        // and the interaction with the daughter field.

        // All the numbers of fields are 0 by default, so we need only
        // to specify that we want two scalar fields.
        // See the model with gauge fields to have an example of how to turn
        // them on and specify interactions.
    };

  #define MODELNAME higgs_palatini
  // Here we define the name of the model. This should match the name of your file.

  template<class R>
  using Model = MakeModel(R, ModelPars);
  // In this line, we define an appropriate generic model, with the correct
  // number of fields, ready to be customized.
  // If you are curious about what this is doing, the macro is defined in
  // the "CosmoInterface/abstractmodel.h" file.

  class MODELNAME : public Model<MODELNAME>
  // Declaration of our model. It inherits from the generic model defined above.
  {
 //...
private:

   double Lambda, xi, sqrtxi;
// Here are the declaration of the model specific parameters. They are 'private'
// to force you using them only within your model and not outside.

// Some parameters which are declared in the class "Model" and which are useful (they are all 'public'):

// fldS0, piS0 : arrays which should contain the initial homogeneous values of
//               the scalar fields
//
// alpha, fStar, omegaStar : time and field rescaling to go to program units.
//
// fldS : The actual object which contains the scalar fields.

  public:

    MODELNAME(ParameterParser& parser, RunParameters<double>& runPar, std::shared_ptr<MemoryToolBox> toolBox): //Constructor of our model.
    Model<MODELNAME>(parser,runPar.getLatParams(), toolBox, runPar.dt, STRINGIFY(MODELLABEL)) //MODELLABEL is defined in the cmake.
    {
    
      	/////////
      	// Initial homogeneous components of the fields (read from parameters file, or specified here if not)
      	/////////

        fldS0 = parser.get<double, 1>("initial_amplitudes");
        piS0  = parser.get<double, 1>("initial_momenta", {0});
        
        /////////
        // Independent parameters of the model (read from parameters file)
        /////////       
            
    	xi     = parser.get<double>("xi");
        Lambda = parser.get<double>("Lambda");
        
        /////////
        // Dependent parameters of the model (dependent on previous ones)
        /////////  
    
    	sqrtxi = sqrt(xi);
        	
        /////////
        // Rescaling for program variables
        /////////
        
        fStar      = 2.435e18;// use the (reduced) planck mass instead of fldS0[0];
       	alpha      = 0;
       	omegaStar  = sqrt(Lambda) * fStar / (2 * xi);
    
       	setInitialPotentialAndMassesFromPotential();
   	}

   /////////
   // Program potential (add as many functions as terms are in the potential)
   /////////  
   
    auto potentialTerms(Tag<0>) // Inflaton potential energy
    {
      return  pow(tanh( sqrtxi * fldS(0_c)), 4);
    }

    
   /////////
   // Derivatives of the program potential with respect fields (add one function for each field)
   ///////// 

    auto potDeriv(Tag<0>) // Derivative with respect to the inflaton
    {
        auto tanarg = sqrtxi * fldS(0_c);
        
        auto num = 4 * sqrtxi * pow( tanh( tanarg ), 3);
        auto den = pow(cosh(tanarg), 2);
        
        return num / den;
    }

    
   /////////
   //  Second derivatives of the program potential with respect fields (add one function for each field)
   /////////

    auto potDeriv2(Tag<0>) // Second derivative with respect to the inflaton
    {
        auto tanarg = sqrtxi * fldS(0_c);
        
        auto num1   = 4 * xi * ( 4 - cosh(2*tanarg) );
        auto num2   = pow(tanh(tanarg), 2);
        auto den1   = pow(cosh(tanarg), 4);
        return num1 * num2  / den1;

    }
  

    };
}

#endif //HIGGS_PALATINI_H
