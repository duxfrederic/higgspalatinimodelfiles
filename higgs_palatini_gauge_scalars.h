#ifndef HIGGS_PALATINI_GAUGE_SCALARS_H //Usual macro guard to prevent multiple inclusion
#define HIGGS_PALATINI_GAUGE_SCALARS_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Frédéric Dux, Year: 2020

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
        static constexpr size_t NScalars = 4;
        // In our phi4 example, we only want 2 scalar fields.
        static constexpr size_t NPotTerms = 4;
        // Our potential naturaly splits into two terms: the inflaton potential
        // and the interaction with the daughter field.

        // All the numbers of fields are 0 by default, so we need only
        // to specify that we want two scalar fields.
        // See the model with gauge fields to have an example of how to turn
        // them on and specify interactions.
    };

  #define MODELNAME higgs_palatini_gauge_scalars
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

  double Lambda, xi, g, gz;
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

        fldS0 = parser.get<double, 1>("initial_amplitudes", {0});
        piS0  = parser.get<double, 1>("initial_momenta", {0});
        
        /////////
        // Independent parameters of the model (read from parameters file)
        /////////       
            
 	xi     = parser.get<double>("xi");
        Lambda = parser.get<double>("Lambda");
        g      = parser.get<double>("g"); 
        gz     = parser.get<double>("gz"); 
        
        /////////
        // Dependent parameters of the model (dependent on previous ones)
        /////////  
    
        	
        /////////
        // Rescaling for program variables
        /////////
        
        fStar      = 2.435e18;// use the (reduced) planck mass instead of fldS0[0];
       	alpha      = 0;
       	omegaStar  = sqrt(Lambda) * 2.435e18 / (2 * xi);
    
       	setInitialPotentialAndMassesFromPotential();
   	}

   /////////
   // Program potential (add as many functions as terms are in the potential)
   /////////  
   
    auto potentialTerms(Tag<0>) // Inflaton potential energy
    { 
      auto sharg = sqrt(xi) * fldS(0_c);
      return pow<4>(tanh(sharg));
    }
    
    auto potentialTerms(Tag<1>) // W+ potential energy
    {
      auto sharg = sqrt(xi) * fldS(0_c);
      return 0.5 * xi / Lambda * pow<2>(tanh(sharg)) / pow<2>(cosh(sharg)) * pow(g,2) * pow<2>(fldS(1_c));
    }

    auto potentialTerms(Tag<2>) // W- potential energy
    {
      auto sharg = sqrt(xi) * fldS(0_c);
      return 0.5 * xi / Lambda * pow<2>(tanh(sharg)) / pow<2>(cosh(sharg)) * pow(g,2) * pow<2>(fldS(2_c));
    }

    auto potentialTerms(Tag<3>) // Z0 potential energy
    {
      auto sharg = sqrt(xi) * fldS(0_c);
      return 0.5 * xi / Lambda * pow<2>(tanh(sharg)) / pow<2>(cosh(sharg)) * pow(gz,2) * pow<2>(fldS(3_c));
    }



    // Advanced note (ignore if you are satisfied with the above) :
    // - The 'auto' return type is important because the object returned is
    // not say an array containing  the value of the expression but the expression itself, which can and will be
    // evaluate later on. The type of the  expression itself depends on the expression and can be intricated. See
    // manual for more  details.
    // - The syntax 0_c is equivalent to Tag<0>(),
    // i.e. creating  an object of type 0. This operator '_c' is a modern C++ user-defined type literal,
    // taken from Boost and located in fcn/util/rangeiteration/tagliteral.h .



   /////////
   // Derivatives of the program potential with respect fields
   // (add one function for each field).
   /////////

    auto potDeriv(Tag<0>) // Derivative with respect to the inflaton
    {
      auto sharg = sqrt(xi) * fldS(0_c);
      auto globalfact =  0.5 * sqrt(xi) / Lambda * tanh(sharg) / pow<4>(cosh(sharg));
      auto inp11 = xi * (3 - cosh(2*fldS(0_c)));
      auto inp12 =   g*g   * ( pow<2>(fldS(1_c)) + pow<2>(fldS(2_c)) ) 
                   + gz*gz *   pow<2>(fldS(3_c));
      auto inp2 = 8 * Lambda * pow<2>(sinh(sharg));
      return globalfact * ( inp11 * inp12 + inp2);
    }
    
    auto potDeriv(Tag<1>) // Derivative with respect to W+
    {
      auto sharg = sqrt(xi) * fldS(0_c);
      return xi / Lambda * pow<2>(tanh(sharg)) / pow<2>(cosh(sharg)) * g*g * fldS(1_c);
    }
    
    auto potDeriv(Tag<2>) // Derivative with respect to W-
    {
      auto sharg = sqrt(xi) * fldS(0_c);
      return xi / Lambda * pow<2>(tanh(sharg)) / pow<2>(cosh(sharg)) * g*g * fldS(2_c);
    }

    auto potDeriv(Tag<3>) // Derivative with respect to Z0
    {
      auto sharg = sqrt(xi) * fldS(0_c);
      return xi / Lambda * pow<2>(tanh(sharg)) / pow<2>(cosh(sharg)) * gz*gz * fldS(3_c);
    }

    /////////
   //  Second derivatives of the program potential with respect fields
   // (add one function for each field)
   /////////

    auto potDeriv2(Tag<0>) // Second derivative with respect to the inflaton
    {
        // term 1
        auto sharg = sqrt(xi) * fldS(0_c);
        auto globalfactor = xi/ Lambda / pow<4>(cosh(sharg));

        auto bosonsum =   g*g * (pow<2>(fldS(1_c)) + pow<2>(fldS(2_c))) 
                        + gz*gz * pow<2>(fldS(3_c));

        auto term1 = ( cosh(2*sharg) - 10 * pow<2>(tanh(sharg)) ) * xi * bosonsum; 
        auto term2 = - 4 * Lambda * (cosh(2*sharg)-4) * pow<2>(tanh(sharg));

        return globalfactor * (term1 + term2);
    }
  
    auto potDeriv2(Tag<1>) // Derivative with respect to W+
    {
      auto sharg = sqrt(xi) * fldS(0_c);
      return xi / Lambda * pow<2>(tanh(sharg)) / pow<2>(cosh(sharg)) * g*g;
    }
    
    auto potDeriv2(Tag<2>) // Derivative with respect to W-
    {
      auto sharg = sqrt(xi) * fldS(0_c);
      return xi / Lambda * pow<2>(tanh(sharg)) / pow<2>(cosh(sharg)) * g*g;
    }

    auto potDeriv2(Tag<3>) // Derivative with respect to Z0
    {
      auto sharg = sqrt(xi) * fldS(0_c);
      return xi / Lambda * pow<2>(tanh(sharg)) / pow<2>(cosh(sharg)) * gz*gz;
    }

    };
}

#endif 
