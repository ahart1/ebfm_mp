// Aggregate and Ecosystem MSY 

// estimate aggregate or ecosystem MSY using this code, which is estimated depends on wether data provided is aggregated by aggregate group or entire ecosystem, fit 1 model for each group

// SEARCH nrow and ncol to ID areas with problemsindexing (for loops
// figure out how to calculate Bmsy given MSY
#include <TMB.hpp>
  
  // Objective function which returns the value to minimize (often the negative log likelihood)
  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    ///// Data section /////
    DATA_IVECTOR(Biomass); //  Vector of biomass for aggregate group of interest across all simulations
    // DATA_IVECTOR(data2);
    
    ///// Parameter section /////
    PARAMETER(dummy); // Include dummy variable to debug code without full estimation of likelihood
    // PARAMETER(logit_r_par); // Growth parameter, may be single value if ecosystem ref point estimated instead of aggregate ref points
    PARAMETER(log_K_par); // Carrying cpacity parameter
    
    // Retransform variables so not in log space 
    // Type r_par = exp(logit_r_par)/(1+exp(logit_r_par));
    Type K_par = exp(log_K_par);


    // Local variables
    Type Bmsy; 
    
     
    
    
    ///// Calculate parameters (could also do in R and pass in as DATA) /////

    
    ///// Initialize objective function at zero /////
    Type obj_fun; // Always need objective function object
    obj_fun = 0;


    ///// Code which contributes to objective function /////
    // Generally make a prediction and compare that prediction to data (likelihood)
    // minimize comparison so you pick parameter values to best predict data
    // obj_fun += likelihood_component1 + likelihood_component2...
    obj_fun += dummy*dummy; // dummy objective function used to debug (make sure code compiles) without estimating all parameters
            
    // obj_fun += r_par*Biomass*(1 - (Biomass/K_par)) - Frate*Biomass; // want to minimize the change in biomass with respect to time (dBt/dt)


    // Bt+1 = Bt + rBt - rBt^2/K - Ct   OR   Bt+1 = Bt + rBt(1-Bt/K) - Ct // Schaefer production model
    // dBt/dt = rBt(1-Bt/K)             OR   dBt/dt = rBt - rBt^2/K       // derivative with respect to time
    // d(dB/dt)/dB = r(1- 2Bt/K)        OR   d(dB/dt)/dB = r - 2rBt/K     // solve fro Bmsy by setting to zero, since r is a multiplier, don't estimate here
    for (int isim=0; isim<Biomass.size(); isim++){
         obj_fun += (1 - Biomass(isim)/K_par)*(1 - Biomass(isim)/K_par); // r is a a multiplier constant which isn't used to estimate Bmsy so not estimated here, squared to ensure never negative, want obj_fun to be as close to zero as possible
    }
    
    
    ///// ADReport reports deviation /////
    // ADREPORT(r_par); // Report variables or parameters with deviation
    ADREPORT(K_par);

    
    ///// Report /////
    // REPORT(variable_or_parameter_name_here); // Report variable or parameter value
    
    
    ///// Return objective function /////
    return(obj_fun);
    
    ///// Advice if not compiling /////
    // Check that all lines end in ;
    // Check that * between multiplied objects in equations
    // Check indexing on objects in equations
    // Check indexing of storage objects (e.g. StorageVector = DataVector1[i]*DataVector2[i] won't work since this produces a single value but there is no indexing of StorageVector[i])
    // Check indexing of for() loops (indexing starts at 0)
    // Check spelling of objects in equations
    // Try commenting out parts of the line causing the error and gradually add them back in until it causes an error again to ID the source
  }

  