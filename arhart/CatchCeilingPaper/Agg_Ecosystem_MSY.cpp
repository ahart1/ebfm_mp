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
    DATA_MATRIX(Biomass); // Matrix of biomasses where each column is a species group (or 1 for whole ecosystem) and rows are different simulations 
    DATA_INTEGER(ncol_Biomass); // Number of columns in Biomass matrix
    DATA_INTEGER(nrow_Biomass); // Number of rows in Biomass matrix
    // DATA_IVECTOR(data2);
    
    ///// Parameter section /////
    PARAMETER(dummy); // Include dummy variable to debug code without full estimation of likelihood
    PARAMETER_VECTOR(logit_r_par); // Growth parameter(s), may be single value if ecosystem ref point estimated instead of aggregate ref points
    PARAMETER_VECTOR(log_K_par); // Carrying cpacity parameter(s)
    PARAMETER_VECTOR(logit_Frate); // Fishing mortality rate(s)

    // Retransform variables so not in log space 
    vector<Type> r_par = exp(logit_r_par)/(1+exp(logit_r_par));
    vector<Type> K_par = exp(log_K_par);
    vector<Type> Frate = exp(logit_Frate)/(1+exp(logit_Frate));


    // Local variables
    vector<Type> Bmsy(nrow_Biomass); 
    
     
    
    
    ///// Calculate parameters (could also do in R and pass in as DATA) /////

    
    ///// Initialize objective function at zero /////
    Type obj_fun; // Always need objective function object
    obj_fun = 0;
    Type temp = 0;

    ///// Code which contributes to objective function /////
    // Generally make a prediction and compare that prediction to data (likelihood)
    // minimize comparison so you pick parameter values to best predict data
    // obj_fun += likelihood_component1 + likelihood_component2...
    obj_fun += dummy*dummy; // dummy objective function used to debug (make sure code compiles) without estimating all parameters
       
    for(int igroup=0; igroup<ncol_Biomass; igroup++){ // Loop over columns (may be aggregate groups or 1 column for whole ecosystem) = estimate for all agg groups at once
      for (int irow=0; irow<nrow_Biomass; irow++){ // Loop over all rows (in this case model simulations) = estimate parameters given info from all simulations at once
        
        obj_fun += r_par(igroup)*Biomass(irow,igroup)*(1 - (Biomass(irow,igroup)/K_par(igroup))) - Frate(igroup)*Biomass(irow,igroup); // want to minimize the change in biomass with respect to time (dBt/dt)
        //obj_fun += temp;

      }
    }
    // !!!!!!! estimate MSY for set biomass (avg over last 6 years for each simulation)
    // Then calculate Bmsy as K/2
    // Do this for all the 
    
    
    ///// ADReport reports deviation /////
    ADREPORT(r_par); // Report variables or parameters with deviation
    ADREPORT(K_par);
    ADREPORT(Frate);

    
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
  