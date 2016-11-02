// SCHAEFER MODEL FOR MS_PROD OUTPUT
// SINGLE SPECIES
// GAVIN FAY
// 9 Feb 2012

DATA_SECTION

 init_int Nsp
 init_int rFase
 init_vector r_init(1,Nsp)
 init_int KFase
 init_vector K_init(1,Nsp)
 init_int zFase
 init_vector z_init(1,Nsp)
 init_int thetaFase
 init_vector theta_init(1,Nsp)
 init_int Fyear
 init_int Lyear
 init_matrix ObsCat(Fyear,Lyear,1,Nsp)
 init_int NBio
 !!int ncol=Nsp+1;
 init_matrix ObsBio(1,NBio,1,ncol)
 init_matrix ObsCV(1,NBio,1,ncol)


 int i
 int iyr
 int iobs


PARAMETER_SECTION

 init_number dummy(-1);

 init_bounded_vector log_r(1,Nsp,-10.1,0.4,rFase)
 init_bounded_vector log_K(1,Nsp,-10,20.1,KFase)
 init_bounded_vector log_z(1,Nsp,0.00001,10,zFase)
 init_bounded_vector logit_theta(1,Nsp,-10,0.,thetaFase)

 vector r(1,Nsp)
 vector K(1,Nsp)
 vector z(1,Nsp)
 vector theta(1,Nsp)
 vector q(1,Nsp)
 vector sigma(1,Nsp) 

 sdreport_matrix Bio(Fyear,Lyear+1,1,Nsp)

 objective_function_value objfun

PRELIMINARY_CALCS_SECTION

 log_r = log(r_init);
 log_K = log(K_init);
 log_z = log(z_init);
 logit_theta = log(elem_div(theta_init,(1.-theta_init)));

PROCEDURE_SECTION

 dvariable fpen1;
 objfun += square(dummy);

 //parameter values
 r = mfexp(log_r);
 K = mfexp(log_K);
 z = mfexp(log_z);
 theta = elem_div(mfexp(logit_theta),(1.+mfexp(logit_theta)));

 //cout << r << endl;
 //cout << K << endl;
 //cout << theta << endl;

 //get_biomass
 Bio(Fyear) = elem_prod(theta,K);
 for (i=1;i<=Nsp;i++)
  {
 for (iyr=Fyear+1;iyr<=Lyear+1;iyr++)
  {
   Bio(iyr,i) = posfun(Bio(iyr-1,i)*(1.+r(i)*(1.-pow((Bio(iyr-1,i)/K(i)),z(i)-1))),1.,fpen1);
   dvariable sr = 1.-ObsCat(iyr-1,i)/Bio(iyr,i);
   dvariable kcat=ObsCat(iyr-1,i);
   objfun+=1000*fpen1;
   if(sr< 0.001)
    {
     dvariable fpen=0.;
     kcat=Bio(iyr-1,i)*posfun(sr,0.001,fpen);
     objfun+=10000*fpen;
     // cout << " kludge "<<iy <<" "<<kcat<<" "<<cat(iy)<<" "<<fpen<<endl;
    }
    Bio(iyr,i)-=kcat;
 //   cout << i << " " << iyr << " " << Bio(iyr,i) << endl;
  }
 }
 

 //objective function value
 for (i=1;i<=Nsp;i++)
  {
   dvariable temp1 = 0;
  //dvariable temp2 = 0;
   int temp3 = 0;
  // for (iobs=1;iobs<=NBio;iobs++)
  //  {
    //if (ObsBio(iobs,1+1)!=-99)
  //   {
  //    temp1 += log(ObsBio(iobs,i+1)/(0.5*(Bio(ObsBio(iobs,1),i)+Bio(ObsBio(iobs,1)+1,i))));
  //  temp2 += 1./square(ObsCV(iobs,i));
   //   temp3 += 1;
   //  }
   // }
  //q(i) = mfexp(temp1/temp3);
  for (iobs=1;iobs<=NBio;iobs++)
   {
    if (ObsBio(iobs,i+1)!=-99)
  //   objfun += log(ObsCV(iobs,i+1))+0.5*square(log(ObsBio(iobs,i+1)/(q(i)*Bio(ObsBio(iobs,1),i))))/square(ObsCV(iobs,i+1));
    {
    temp1 += square(log(ObsBio(iobs,i+1)/(0.5*(Bio(ObsBio(iobs,1),i)+Bio(ObsBio(iobs,1)+1,i))))); 
    temp3 += 1;
   }
   }
   cout << temp1 << " " << temp3 << endl;
   sigma(i) = pow(temp1/temp3,0.5);
//   cout << sigma(i) << endl;
   objfun += temp3*log(sigma(i)) + temp3/2;
//     objfun += log(ObsCV(iobs,i+1))+0.5*square(log(ObsBio(iobs,i+1)/(0.5*(Bio(ObsBio(iobs,1),i)+Bio(ObsBio(iobs,1)+1,i)))))/square(ObsCV(iobs,i+1));
     //objfun += 0.5*square(log(ObsBio(iobs,i+1)/(0.5*(Bio(ObsBio(iobs,1),i)+Bio(ObsBio(iobs,1)+1,i)))))/square(ObsCV(iobs,i+1));
//    cout << i << " " << iobs << " " << objfun << endl;
  // }
  }
  cout << objfun << endl;


REPORT_SECTION

 report << r << endl;
 report << K << endl;
 report << z << endl;
 report << theta << endl;
 report << q << endl;
 report << sigma << endl;
 for (iyr=Fyear;iyr<=Lyear+1;iyr++)
  report << iyr << " " << Bio(iyr) << endl;

 