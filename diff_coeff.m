% Gaseous diffusion coefficient in m2/s
% according to Fuller - form VDI
% p. 150


function D=diff_coeff(T,p)
  M_co2=44.01; %g/mol
  M_co=28.01;
  M_h2=2.016;
  M_ch4=16.04;
  M_h2o= 18.015; 
  M_n2= 28.0134; 

    
  %Diffusion Volume taken form Table 9 - VDI p. 150 
  d_nu_co2=26.9;
  d_nu_co=18;
  d_nu_h2=6.12;
  d_nu_ch4=15.9+4*2.31;
  d_nu_h2o=13.1;
  d_nu_n2=18.5;
  

  %m2/s
  D_ch4=(0.00143*T^1.75 * (M_ch4^(-1) + M_co2^(-1))^0.5)/...
      (p*10^(-5)*sqrt(2)* (d_nu_ch4^(1/3)+d_nu_co2^(1/3))^2)*10^(-4); 

  D_co2=D_ch4;

  D_h2=(0.00143*T^1.75 * (M_h2^(-1) + M_co2^(-1))^0.5)/...
      (p*10^(-5)*sqrt(2)* (d_nu_h2^(1/3)+d_nu_co2^(1/3))^2)*10^(-4);  

  D_co=(0.00143*T^1.75 * (M_co^(-1) + M_co2^(-1))^0.5)/...
      (p*10^(-5)*sqrt(2)* (d_nu_co^(1/3)+d_nu_co2^(1/3))^2)*10^(-4); 

  D_h2o=(0.00143*T^1.75 * (M_h2o^(-1) + M_co2^(-1))^0.5)/...
      (p*10^(-5)*sqrt(2)* (d_nu_h2o^(1/3)+d_nu_co2^(1/3))^2)*10^(-4); 

  D_n2=(0.00143*T^1.75 * (M_h2^(-1) + M_co2^(-1))^0.5)/...
      (p*10^(-5)*sqrt(2)* (d_nu_n2^(1/3)+d_nu_co2^(1/3))^2)*10^(-4); 

  D=[D_ch4 D_co2 D_co D_h2 D_h2o D_n2];
end