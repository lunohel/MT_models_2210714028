%density function
% VDI p. 131

%fc: Component fraction 
function rho = density_function(fc,T,p)

    v_mol=8.314*T/p; %m3/mol
    
    %Order: CH4 CO2 CO H2 H2O N2
    rho_c=[16.04/v_mol, 44.01/v_mol, 28.01/v_mol, 2.016/v_mol,  ...
         18.015/v_mol, 28.0134/v_mol]; %g/m^3

    rho=rho_c*fc'*10^(-3); %kg/m^3

end