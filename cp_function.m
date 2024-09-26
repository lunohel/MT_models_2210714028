%Cp_gases Unit: J/(kg K)
% Equation and Data from VDI Heat Atlas p. 302

function cp = cp_function(fc,T)


num_components=6;
       

%Order:    CH4          CO2       CO          H2         H2O         N2
A = [ 1530.8043 , 514.5073 , 407.9796  , 392.8422 , 706.3032 , 432.2027 ];
B = [ 4.2038    , 3.4923   , 3.5028    , 2.4906   , 5.1703   , 3.5160   ];
C = [ -16.6150  , -0.9306  , 2.8524    , -3.6262  , -6.0865  , 2.8021   ];
D = [ -3.5668   , -6.0861  , -2.3018   , -1.9624  , -6.6011  , -4.1924  ];
E = [ 43.0563   , 54.1586  , 32.9055   , 35.6197  , 36.2723  , 42.0153  ];
F = [ -86.5507  , -97.5157 , -100.1815 , -81.3691 , -63.0965 , -114.2500];
G = [ 65.5986   , 70.9687  , 106.1141  , 62.6668  , 46.2085  , 111.1019 ];
 
M_ch4= 16.04; %g/mol
M_co2= 44.01; %g/mol
M_co = 28.01; %g/mol
M_h2 = 2.016; %g/mol
M_h2o= 18.015; %g/mol
M_n2= 28.0134; %g/mol

R_ch4=8.314/(M_ch4*10^(-3)); %J/(kg K)
R_co2=8.314/(M_co2*10^(-3)); %J/(kg K)
R_co=8.314/(M_co*10^(-3));   %J/(kg K)
R_h2=8.314/(M_h2*10^(-3));   %J/(kg K)
R_h2o=8.314/(M_h2o*10^(-3)); %J/(kg K)
R_n2=8.314/(M_n2*10^(-3)); %J/(kg K)

R_spez=[R_ch4 R_co2 R_co R_h2 R_h2o R_n2];

%calculate cp of each component

cp_sc=zeros(1,num_components);

for i=1:num_components
    
    cp_sc(i) = B(i) + (C(i)-B(i))*(T/(A(i)+T))^2 * ...
         (1- A(i)/(A(i)+T)* (D(i) + E(i)* T/(A(i)+T) + ...
           F(i)*(T/(A(i)+T))^2 + G(i)*(T/(A(i)+T))^3));
       
    cp_sc(i) = cp_sc(i)*R_spez(i);
    
end

cp = cp_sc*fc'; %J/(kg K)

end





