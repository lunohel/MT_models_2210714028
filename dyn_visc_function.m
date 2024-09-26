% Dynamic viscosity Unit: Pa s
% Equation and Data from VDI Heat Atlas p. 302
% Wilke mixing rule for low pressures (VDI p. 145)
% Does not work for pure substances!

function dyn_visc = dyn_visc_function(fc,T)


num_components=6;

%Order:     CH4        CO2       CO         H2      H2O         N2 
A =     [ -0.07759 ,-0.18024, 0.01384 , 0.18024 , 0.64966  , -0.01020 ];
B =     [ 0.50484  ,0.65989 , 0.74306 , 0.27174 , -0.15102 , 0.74785  ];
C =     [ -0.43101 ,-0.37108, -0.62996, -0.13395, 1.15935  , -0.59037 ];
D =     [ 0.03118  ,0.01586 , 0.03948 , 0.00585 , 0.10080  , 0.03230  ];
E =     [ -0.00981 ,-0.00300, -0.01032, -0.00104, 0.03100  , -0.00673 ];
 
M_ch4=16.04; %g/mol
M_co2=44.01;
M_co=28.01;
M_h2=2.016;
M_h2o= 18.015; 
M_n2= 28.0134; 

M=[M_co2 M_co M_h2 M_ch4 M_h2o, M_n2];

dyn_visc_sc=zeros(1,num_components);

for i=1:num_components
     dyn_visc_sc(i)=A(i)*10^(-5) + B(i)*10^(-7)*T + C(i)*10^(-10)*T^2 +...
              D(i)*10^(-12)*T^3 + E(i)*10^(-15)*T^4;
end

F=zeros(num_components,num_components);
yF=zeros(num_components,1);
whole_term=zeros(num_components,1);
    
for i=1:num_components
    for j=1:num_components

        % to avoid calclulating the component with itself:
        if j~=i

            % if component j does not exist, F(i,j) is not calculated:
            if fc(j)~=0
            F(i,j)=(1+(dyn_visc_sc(i)/dyn_visc_sc(j))^(0.5)*(M(j)/M(i))^(0.25))^2/...
                sqrt(8+(1+M(i)/M(j)));
            end
        end
    end
    yF(i)=sum(fc(i)*F(i,:));

    %Avoid denominator beeing 0 if one component doesn't exist
    if yF(i)~=0
    whole_term(i)=(fc(i)*dyn_visc_sc(i))/yF(i);
    end
end

dyn_visc=sum(whole_term);