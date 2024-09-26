% Thermal conductivity Unit: W/(m K)
% Equation and Data from VDI Heat Atlas p. 302

function lambda = lambda_function(fc,T)


num_components=6;
%Order:  CH4          CO2       CO          H2           H2O         N2
A = [ 8.154     ,    -3.882   , -0.783    , 0.651     , 13.918   , -0.133   ];
B = [ 0.008     ,    0.053    , 0.103     , 0.767     , -0.047   , 0.101    ];
C = [ 0.351530  ,    0.071460 , -0.067590 , -0.687050 , 0.258066 , -0.060650];
D = [ -0.338650 ,    -0.070310, 0.039450  , 0.506510  , -0.183149, 0.033610 ];
E = [ 0.140920  ,    0.018090 , -0.009470 , -0.138540 , 0.055092 , -0.007100];

M_ch4=16.04; %g/mol
M_co2=44.01;
M_co=28.01;
M_h2=2.016;
M_h2o=18.015; 
M_n2=28.0134; 

M=[M_ch4 M_co2 M_co M_h2 M_h2o, M_n2];

lambda=zeros(1,num_components);


lambda=A*10^(-3) + B*10^(-3)*T + C*10^(-6)*T^2 +...
         D*10^(-9)*T^3 + E*10^(-12)*T^4;



F=zeros(num_components,num_components);
yF=zeros(num_components,1);
whole_term=zeros(num_components,1);
    
for i=1:num_components
    for j=1:num_components

        % to avoid calclulating the component with itself:
        if j~=i

            % if component j does not exist, F(i,j) is not calculated:
            if fc(j)~=0
            F(i,j)=(1+(lambda(i)/lambda(j))^(0.5)*(M(j)/M(i))^(0.25))^2/...
                sqrt(8+(1+M(i)/M(j)));
            end
        end
    end

     yF(i)=sum(fc(i)*F(i,:));

    %Avoid denominator beeing 0 if one component doesn't exist
    if yF(i)~=0
    whole_term(i)=(fc(i)*lambda(i))/yF(i);
    end
end

lambda=lambda*fc';
end