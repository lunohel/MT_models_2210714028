% 1D coupon reactor model - Luca Nohel 09.24.
clear all
close all

% RWGS on = 1 / off = 0;
rwgs=1;

% Scaling factors from kinetic modeling (take from thesis) (k K2 K3)
sf=[1 1 1];

T=1123; % K Reactor/Gas Temperature
p=101325; % Pa

% Catalyst mass taken from experiment
mass_cat=0; % g of Catalyst on coupon

% Inlet volume flows
V=[100 200 800]; %ml/min at standard conditions
V_cold=V;
V=V*(T/293)/60*10^(-6); %m3/s at T reaction (calculated form ideal gas)

% Number of segments 
nz=1000;

% --------------------------- DR Kinetics ---------------------------------

% Dry Reforming
k=1.06111e-08*exp(-1.0669e+04/(8.314*T))*sf(1); %mol/(g s Pa) 
% the g in the unit for reaction rate comes from catalyst loading
% (...per gram of catalyst)
k1=1.72*10^(-5)*exp(1.6610e+04/(8.314*T)); %Pa-1
k4=1.19*10^(-7)*exp(2.1213e+04/(8.314*T));
k2=2.79*10^(-5)*exp(6.3178e+04/(8.314*T))*sf(2);
k3=5.1*10^(-7)*exp(1.3305e+04/(8.314*T))*sf(3);


%-------------------------------------------

% Gibbs energy in J/mol and equilibrium constant
dG0_ch4=-50.5*10^3;
dG0_co=-137.2*10^3;
dG0_co2=-394.4*10^3;
dG0_h2o=-229*10^3;

dG_dr= -dG0_ch4 -dG0_co2 +2*dG0_co; 
dG_rwgs= -dG0_co2 + dG0_h2o + dG0_co;

% Eq constant at 273K
K273_dr=exp(-dG_dr/(8.314*273));
K273_rwgs=exp(-dG_rwgs/(8.314*273));

% Enthalpy of reaction in J/mol
dH0_ch4=-74.6*10^3;
dH0_co=-110.5*10^3;
dH0_co2=-393.5*10^3;
dH0_h2o=-242*10^3;

dH_dr=-dH0_ch4 - dH0_co2 + 2*dH0_co;
dH_rwgs=-dH0_co2 + dH0_co + dH0_h2o;

% Equilibrium constant T dependency
K_eqdr=exp(dH_dr/8.314*(1/273-1/T))*K273_dr; 
K_eqdr=552.946; %Aspen Plus REquil
K_eqrwgs=exp(dH_rwgs/8.314*(1/273-1/T))*K273_rwgs; 
K_eqrwgs=1.08057; %Aspen Plus REquil


% Volume fractions inlet gas (take form thesis)
ch4=0;
co2=0;
n2=0;
co=0;
h2=0;
h2o=0;

M_ch4= 16.04; %g/mol
M_co2= 44.01; %g/mol
M_co = 28.01; %g/mol
M_h2 = 2.016; %g/mol
M_h2o= 18.015; %g/mol
M_n2= 28.0134; %g/mol

M=[M_ch4 M_co2 M_co M_h2 M_h2o M_n2];
M=M.*10^(-3); %kg/mol

p_ch4=p*ch4;
p_co2=p*co2;
p_co=p*co;
p_h2=p*h2;
p_n2=p*n2;
p_h2o=p*h2o;

% Geometry
A_r=0.019*0.003; %m2
coupon_l=0.047; %m
dz=coupon_l/nz; %m

%segment
n_segment=p*(A_r*coupon_l/nz)/(8.314*T);

n_ch4=zeros(nz,length(V));
n_co2=zeros(nz,length(V));
n_h2=zeros(nz,length(V));
n_co=zeros(nz,length(V));
n_n2=zeros(nz,length(V));
n_h2o=zeros(nz,length(V));

n_ch4(1,:)=n_segment*ch4;
n_co2(1,:)=n_segment*co2;
n_co(1,:)=n_segment*co;
n_h2(1,:)=n_segment*h2;
n_n2(1,:)=n_segment*n2;
n_h2o(1,:)=n_segment*h2o;

v=zeros(nz,length(V));
v(1,:)=V/A_r; %m/s Velocity through reactor
v_inlet=v;

h2co=zeros(1,length(V));

for l=1:length(V)

    m_dot=V(l)*density_function([ch4 co2 co h2 h2o n2],T,p);

    delta_z=0;
    % Time step per segment
    dt=(dz/v(1,l));%s

    for i=1:nz-1

    p_ch4=p*ch4;
    p_co2=p*co2;
    p_co=p*co;
    p_h2=p*h2;
    p_n2=p*n2;
    p_h2o=p*h2o;


    %---------------------------- Dry Reforming ---------------------------
    
    %Difference from equilibrium needs to be calculated in atm
    % 1atm = 1.01325bar
    gamma=1-((p_co*10^(-5)/1.01325)^2*(p_h2*10^(-5)/1.01325)^2)...
        /(K_eqdr*p_ch4*10^(-5)/1.01325*p_co2*10^(-5)/1.01325);
     
    % Rate of change of CH4 - unit: mol/(g s)
    R_ch4=-(k*p_ch4*p_co2*gamma)/((1+k1*p_ch4+k4*p_h2)*...
        (1+k2*p_co2+k3*p_co));

    r_dr=(R_ch4/-1)*(mass_cat/nz); %mol/s
    
    % Amount of substance increases with DR (ar... after reaction)
    n_ch4_ar=n_ch4(i,l)-r_dr*dt;
    n_co2_ar=n_co2(i,l)-r_dr*dt;
    n_co_ar=n_co(i,l)+2*r_dr*dt;
    n_h2_ar=n_h2(i,l)+2*r_dr*dt;
    n_n2_ar=n_n2(i,l);
    n_h2o_ar=n_h2o(i,l);

    n_sum_ar=n_ch4_ar+n_co2_ar+n_co_ar+n_h2_ar+n_n2_ar+n_h2o_ar;
    
    ch4=n_ch4_ar/n_sum_ar;
    co2=n_co2_ar/n_sum_ar;
    co=n_co_ar/n_sum_ar;
    h2=n_h2_ar/n_sum_ar;
    n2=n_n2_ar/n_sum_ar;
    h2o=n_h2o_ar/n_sum_ar;

    v(i+1,l)=(m_dot/density_function([ch4 co2 co h2 h2o n2],T,p))/A_r;
    dt=(dz/v(i+1,l));

    p_ch4=p*ch4;
    p_co2=p*co2;
    p_co=p*co;
    p_h2=p*h2;
    p_n2=p*n2;
    p_h2o=p*h2o;

    %------------------------------- RWGS ---------------------------------
    
    % Assume that RWGS is in equilibrium as the rate is much larger than
    % the DR rate
    
    if rwgs == 1
    % atm because of K_eq calculations
    w=(p_co*10^(-5)/1.01325+p_h2o*10^(-5)/1.01325+K_eqrwgs*...
        (p_co2*10^(-5)/1.01325+p_h2*10^(-5)/1.01325))/(1-K_eqrwgs);
    q=(p_co*10^(-5)/1.01325*p_h2o*10^(-5)/1.01325-K_eqrwgs*...
        p_h2*10^(-5)/1.01325*p_co2*10^(-5)/1.01325)/(1-K_eqrwgs);

    % from atm to Pa
    x1=(-w/2+sqrt((w/2)^2-q))*10^(5)*1.01325;
    x2=(-w/2-sqrt((w/2)^2-q))*10^(5)*1.01325;
    
    if p_h2-x1>=0
        p_co2=p_co2-x1;
        p_co=p_co+x1;
        p_h2=p_h2-x1;
        p_h2o=p_h2o+x1;
    else
        p_co2=p_co2-x2;
        p_co=p_co+x2;
        p_h2=p_h2-x2;
        p_h2o=p_h2o+x2;
    end

    ch4=p_ch4/p;
    co2=p_co2/p;
    co=p_co/p;
    h2=p_h2/p;
    n2=p_n2/p;
    h2o=p_h2o/p;

    end

    n_ch4(i+1,l)=n_segment*ch4;
    n_co2(i+1,l)=n_segment*co2;    
    n_co(i+1,l)=n_segment*co;
    n_h2(i+1,l)=n_segment*h2;
    n_n2(i+1,l)=n_segment*n2;
    n_h2o(i+1,l)=n_segment*h2o;

    end
    
    h2co(l)=h2/co;
    
end

for l=1:length(V)
X_ch4(l)=(v(1,l)*n_ch4(1,l)-v(end,l)*n_ch4(end,l))./(v(1,l)*n_ch4(1,l));
end

X_ch4
h2co