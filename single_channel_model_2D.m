%2D single-channel model - Luca Nohel 09.2024
clear all
close all

% Indexing order of components: CH4 CO2 CO H2 H2O N2

%% Model Setup and Inputs

%Catalyst Properties
load cat_parameters.mat

%-------------------- Numerical setup -------------------------------------

%Turn on DR reaction? yes=1, no=0
dr=1;

%Turn on RWGS reaction? yes=1, no=0
rwgs=1;

%Want to calculate preheaters? (Turn off reactions)
heatup=0;

%Plot mass conservation? yes=1, no=0
conscond=1;


% Amount of segements in z-Direction
if heatup==1
    nz=200000;
    cat_l=0.05;
else
    nz=120000; 
end

% Number of cells per row in one segment
% Is now dependent on the number of cells used in velocity caluclations
nc=17; %needs to be uneven

%Diffusion coefficient scaling factor
diff_sf=1;

%Results from kinetic parameter optimization (take from thesis: (k1 K2 K3))
ks=[1 1 1];

%-------------------- Initial compositions and volumes --------------------

% Width of one cell in m
gs=sqrt(opencs_cell)/nc;

%Make the script callable from execution script
if exist('T_gas_in','var')==0 && exist('T_cat_in','var')==0
T_gas_in=1123; %850°C from meeting 19.3.24
T_cat_in=1123;
end

V_gas=30*T_gas_in/1123; %m3/h - calculated from 30 m3/h at 850°C
V_gas_channel=V_gas/(num_cells*3600); %m3/(s cell)
vel_gas_channel=(V_gas_channel)/opencs_cell; %m/s

p=101325; %Pa - Normal pressure

% Catalyst loading
cat_load=90000; %g/m3

%Initial gas composition
noc=6; %number of components

% Add molar fractions for feed composition (take from thesis)
ch4=0;
co2=0;
co=0;
h2=0;
h2o=0;
n2=0;

%Initial composition as a vector
com=[ch4 co2 co h2 h2o n2];

%Molar mass of components in g/mol
M_ch4=16.04;
M_co2=44.01; 
M_co=28.01;
M_h2=2.016;
M_h2o= 18.015; 
M_n2= 28.0134; 

M=[M_ch4 M_co2 M_co M_h2 M_h2o M_n2];

%Massflow constant:
ndot=((30/3600)*p)/(8.314*1123); %mol/s

massflow_init=(com*ndot)*M'*10^(-3);

%Stochiometric coefficients
nu_dr=[-1 -1 2 2 0 0];
nu_rwgs=[0 -1 1 -1 1 0];

%-------------------- Time and gird step ----------------------------------

%Segment length in m
dz=cat_l/nz; 
% z-Coordinate
z_coord=linspace(0,cat_l,nz); 

% Inital average time step needed for heat up
dt_m=dz/vel_gas_channel;

%Time steps for reaction and heat transfer in s
dt=zeros(nc,nz-1); 

%------------------- Thermodynamics ---------------------------------------

% Gibbs energy in J/mol and equilibrium constant
dG0_ch4=-50.5*10^3;
dG0_co=-137.2*10^3;
dG0_co2=-394.4*10^3;
dG0_h2o=-229*10^3;

dG_dr= -dG0_ch4 -dG0_co2 +2*dG0_co; 
dG_rwgs= -dG0_co2 + dG0_h2o + dG0_co;

% Eq constant at 273K
K298_dr=exp(-dG_dr/(8.314*298));
K298_rwgs=exp(-dG_rwgs/(8.314*298));

% Enthalpy of reaction in J/mol
dH0_ch4=-74.6*10^3;
dH0_co=-110.5*10^3;
dH0_co2=-393.5*10^3;
dH0_h2o=-242*10^3;

dH_dr=-dH0_ch4 - dH0_co2 + 2*dH0_co;
dH_dr=2.598946524696951e+05; %Hess' law
dH_rwgs=-dH0_co2 + dH0_co + dH0_h2o;
dH_rwgs=3.332939251623346e+04; %Hess' law

% Equilibrium constant T dependency
K_eqdr=exp(dH_dr/8.314*(1/298-1/T_gas_in))*K298_dr; 
K_eqdr=552.946; %Aspen
K_eqrwgs=exp(dH_rwgs/8.314*(1/298-1/T_gas_in))*K298_rwgs; 
K_eqrwgs=1.08057; %Aspen


%------------------------------- Matrices ----------------------------------

%Gas Temperature in K
if exist('T_gas_new','var')==0
T_gas=zeros(nc,nc,nz);
T_gas(:,:,1)=T_gas_in;
else
T_gas=zeros(nc,nc,nz);
T_gas(:,:,1)=T_gas_new;
end

%Velocity (needed for weighted average in results processing)
u=zeros(nc,nc,nz-1);

%Catalyst Temperature in K assumed to be uniform across walls
T_cat=zeros(nc,nz);
T_cat(:,1)=T_cat_in;

%Amount of substance in one channel and segment in mol/segment spilt
%equally to all cells - Ideal Gas
sum_moles=0;
for ncom=1:noc
    %Molar fractions for each component; using initial composition 
    f{ncom}=zeros(nc,nc,nz);
    f{ncom}(:,:,1)=com(ncom);
    %Molar fraction in whole segment needed for velocity profile
    fs(ncom)=com(ncom);
    
    %Partial pressures for each component; using initial composition 
    pp{ncom}=zeros(nc,nc,nz);
    pp{ncom}(:,:,1)=p*com(ncom);

    %Moles of each component per cell;
    n{ncom}=zeros(nc,nc,nz);
    n{ncom}(:,:,1)=((pp{ncom}(1,1,1)*(dz*opencs_cell))/(8.314*T_gas_in))/nc^2;
    %Sum of moles in one segment; Isobaric --> sum_moles = constant if T
    %constant!
    sum_moles=sum_moles+sum(n{ncom}(:,:,1),'all');
end

%Needed for heat up calculations
sum_moles_temperature=sum_moles;

sum_moles_cell=sum_moles/nc^2;

%------------------------ Heat balance set-up -----------------------------

% Reaction rate in mol/s needed for heat of reaction
r_dr=zeros(nc,nz-1);

% Heat of reaction in J
Q_r=zeros(nc,nz-1);

% Electric heat
Q_el=zeros(nc,nz-1);

% Overall energy balance cat (-Q_r+Q_el)
Q_cat=zeros(nc,nz-1);

% Heatup-power, sum of convection heat power
P_heatup=zeros(nc,nz-1);

% Power demand of reaction in W
P_r=zeros(nc,nz-1);

% Electric heating power in W
P_el=zeros(nc,nz-1);

check=zeros(1,nz-1);

%-------------------- Introduce Live Plots --------------------------------
%when stated before
if conscond==1
a=figure;
mc = animatedline;
set(a, 'WindowStyle', 'Docked');
grid on 
box on
xlim([0 cat_l])
xlabel('Catalyst length in m')
ylabel('Relative change of mass')
end

%% Calculations

% i is represents main loop with i indexing the segment
for i=1:nz-1

    %--------------------- New volume flow from new n ---------------------                           

    for ncom=1:noc
        fs(ncom)=sum(n{ncom}(:,:,i),'all')/sum_moles;
    end
     
    if i~=1 && heatup==1
        T_mean=(sum(u(:,:,i-1).*T_gas(:,:,i-1),"all")/sum(u(:,:,i-1),"all"));
    else
        T_mean=T_gas_in;
    end

    u(:,:,i)=vel_profile(T_mean,p,fs,nz); 

    %----------------------------------------------------------------------
    %-------------------------- Mass Transport ----------------------------    
    %----------------------------------------------------------------------

    %---------------------------- Diffusion -------------------------------
    
    %Diffusion coefficients in m2/s
    D=diff_coeff(mean(T_gas(:,:,i),'all'),p);
    D=D*diff_sf;

    %Diffusion matrix, so that new data does no influence diff.
    %calucaltions until they are completed:
    for ncom=1:noc
        n_diff{ncom}=zeros(nc,nc);
    end

    for j=2:nc-1
        for k=2:nc-1 
            for ncom=1:noc
                n_diff{ncom}(j,k)=n{ncom}(j,k,i)+(D(ncom)*dz)/(u(j,k,i)*gs^2)...
                    *(n{ncom}(j+1,k,i)+n{ncom}(j-1,k,i)+n{ncom}(j,k+1,i)...
                        +n{ncom}(j,k-1,i)-4*n{ncom}(j,k,i));
            end
        end
        
        % %Neumann Walls: dC/dx=0 only calculated for one side (symmetry)
        for ncom=1:noc
        n_diff{ncom}(j,1)=n{ncom}(j,1,i)+(D(ncom)*dz)/(u(j,1,i)*gs^2)*...
            (n{ncom}(j+1,1,i)+n{ncom}(j-1,1,i)+2*n{ncom}(j,1+1,i)-...
                4*n{ncom}(j,1,i));
        end
    end

    %Corner values: Double Neumann dC/dx=0 & dC/dy=0
    for ncom=1:noc
    n_diff{ncom}(1,1)=n{ncom}(1,1,i)+(D(ncom)*dz)/(u(1,1,i)*gs^2)*...
        (2*n{ncom}(1+1,1,i)+2*n{ncom}(1,1+1,i)-4*n{ncom}(1,1,i));
    %Fill whole matrix from symmetry
    n_diff{ncom}(end,1)=n_diff{ncom}(1,1);
    n_diff{ncom}(:,end)=n_diff{ncom}(:,1);
    n_diff{ncom}(1,:)=n_diff{ncom}(:,1);
    n_diff{ncom}(end,:)=n_diff{ncom}(:,1);

    %Replace amounts in segment i with diffusion matrix
    n{ncom}(:,:,i)=n_diff{ncom};
    end
    
    %---------------------------- Reaction --------------------------------
    
    %Amounts of substance per cell +- the amounts produced/used; ar stands
    %for "after reaction"
    for ncom=1:noc
        n_ar{ncom}=zeros(nc,nc);
    end

    %------------------------------ DRM -----------------------------------

    %Reaction happens at the outer cells
    for j=1:nc

    %DR on/off
    if dr==1
    k_r=1.06111e-08*exp(-1.0669e+04/(8.314*T_gas(1,j,i)))*ks(1); %mol/(g s Pa2) 
    % the g in the unit for reaction rate comes from catalyst loading
    % (...per gram of catalyst)
    k1=1.72*10^(-5)*exp(3970*4.184/(8.314*T_gas(1,j,i))); %Pa-1
    k4=1.19*10^(-7)*exp(5070*4.184/(8.314*T_gas(1,j,i)));
    % Adsorption of CO2 and CO onto CeZrO2 optimized
    k2=2.79*10^(-5)*exp(15100*4.184/(8.314*T_gas(1,j,i)))*ks(2);
    k3=5.1*10^(-7)*exp(3180*4.184/(8.314*T_gas(1,j,i)))*ks(3);
    
    gamma=1-((pp{3}(1,j,i)*10^(-5)/1.01325)^2*(pp{4}(1,j,i)*10^(-5)/1.01325)^2)...
        /(K_eqdr*(pp{1}(1,j,i)*10^(-5)/1.01325)*(pp{2}(1,j,i)*10^(-5)/1.01325));

    % Rate of change of CH4 - unit: mol/(g s)
    R_ch4=-(k_r*(pp{1}(1,j,i))*(pp{2}(1,j,i))*gamma)...
        /((1+k1*(pp{1}(1,j,i))+k4*(pp{4}(1,j,i)))*(1+k2*(pp{2}(1,j,i))+k3*(pp{3}(1,j,i))));
    
    % Rate of change of CH4 , unit: mol/s 
    % for each cell at the wall; at the corners, there are two sides with
    % cat --> double the amount of cat available in corners
        if j==1 || j==nc
            R_ch4=R_ch4*(opencs_cell*dz*cat_load)/((nc)*2);
        else
            R_ch4=R_ch4*(opencs_cell*dz*cat_load)/((nc)*4);
        end

    % Reaction rate (Rate of change of CH4/stochiometric coeff.) in mol/s
    r_dr(j,i)=R_ch4/-1;
    
    %DR off
    else
   
    r_dr(j,i)=0;

    end

    % dt for reaction is evaluated from the velocity in each cell - also
    % needed for heat transfer
    dt(j,i)=dz/u(j,1,i);
    
    % New amount of all substances after reaction step at all walls: 
    % These are the amounts of substance in each cell +/- the amounts 
    % produced/used during the reaction
        sum_moles_ar=0;
        for ncom=1:noc
        n_new=n{ncom}(j,1,i)+nu_dr(ncom)*r_dr(j,i)*dt(j,i);
        n_ar{ncom}(j,1)=n_new; 
        n_ar{ncom}(1,j)=n_new;
        n_ar{ncom}(end,j)=n_new; 
        n_ar{ncom}(j,end)=n_new;
        n_ar{ncom}(2:nc-1,2:nc-1)=n{ncom}(2:nc-1,2:nc-1,i);

        %Sum of moles after reaction
        sum_moles_ar=sum_moles_ar+sum(n_ar{ncom}(:,:),'all');
        end
    end
    
    for j=1:nc
        for k=1:nc
            sum_moles_cell_ar=0;
            for ncom=1:noc
            %Sum of moles for each cell after the reaction 
            sum_moles_cell_ar=n_ar{ncom}(j,k)+sum_moles_cell_ar;
            end

            %Component fractions
            for ncom=1:noc
            f{ncom}(j,k,i+1)=n_ar{ncom}(j,k)/sum_moles_cell_ar;
            %Partial pressures in Pa for reaction in i+1
            pp{ncom}(j,k,i+1)=f{ncom}(j,k,i+1)*p;
            end
            % ----------------------- RWGS --------------------------------
            % using the equilibrium approach
            if rwgs == 1
                w=(pp{3}(j,k,i+1)*10^(-5)/1.01325+pp{5}(j,k,i+1)*10^(-5)/1.01325+K_eqrwgs...
                    *(pp{2}(j,k,i+1)*10^(-5)/1.01325+pp{4}(j,k,i+1)*10^(-5)/1.01325))/(1-K_eqrwgs);
                q=(pp{3}(j,k,i+1)*10^(-5)/1.01325*pp{5}(j,k,i+1)*10^(-5)/1.01325-K_eqrwgs...
                    *pp{4}(j,k,i+1)*10^(-5)/1.01325*pp{2}(j,k,i+1)*10^(-5)/1.01325)/(1-K_eqrwgs);
            
                x(j,i)=(-w/2+sqrt((w/2)^2-q))*10^(5)*1.01325;
                                
                if pp{4}(j,k,i+1)-x(j,i)>=0
                    pp{2}(j,k,i+1)=pp{2}(j,k,i+1)-x(j,i);
                    pp{3}(j,k,i+1)=pp{3}(j,k,i+1)+x(j,i);
                    pp{4}(j,k,i+1)=pp{4}(j,k,i+1)-x(j,i);
                    pp{5}(j,k,i+1)=pp{5}(j,k,i+1)+x(j,i);
                else
                    x(j,i)=(-w/2-sqrt((w/2)^2-q))*10^(5)*1.01325;
                    pp{2}(j,k,i+1)=pp{2}(j,k,i+1)-x(j,i);
                    pp{3}(j,k,i+1)=pp{3}(j,k,i+1)+x(j,i);
                    pp{4}(j,k,i+1)=pp{4}(j,k,i+1)-x(j,i);
                    pp{5}(j,k,i+1)=pp{5}(j,k,i+1)+x(j,i);
                end
            
                f{2}(j,k,i+1)=pp{2}(j,k,i+1)/p;
                f{3}(j,k,i+1)=pp{3}(j,k,i+1)/p;
                f{4}(j,k,i+1)=pp{4}(j,k,i+1)/p;
                f{5}(j,k,i+1)=pp{5}(j,k,i+1)/p;
            
            end

            %Calculate the amounts of substance for each cell after reaction
            %fulfilling constant amounts of substance in each segment
            for ncom=1:noc
                 n{ncom}(j,k,i+1)=sum_moles_cell*f{ncom}(j,k,i+1);
            end
        end
    end

    %----------------------------------------------------------------------    
    %--------------------------- Heat Balance -----------------------------
    %----------------------------------------------------------------------

    % ---------------------------- Catalyst -------------------------------
    %Reactions
    for j=1:nc
    
        %Extent of reaction for DR:
        eor_dr=r_dr(j,i)*dt(j,i);
    
        %Extent of reaction for RWGS: reaching the equilibrium by solving 
        % for x and correcting the concentrations of the substances: 
        % x/p is the vol% of the amounts of substance reacting 
        if rwgs==1
        eor_rwgs=x(j,i)/p*sum_moles_cell;
        else
        eor_rwgs=0;
        end
        
        %Heat of reaction for every wall cell - catalyst side
        if j==1 || j==nc
        Q_r(j,i)=2*eor_dr*dH_dr+2*eor_rwgs*dH_rwgs; %J
        else 
        Q_r(j,i)=4*eor_dr*dH_dr+4*eor_rwgs*dH_rwgs; %J
        end

        %Heat power of reaction for cells (4 identical cells are summed)
        P_r(j,i)=Q_r(j,i)/dt(j,i); %W

        if heatup==1
        T_cat(:,:)=T_cat_in;        
        else
        % if Q_el is to be determined, it is initially set to 0 and then
        % calculated via the energy balance
        Q_el(j,i)=Q_r(j,i);
        P_el(j,i)=P_r(j,i);
        end
        
        % Cat heat balance
        Q_cat(j,i)=-Q_r(j,i)+Q_el(j,i);

    end

    %----------------------------- Radial conduction --------------------------
    % Save current gas temperature so that volume increase can be
    % calculated
    T_prev=mean(T_gas(:,:,i),'all');

    fc=zeros(1,noc);
    T_rc=T_gas(:,:,i);
    dt_rc=dz./u(:,:,i);


    for j=2:nc-1
        for k=2:nc-1
            for ncom=1:noc
                fc(ncom)=f{ncom}(j,k,i);
            end
            %Average values for between all concerned cells
            lambda=(lambda_function(fc,T_gas(j,k,i))+lambda_function(fc,T_gas(j+1,k,i))+...
                lambda_function(fc,T_gas(j-1,k,i))+lambda_function(fc,T_gas(j,k+1,i))+...
                lambda_function(fc,T_gas(j,k-1,i)))/5;
            den=(density_function(fc,T_gas(j,k,i),p)+density_function(fc,T_gas(j+1,k,i),p)+...
                density_function(fc,T_gas(j-1,k,i),p)+density_function(fc,T_gas(j,k+1,i),p)+...
                density_function(fc,T_gas(j,k-1,i),p))/5;
            cp=(cp_function(fc,T_gas(j,k,i))+cp_function(fc,T_gas(j+1,k,i))+cp_function(fc,T_gas(j-1,k,i))...
                +cp_function(fc,T_gas(j,k+1,i))+cp_function(fc,T_gas(j,k+1,i)))/5;
            
            dT=lambda/(den*cp*gs^2)*(T_gas(j+1,k,i)+T_gas(j-1,k,i)+T_gas(j,k+1,i)+...
                T_gas(j,k-1,i)-4*T_gas(j,k,i))*dt_rc(j,k);
            T_rc(j,k)=T_rc(j,k)+dT;  
        end
    end
            T_gas(:,:,i+1)=T_rc;

    
    
    %-------------------- Convective heat transfer ------------------------
    
    %Only if there is a temperature difference between gas and cat --> save
    %computation demand
    if T_cat(1,i)~=T_gas(1,1,i)
   
    %Get an array fc containing fractions of each components in the wall
    %cells - needed for heat transfer between wall cells and cat
    fc=zeros(nc,noc);
    %Fractions in each wall cell (position of cell, component)
    for ncom=1:noc
        fc(1:nc,ncom)=f{ncom}(1:nc,1,i);
    end
    %Average fraction over whole segment
    for ncom=1:noc
        fs(ncom)=sum(n{ncom}(:,:,i),'all')/sum_moles;
    end

    
    mdot=zeros(nc,1);
    cp_gas=zeros(nc,1);
   
    for j=1:nc
        for ncom=1:noc
            mdot(j)=(n{ncom}(j,1,1)*M(ncom)*10^(-3))/dt(j,i)+mdot(j);
        end
            cp_gas(j)=cp_function(fc(j,:),T_gas(j,1,i));
       
            A=sqrt(opencs_cell)/nc*dz;

            A_cell=dz*sqrt(opencs_cell)/nc;


       %Calculate the needed dimensionless numbers (Re and Prantl contains d_h:
       %needs to be calculatecd for whole channel, not for every cell)
       Re=(u(j,1,i)*dz*density_function(fc(j,:),T_gas(j,1,i),p))/dyn_visc_function(fc(j,:),T_gas(j,1,i));
       Pr=(cp_function(fc(j,:),T_gas(j,1,i))*dyn_visc_function(fc(j,:),T_gas(j,1,i)))/lambda_function(fc(j,:),T_gas(j,1,i));
       Nu=0.332*Re^(1/2)*Pr^(1/3);
       lambda=lambda_function(fc(j,:),T_gas(j,1,i));
       alpha=(Nu*lambda)/dz;

        totadv=0;
        totconv=0;
        totcond=0;
        iterations=0;
        Hbdif=1;
        tol=10^(-14);
        start=273;
        end1=T_cat_in-tol;
        To=(start+end1)/2;
        while abs(Hbdif)>tol || isnan(Hbdif)
           To=(start+end1)/2;
           
               advection=mdot(j)*cp_gas(j)*(T_gas(j,1,i)-To);   
               
               if j==1 || j==nc
               convection=alpha*A*(To-T_gas(j,1,i))/log((T_cat(j,i)-(T_gas(j,1,i)))/(T_cat(j,i)-To))*2;
               else
               convection=alpha*A*(To-T_gas(j,1,i))/log((T_cat(j,i)-(T_gas(j,1,i)))/(T_cat(j,i)-To));
               end

               if j==1 || j==nc
               To_cond=To;
               lambda1=(lambda_function(fc(j,:),T_gas(j,1,i))+lambda_function(fc(j,:),T_gas(j,2,i)))/2; 
               lambda2=(lambda_function(fc(j,:),T_gas(j,1,i))+lambda_function(fc(j,:),T_gas(2,1,i)))/2;
               conduction=-lambda1*A_cell/gs*(To_cond-T_gas(1,2,i))-lambda2*A_cell/gs*(To_cond-T_gas(1,2,i));
               else
               To_cond=To;
               lambda1=(lambda_function(fc(j,:),T_gas(j,1,i))+lambda_function(fc(j,:),T_gas(j-1,1,i)))/2; 
               lambda2=(lambda_function(fc(j,:),T_gas(j,1,i))+lambda_function(fc(j,:),T_gas(j+1,1,i)))/2;
               lambda3=(lambda_function(fc(j,:),T_gas(j,1,i))+lambda_function(fc(j,:),T_gas(j,2,i)))/2;      
               conduction=-lambda1*A_cell/gs*(To_cond-T_gas(j-1,1,i))-lambda2*A_cell/gs*(To_cond-T_gas(j+1,1,i))...
                  -lambda3*A_cell/gs*(To_cond-T_gas(j,2,i));
               end

           convcheck=-convection+conduction;

           Hbdif=(advection+convection+conduction);
             

           if Hbdif>=0 || isnan(convection)
              start=To;
           else
              end1=To;
           end
           iterations=iterations+1;

        if abs((To-T_cat_in)/T_cat_in)<=tol
            advection=0;
            Hbdif=0;
            if j==1 || j==nc
            To_cond=To;
            lambda1=(lambda_function(fc(j,:),T_gas(j,1,i))+lambda_function(fc(j,:),T_gas(j,2,i)))/2; 
            lambda2=(lambda_function(fc(j,:),T_gas(j,1,i))+lambda_function(fc(j,:),T_gas(2,1,i)))/2;
            conduction=-lambda1*A_cell/gs*(To_cond-T_gas(1,2,i))-lambda2*A_cell/gs*(To_cond-T_gas(2,1,i));
            convection=conduction;
            else
            To_cond=To;
            lambda1=(lambda_function(fc(j,:),T_gas(j,1,i))+lambda_function(fc(j,:),T_gas(j-1,1,i)))/2; 
            lambda2=(lambda_function(fc(j,:),T_gas(j,1,i))+lambda_function(fc(j,:),T_gas(j+1,1,i)))/2;
            lambda3=(lambda_function(fc(j,:),T_gas(j,1,i))+lambda_function(fc(j,:),T_gas(j,2,i)))/2;
            conduction=-lambda1*A_cell/gs*(To_cond-T_gas(j-1,1,i))-lambda2*A_cell/gs*(To_cond-T_gas(j+1,1,i))...
                  -lambda3*A_cell/gs*(To_cond-T_gas(j,2,i));
            convection=conduction;
            end
        end
           
        end  

        totcond=totcond+conduction;
        totadv=totadv+advection;
        totconv=totconv+convection;
        check(i)=totcond+totadv+totconv;
        
        if j==1 || j==nc
             P_heatup(j,i)=abs(convection)*2;
        else
             P_heatup(j,i)=abs(convection)*4;
        end

        T_rc(j,1)=To;
        T_rc(j,end)=To;
        T_rc(1,j)=To;
        T_rc(end,j)=To;
        if j==1
        T_rc(end,end)=To;    
        end

        

    end    
        
            %New number of moles per substance from temperature increase
            for j=1:nc
                for k=1:nc
                    for ncom=1:noc
                        f{ncom}(j,k,i+1)=f{ncom}(j,k,i);
                        n{ncom}(j,k,i+1)=sum_moles_cell*f{ncom}(j,k,i+1);
                    end
                end
            end
    end

    if heatup==1
        T_gas(:,:,i+1)=T_rc;

         %New volume and velocity form temperature change 
            deltaT_gas=mean(T_gas(:,:,i+1),'all')-T_prev;
        
            V_dot_new=((mean(T_gas(:,:,i+1),'all')*sum_moles_temperature*8.314)/p)/dt_m;
            vel_gas_channel=V_dot_new/opencs_cell;
        
            %New sum of moles coming from temperature increase
            sum_moles=(p*opencs_cell*dz)/(8.314*mean(T_gas(:,:,i+1),'all'));
            sum_moles_cell=sum_moles/nc^2;
        
    else 
        % If no heat is tranferred from cat to gas, tempearture just
        % remains constant
        T_gas(:,:,i+1)=T_gas(:,:,i);
        T_cat(:,i+1)=T_cat(:,i);
    end
    

    %-------------------------- Live Plots --------------------------------

    if conscond==1 && ismember(i,[2:20:nz])
    massflow_comp=zeros(1,noc);
    massflow_prev=zeros(1,noc);
        for ncom=1:noc
            massflow_comp(ncom)=sum((n{ncom}(:,:,i)*M(ncom))./(dz./u(:,:,i)),"all");
            massflow_prev(ncom)=sum((n{ncom}(:,:,1)*M(ncom))./(dz./u(:,:,1)),"all");
        end
        massflow_comp=sum(massflow_comp);
        massflow_prev=sum(massflow_prev);

    massflow=(sum(massflow_comp,"all")*num_cells*10^(-3))/(sum(massflow_prev,"all")*num_cells*10^(-3));
    
    addpoints(mc,z_coord(i),massflow);
    drawnow
    end

 end

