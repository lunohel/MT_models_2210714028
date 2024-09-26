
%fs: Component fraction, whole segment average
function u=vel_profile(T,p,fs,nz)
  
load cat_parameters.mat

if T<1123
cat_l=0.05;
end

ch4=fs(1);
co2=fs(2);
co=fs(3);
h2=fs(4);
h2o=fs(5);
n2=fs(6);

nm=18;
%Two cells are 
dx=sqrt(opencs_cell)/(nm);
dy=dx;

%Molar mass of components in g/mol
    M_ch4=16.04;
    M_co2=44.01; 
    M_co=28.01;
    M_h2=2.016;
    M_h2o= 18.015; 
    M_n2= 28.0134; 
    
    M=[M_ch4 M_co2 M_co M_h2 M_h2o M_n2];

    dz=cat_l/nz;
    
    sum_moles=p*(opencs_cell*dz)/(8.314*T);

    mdot=(30/3600*T/1123)/num_cells*density_function([0.25,0.75,0,0,0,0],T,p);

    vel_gas_channel=(mdot/density_function(fs,T,p))/opencs_cell;

v_old=zeros(nm,nm);

%inital guess, just average velocity, walls have no-slip condition
v_old(2:end-1,2:end-1)=vel_gas_channel;  

v_new=v_old;

den=density_function(fs,T,p);
visc=dyn_visc_function(fs,T);

Re=den*sqrt(opencs_cell)*...
    vel_gas_channel/visc;
    

dp=64/Re*cat_l/sqrt(opencs_cell)*den/2*vel_gas_channel^2;

alpha=1/visc*(dp/cat_l);

beta=2/dx^2+2/dy^2;

tol=10^(-9);
iter_max=40000;
iter=1;
convergence=0;


while convergence==0 && iter<=iter_max
    counter=0;
    for i=2:nm-1
        for j=2:nm-1 

            v_new(i,j)=((v_old(i,j+1)+v_old(i,j-1))/dx^2+...
                (v_old(i+1,j)+v_old(i-1,j))/dy^2+alpha)/beta;   

            test=(dx^2*dp)/(4*visc*cat_l)+(v_old(i,j+1)+v_old(i,j-1)+v_old(i+1,j)+v_old(i-1,j))/4;

            err=v_new(i,j)-v_old(i,j);
            err=(v_new(i,j)-v_old(i,j))/v_old(i,j);
            if abs(err)<=tol
            counter=counter+1;
            end 
        end
    end
    if counter==(i-1)*(j-1)
          convergence=1;
    end
    v_old=v_new;
    iter=iter+1;
end

u=v_new;

for i=1:nm-1
     for j=1:nm-1
        u_av(i,j)=sum(u(i:i+1,j:j+1),"all")/(4);
     end   
end



save velocity_guess.mat u

% u_av=zeros(nm/3,nm/3); 
% 
% ix=1:3:nm;
% iy=1:3:nm;
%  for i=1:nm/3
%      for j=1:nm/3
%         u_av(i,j)=sum(u(ix(i):ix(i)+2,iy(j):iy(j)+2),"all")/9;
%      end   
%  end

 u=u_av;

end