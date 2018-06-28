function spectral_2FGEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-field gyrofluid model to study the Dimits shift.                            %
% Based on a trunctation of Bruce Scott's GEM model.                            %
%                                                                               %
% Can use Wang and Hahm's formula for the neoclassical polarization of          %
% zonal flows                                                                   %
%                                                                               %                            
% Refs:                                                                         %
%   Dimits et al. Phys. Plasmas (2000)                                          %
%   Scott Phys. Plasmas (2005)                                                  %
%   Wang and Hahm, Phys. Plasmas (2007)                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
mkdir('newCode_tau1'); cd('newCode_tau1');
mkdir('plots');
delete('log.txt');

%CBC
RLn = 2.2;
rR = 0.18;
q = 1.4;
shat = 0.8;
ep = 2*rR;
Lna = 1/(RLn*ep);
LB = 2*(1-shat)*Lna;
Ld = 0.10;
eta=8.0;
scale=1;            % quick scale of the linear terms
mu=0;
nu=1e-6;
neo=true;

%LT=3.0; %scale*sqrt(2/pi);               % inverse temperature gradient scale
%ep=scale*0.3;            % inverse aspect ratio
%Ld=0.3;
h=3;                % hyperviscosity factor
hz=1;                % hyperviscosity factor
tau=1;          % ion temperature
hm=1;
%eta=LT;
LX=2*pi*20;      % X scale
LY=2*pi*20;      % Y scale
NX_real=256;     % resolution in x
NY_real=256;     % resolution in y
dt=1e-4;        % time step. Should start small as CFL updated can pick up the pace
max_dt=1e-2;    % Max time step after CFL adjustment
zonal=true;     % include zonal flow physics
TF=75000.0;   % final time
iF=1500000;  % final iteration, whichever occurs first
i_report=500;
TSCREEN=2000; % sreen update interval time (NOTE: plotting is usually slow)
initial_condition='random w';   %'simple vortices' 'vortices' 'random' or 'random w' 
AB_order=3; % Adams Bashforth order 1,2,3, or 4 (3 more accurate, 2 possibly more stable)
linear_term='exact'; % CN, BE, FE, or exact
simulation_type='NL'; % NL, QL, or L (nonlinear, quasilinear and linear respectively)
padding = true; % 3/2 padding, otherwise 2/3 truncation. Probably should always be true
save_plots = true; % save plots to file
cfl_cadence=5;
safety=0.3;
iRST=20000;
diagnostics=false;
pade=false;


rng('shuffle');
s=rng;

% print log file.
diary('log.txt') 
diary on;

fig1=figure(1);
set(fig1,'Position',[0, 0, 1200, 900]);

fig2=figure(2);



% ensure parameters get printed  to log.
%fprintf('mu: %.05e  nu:%.05e lT:%.05e lW:%.05e \n',mu,nu,landau_T,landau_n);
%fprintf('muh:%.05e L_T:%.05e ep:%.05e tau:%.05e\n',muh, LT,ep,tau);
%fprintf('nuz:%.05e\n',nuz);
%fprintf('LX:%.02f LY:%.02f NX:%d NY:%d\n',LX, LY, NX_real, NY_real);
%fprintf('h:%d hz:%d scale:%d Tf:%.01f iF:%d\n',h, hz, scale, TF, iF);
%fprintf('Nonlinear:%s padding:%d zonal_flows:%d\n',simulation_type, padding,zonal);
%fprintf('random seed:%d AB order:%d CFL step:%d\n',s.Seed, AB_order, cfl_cadence);

energyFile=0; 
if(strcmp(initial_condition,'restart'))
  energyFile = fopen('energy.dat','a');
else
  energyFile = fopen('energy.dat','w');
end

fprintf(energyFile,'## [1]time   [2]energy [3]enstrophy [4]ZF en [5]DW en [6]flux\n');

NX=NX_real;
NY=NY_real;
if(padding)
  NX=3*NX/2;
  NY=3*NY/2;
end

dx=LX/NX;
dy=LY/NY;
dkx=2*pi/LX;
dky=2*pi/LY;
LXnum=0:dx:LX;
LYnum=0:dy:LY;

minkx= -(NX_real/2 - 1)*dkx;
maxkx=  (NX_real/2 - 1)*dkx;
minky= -(NY_real/2 - 1)*dky;
maxky=  (NY_real/2 - 1)*dky;
kxnum= minkx:dkx:maxkx;
kynum= minky:dky:maxky;

I=sqrt(-1);

kx=dkx*ones(1,NY)'*(mod((1:NX)-ceil(NX/2+1),NX)-floor(NX/2));% matrix of wavenumbers in x direction 
ky=dky*(mod((1:NY)'-ceil(NY/2+1),NY)-floor(NY/2))*ones(1,NX);% matrix of wavenumbers in y direction 

% Cutting of frequencies using the 2/3 rule. 
% Discard asymmetric N/2 term. Otherwise reality is not enforced.
dealias=abs(kx/dkx) <1/3*NX & abs(ky/dky) <1/3*NY;


zonal_part = zeros(NY,NX);
fluct_part = ones(NY,NX);

% No damping on zonal modes. Modified Poisson equation (proper adiabatic electron response)
if(zonal)
  zonal_part(1,:) = zonal_part(1,:) + ones(1,NX);
  fluct_part = fluct_part - zonal_part;   
end


ksquare=kx.^2+ky.^2; 
ksquare(1,1) = 1;% Laplacian in Fourier space
b = tau.*ksquare;


G0 = 1.0./(1.0 + b);
G1 = 1.0./(1.0 + 0.5*b);
G2 = - 0.5* b ./ (1 + 0.5*b).^2;
if(~pade)
  bsm = (b<100).*b;
  G0 = (b<100).*(exp(-bsm).*besseli(0,bsm)) + (b>=100)./sqrt(2*pi*b);
  G1 = sqrt(G0);
  G2 = (b<100).*(0.5*bsm.*exp(-bsm).*(besseli(1,bsm) - besseli(0,bsm))./G1) - 0.25*(b>=100)./(2*pi*b).^0.25;
  clearvars bsm
end

ksquare_poisson = (fluct_part + (1-G0)/tau);
if(tau == 0)
  ksquare_poisson = fluct_part + ksquare;
else


if(neo)
  ksqbanana = b*q^2/rR^2;
  Gtr = 0.915966./sqrt(pi*ep*ksqbanana);
  Gp = 1.0./(2*ep*sqrt(2*pi*ksqbanana));
  Gtrp = 2*Gtr/pi;
  Gpp = 2*Gp/pi;
  ichi = 1.0./(1.83.*ep^1.5.*ksqbanana) + (1 + sqrt(8*ep)/pi.*Gtrp + (1-sqrt(8*ep)/pi).*Gpp)./(1.0+b) ...
      + sqrt(0.5*pi^3.*b).*(1 + sqrt(8*ep)/pi.*Gtr + (1-sqrt(8*ep)/pi).*Gp).*b./(1+b);
  chi=1.0./ichi;

  ksquare_poisson(1,:) =  (1-G0(1,:))/tau + chi(1,:);
  clearvars kbanan Gtr Gtrp Gp Gpp ichi chi;
end

ksquare_poisson(1,1)=1;  % fixed Laplacian in Fourier space for Poisson's equation
clearvars b

%kmu = fluct_part.*(mu+muh./ksquare);             % drag force array
%ky_abs = abs(ky);                           % Landau-like damping
%ksquare_viscous=nu*(ksquare.^h);            % Viscosity array
%ksquare_viscous=nu*(ksquare.^h).*fluct_part;            % Viscosity array
%ksquare_viscous(1,:)=nuz*(ksquare(1,:).^hz);            % Viscosity array




  
end

T11 = 1;
T12 = 0;
T21 = 0;
T22 = 1;
L11 = (I.*ky.*hm./ksquare_poisson + Ld*abs(ky) + nu.*ksquare).*fluct_part;
L12 = (-I*ky.*LB).*fluct_part;
L21 = (I*ky.*eta./ksquare_poisson).*fluct_part;
L22 = (Ld*abs(ky) + nu.*ksquare).*fluct_part;


phi_n = G1./ksquare_poisson;
phi_T = G2./ksquare_poisson;
if(tau > 0)
 % T11 = (fluct_part - (G0-1)/tau)./G1;
 % T12 = -G2./G1;
 % T21 = 0;
 % T22 = 1;
  T11 = 1.0;
  T12 = 0;
  T21 = 0;
  T22 = 1.0;
  L11 = fluct_part.*(I.*ky.*hm.*(G1 + eta*G2).*phi_n - I*ky.*LB.*((2*G1 + G2).*phi_n + 2) + Ld*abs(ky) + nu.*ksquare.^h);
  L12 = fluct_part.*(I.*ky.*hm.*(G1 + eta*G2).*phi_T - I*ky.*LB.*((2*G1 + G2).*phi_T + 2)); %- Ld*abs(ky).*G2./G1;
  L21 = fluct_part.*(I.*ky.*hm.*(G2 + eta*(G1+2*G2)).*phi_n - I*ky.*LB.*((G1 + 4*G2).*phi_n + 1)); %+ Ld*abs(ky).*G2;
  L22 = fluct_part.*(I.*ky.*hm.*(G2 + eta*(G1+2*G2)).*phi_T - I*ky.*LB.*((G1 + 4*G2).*phi_T + 4) + Ld*abs(ky)+ nu.*ksquare.^h);
end

clearvars phi_n phi_t 

tDet = 1.0./(T11.*T22 - T12.*T21);
disc = sqrt((L22.*T11 - L21.*T12 - L12.*T21 + L11.*T22).^2 - 4*(L12.*L21 - L11.*L22).*(T12.*T21 - T11.*T22));
eig_un = 0.5.*tDet.*(-L22.*T11 + L21.*T12 + L12.*T21 - L11.*T22 + disc);
eig_st = 0.5.*tDet.*(-L22.*T11 + L21.*T12 + L12.*T21 - L11.*T22 - disc);

maxfreq = 1/max(max(abs(eig_un(:))),max(abs(eig_st(:))));
if(max_dt > maxfreq) 
  max_dt = maxfreq;
  if(dt > max_dt)
    dt = max_dt;
  end
end

clearvars tDet disc

i=0;
k=0;
t=0.;
n_hat=0; %this is our field.
T_hat=0;
cn1=0;
cn2=0;
cn3=0;
cT1=0;
cT2=0;
cT3=0;
dt1=0;
dt2=0;
dt3=0;
u=0;
v=0;

% Define initial vorticity distribution
switch lower(initial_condition)
  case {'vortices'}
    [i,j]=meshgrid(1:NX,(1:NY));
     psi=exp(-((i*dkx-pi).^2+(j*dky-pi+pi/4).^2)/(0.2))+exp(-((i*dx-pi).^2+(j*dky-pi-pi/4).^2)/(0.2))-0.5*exp(-((i*dkx-pi-pi/4).^2+(j*dky-pi-pi/4).^2)/(0.4));
  case {'simple vortices'}
     [i,j]=meshgrid((1:NX)*(2*pi/NX),(1:NY)*(2*pi/NY));
     psi=1*sin(i).*cos(j).^2;    
  case {'random'}
     psi=1e-3*(2*rand(NX,NY) - 1);%normally 1e-3
     n_hat=ksquare.*fft2(psi);
  case {'random w'}
    n=5e-2*(2*rand(NY,NX)-1);%normally 5e-2
    n_hat=fft2(n);

      
    n = ifft2(n_hat);
    %w = w + 20*cos(15*2*pi*i);
    n_hat = fft2(n);
      
    T=5e-2*(2*rand(NY,NX)-1);%normally 5e-2
    T_hat=fft2(T);
    
    n_hat(1,:)=zeros(1,NX);
    T_hat(1,:)=zeros(1,NX);
      
    %n = ifft2(n_hat);
    %n = n + 0.2*cos(3*2*pi*i + pi);
    %n_hat = fft2(n);
      
  case {'restart'}
   readdump('restart.bin');
  case {'waves'}
    [i,j]=meshgrid((1:NX)*(2*pi/NX),(1:NY)*(2*pi/NY));
    n=20*sin(32*j) + 1e-1*cos(13*i);
    %w=2e1*sin(32*j) + 5e-2*(2*rand(NY,NX)-1);
    n_hat=fft2(n);
    % w_hat(1,:)=zeros(1,NY);
  otherwise
    disp('Unknown initial conditions !!!');
    return
end

if(padding) % ensure there's no energy in the padded region of k-space.
  n_hat=dealias.*n_hat;
  T_hat=dealias.*T_hat;
end
n_hat(1,1)=0; % Gauge condition, though not entirely necessary.
T_hat(1,1)=0;
n_hat0=n_hat; % Keep initial conditions for possible future diagnostics.

% create matrix components
calc_matrix();

if(save_plots) %plot growth rate contours
  plotGrowth();
end

while t<TF && i<iF
  psi_hat = (G1.*n_hat+ G2.*T_hat)./ksquare_poisson;  % Solve Poisson's Equation
  if(any(isnan(n_hat(:)))) % Something really bad has happened.
    fprintf('Divergence at iteration: %d\n',i);
    return;
  end
  
  if(mod(i,i_report)==0)
    fprintf('iteration: %d    dt:%.02e     t:%.03e \n',i,dt,t); 
    diary off; %flush diary
    diary on;
  end
     
  if (mod(i,iRST) == 0)
    dump('restart.bin');
  end
     
  %if(t > nextScreen && with_plotting)
    
  if (mod(i,TSCREEN)== 0) 
    plotfunc();
    nextScreen = t + TSCREEN;
  end
   
  % Compute the potential and get the velocity and gradient of vorticity
  conv_hat = 0;
    
  switch upper(simulation_type)
    case {'NL'} %full non-linear
      uhat = -I*ky.*G1.*psi_hat;
      vhat =  I*kx.*G1.*psi_hat;
      n_xhat = I*kx.*n_hat;
      n_yhat = I*ky.*n_hat;
      T_xhat = I*kx.*T_hat;
      T_yhat = I*ky.*T_hat;
            
            
      % dealiasing here truncates if not padded, other it has no effect
      u  =real(ifft2(dealias.*uhat));      % Compute  y derivative of stream function ==> u
      v  =real(ifft2(dealias.*vhat));      % Compute -x derivative of stream function ==> v
      n_x=real(ifft2(dealias.*n_xhat));      % Compute  x derivative of vorticity
      n_y=real(ifft2(dealias.*n_yhat));
      T_x=real(ifft2(dealias.*T_xhat));
      T_y=real(ifft2(dealias.*T_yhat));
      conv_n     = (u.*n_x + v.*n_y);         % evaluate the convective derivative (u,v).grad(w)   
      conv_T     = (u.*T_x + v.*T_y);         % evaluate the convective derivative (u,v).grad(w) 
            
      if(tau > 0 )
        u2hat = -I*ky.*G2.*psi_hat;
        v2hat =  I*kx.*G2.*psi_hat;          
        u2=real(ifft2(dealias.*u2hat)); 
        v2=real(ifft2(dealias.*v2hat));
        conv_n2 = (u2.*T_x + v2.*T_y);
        conv_T2 = (u2.*n_x + v2.*n_y) + 2*(u2.*T_x + v2.*T_y);
        conv_n = conv_n + conv_n2;
        conv_T = conv_T + conv_T2;
      end
            
      conv_n_hat = fft2(conv_n);              % go back to Fourier space
      conv_T_hat = fft2(conv_T);              % go back to Fourier space
 
    case {'QL'} %quasi-linear
      %for zonal part
      uphat = -ky.*fluct_part.*psi_hat;
      vphat =  kx.*fluct_part.*psi_hat;
        
      %for fluctuation part
      U_hat = kx.*zonal_part.*psi_hat;
      U_xxhat = (kx.^2).*U_hat;
      w_yhat = ky.*fluct_part.*n_hat;       
        
      % dealiasing here truncates if not padded, other it has no effect
      u_f  =real(ifft2(dealias.*uphat));      % Compute  y derivative of stream function ==> u
      v_f  =real(ifft2(dealias.*vphat));      % Compute -x derivative of stream function ==> v
      U    =real(ifft2(dealias.*U_hat));      % Compute zonal velocity
      U_xx =real(ifft2(dealias.*U_xxhat));      % Compute zonal velocity
      n_y  =real(ifft2(dealias.*w_yhat));
      conv_fluct = U.*n_y + u_f.*U_xx;         % evaluate the convective derivative (u,v).grad(w)   
      conv_zonal = u_f.*v_f;         % evaluate the convective derivative (u,v).grad(w)       
      conv_fluct_hat = fft2(conv_fluct);               % go back to Fourier space
      conv_zonal_hat = fft2(conv_zonal);
      conv_hat = fluct_part.*conv_fluct_hat + (kx.^2).*zonal_part.*conv_zonal_hat;
    case {'L'}%full linear 
            %do nothing
    otherwise
      disp('Unknown simulation type.');
      return
  end    
   
  if(padding)
    conv_n_hat = dealias.*conv_n_hat;
    conv_T_hat = dealias.*conv_T_hat;
  end
    
  if(mod(i+1,cfl_cadence)==0) % compute new timestep from CFL condition.
    abs_u=abs(u);
    abs_v=abs(v);
    maxV= max(abs_u(:))/dx + max(abs_v(:))/dy;
    new_dt=0;
    if(maxV>0)
      new_dt=1/maxV;
    else
      new_dt=inf;
    end
    new_dt_corr=min(max(0.5*dt,min(safety*new_dt,1.1*dt)),max_dt);
    if(dt ~= new_dt_corr)
      dt=new_dt_corr;
      % create matrix components
      calc_matrix()
    end
  end
    
  AB1=1.0; AB2=0; AB3=0; AB4=0;
  
  if (i < 1 || AB_order == 1) %F orward-Euler to generate history. Should run a small time-step at first. 
  elseif (i < 2 || AB_order == 2)
    w1=dt1/dt;
    AB1=(1.0 + 0.5*w1);
    AB2=0.5*w1;
  elseif (i < 3 || AB_order == 3)
    w1=dt1/dt;
    w2=dt2/dt;
    AB1=(2.0 + 3.0*w2 + 6.0*w1*(1+w1+w2))/(6.0*w1*(w1+w2));
    AB2=(2.0 + 3.0*w1 + 3.0*w2)/(6.0*w1*w2);
    AB3=(2.0 + 3.0*w1)/(6.0*w2*(w1+w2));
  elseif (i < 4 || AB_order == 4)
    w1=dt1/dt;
    w2=dt2/dt;
    w3=dt3/dt;
    AB1=(3.0 + 8.0*w2 + 4.0*w3 + 6.0*(2*w1^3 + w2*(w2+w3) +2.0*w1*(1.0+w2)*(1.0+w2+w3) + w1^2*(3.0 + 4.0*w2+2*w3)));
    AB1=AB1/(12.0*w1*(w1+w2)*(w1+w2+w3));
    AB2=(3.0 + 6.0*(w1+w2)*(w1+w2+w3) + 4.0*(2.0*(w1+w2)+w3))/(12.0*w1*w2*(w2+w3));
    AB3=(3.0 + 6.0*w1*(w1+w2+w3)+4.0*(2.0*w1 + w2+w3))/(12.0*w2*(w1+w2)*w3);
    AB4=(3.0 + 6.0*w1*(w1+w2)+4.0*(2.0*w1 + w2))/(12.0*w3*(w2+w3)*(w1+w2+w3));
  end
  
    
    
  switch upper(linear_term) %deprecated for now
    case {'CN'} %crank nicolson
     %       L1=(1 + 0.5*dt*lin_growth)./(1 - 0.5*dt*lin_growth);
     %       L2=dt./(1 - 0.5*dt*lin_growth);    
    case {'BE'} %backward euler
     %       L1=1.0./(1 - dt*lin_growth);
     %       L2=dt./(1 - dt*lin_growth);    
    case {'EXACT'}%exact
     %       L1=exp(dt*lin_growth);
     %       L2=dt*L1;   
     %       %conv_hat = exp(-t*lin_growth).*conv_hat;
     %       B=B.*exp((dt1)*lin_growth);
     %       C=C.*exp((dt1+dt2)*lin_growth);
     %       D=D.*exp((dt1+dt2+dt3)*lin_growth);
    otherwise
      disp('Unknown linear handling type.');
      return;
  end
    
  Qn = C11.*n_hat + C12.*T_hat - dt*(AB1*conv_n_hat - AB2*cn1 + AB3*cn2 - AB4*cn3);
  QT = C21.*n_hat + C22.*T_hat - dt*(AB1*conv_T_hat - AB2*cT1 + AB3*cT2 - AB4*cT3);
    

  % Compute Solution at the next step
  % Implicitly solve the linear term with 2nd order Cracnk-Nicholson
  w_hat_new = B11.*Qn + B12.*QT;
  n_hat_new = B21.*Qn + B22.*QT;
  
  t=t+dt;
  i=i+1;
        
  cn3=cn2;
  cn2=cn1;
  cn1=conv_n_hat;

  cT3=cT2;
  cT2=cT1;
  cT1=conv_T_hat;
    
  dt3=dt2;
  dt2=dt1;
  dt1=dt;
    
  n_hat=dealias.*w_hat_new;  
  T_hat=dealias.*n_hat_new;
  n_hat(1,1)=0;
  T_hat(1,1)=0;
    
  if(mod(i,1000)==0) % ensure reality
    n_hat=(fft2(real(ifft2(n_hat))));
    T_hat=(fft2(real(ifft2(T_hat))));
  end
    
end
fclose(energyFile);

%write final state 
fprintf('Simulation ended at time %.03e and iteration %d\n',t,i);
fileID=fopen('final_state.bin','w');
fwrite(fileID,n_hat);
fclose(fileID);

if(strcmp(upper(simulation_type),'L') && diagnostics)
  w0f=w_hat0.*exp(t*lin_growth); 
  comp=abs(n_hat-w0f)./abs(w0f);
  fprintf('Maximum linear error:%e\n',max(comp(:)));
end

  function calc_matrix()      
  % create matrix components
    A11 = T11 + 0.5*dt*L11;
    A12 = T12 + 0.5*dt*L12;
    A21 = T21 + 0.5*dt*L21;
    A22 = T22 + 0.5*dt*L22;
    C11 = T11 - 0.5*dt*L11;
    C12 = T12 - 0.5*dt*L12;
    C21 = T21 - 0.5*dt*L21;
    C22 = T22 - 0.5*dt*L22;
    D= A11.*A22 - A12.*A21;
    B11= A22./D;
    B12=-A12./D;
    B21=-A21./D;
    B22= A11./D;
  end


  function plotGrowth()
    lin_growth = real(eig_un);
    subplot(1,1,1);
    plotg=0;
    if(padding)
      plotg=circshift(dealias.*lin_growth,[NY_real/2,NX_real/2]); 
      plotg=plotg(2:NY_real,2:NX_real);
    else
      plotg=circshift(dealias.*lin_growth,[NY_real/2,NX_real/2]); 
      plotg=plotg(2:NY_real,2:NX_real);
    end
    grow=real(plotg);
    numel(find(grow(:)>0))
    max(grow(:))

    imagesc(kxnum,kynum,real(plotg)), axis equal tight, colorbar;
    set(gca,'Ydir','Normal')
    colormap(jet)
    title('growth rates');    
    xlabel('kx');
    ylabel('ky');
    if(save_plots)
      saveas(gcf,'growth.png');
    end
    drawnow
  end
  function plotfunc()
           % Go back in real space omega in real space for plotting
    diary off; %flush diary
    diary on;
    w_curr=n_hat;
    psi_curr=psi_hat;
    psi=real(ifft2(psi_curr));
    w=real(ifft2(w_curr)); 
    T=real(ifft2(T_hat));
    psi_y=real(ifft2(I*ky.*psi_curr));
        
    enstrophy = 0.5*(T_hat-w_curr).*conj(T_hat-w_curr);
    energy = 0.5*real(conj(psi_curr).*w_curr+ T_hat.*conj(T_hat));
    flux = psi_y.*T;
        
    enstrophy_tot = sum(enstrophy(:))/(NX*NY)^2;
    energy_tot    = sum(energy(:))/(NX*NY)^2; 
    flux_tot = sum(flux(:))/(NX*NY);
    ZF_energy = sum(sum(zonal_part.*energy))/(NX*NY)^2;
    DW_energy = sum(sum(fluct_part.*energy))/(NX*NY)^2;
        
    fprintf(energyFile,'%e %e %e %e %e %e\n',t,energy_tot,enstrophy_tot,ZF_energy,DW_energy, flux_tot);
             
    if(padding)
      enstrophy=circshift(enstrophy,[NY_real/2,NX_real/2]); 
      enstrophy=enstrophy(2:NY_real,2:NX_real);
      energy=circshift(energy,[NY_real/2,NX_real/2]); 
      energy=energy(2:NY_real,2:NX_real);
    else
      enstrophy=circshift(enstrophy,[NY_real/2,NX_real/2]); 
      enstrophy=enstrophy(2:NY_real,2:NX_real);
      energy=circshift(energy,[NY_real/2,NX_real/2]); 
      energy=energy(2:NY_real,2:NX_real);
    end
      
    wlog=max(log10(enstrophy),-10);
    energylog=max(log10(energy),-10);
    colorbar;
    set(0,'CurrentFigure',fig1);
    subplot(2,3,1)
    imagesc(LXnum,LYnum,psi), axis equal tight,c=colorbar();
    set(gca,'Ydir','Normal')
    colormap(jet)
    title(sprintf('psi t=%.02f',t));
    xlabel('x');
    ylabel('y');
    %ax = gca; axpos = ax.Position; cpos = c.Position;cpos(3) = 0.5*cpos(3); c.Position = cpos; ax.Position = axpos;
    subplot(2,3,2)
    imagesc(LXnum,LYnum,w), axis equal tight,c=colorbar();
    title(sprintf('vorticity t=%.02f',t));
    set(gca,'Ydir','Normal')
    xlabel('x');
    ylabel('y');
    %ax = gca; axpos = ax.Position; cpos = c.Position;cpos(3) = 0.5*cpos(3); c.Position = cpos; ax.Position = axpos;
    subplot(2,3,3);
    imagesc(LXnum,LYnum,T), axis equal tight,c=colorbar();
    title(sprintf('temperature t=%.02f',t));
    set(gca,'Ydir','Normal')
    xlabel('x');
    ylabel('y');
    subplot(2,3,4);
    %ax = gca; axpos = ax.Position; cpos = c.Position;cpos(3) = 0.5*cpos(3); c.Position = cpos; ax.Position = axpos;

    imagesc(kxnum,kynum, energylog), axis equal tight,c=colorbar();
    title('log10(energy power)');    
    set(gca,'Ydir','Normal')
    xlabel('kx');
    ylabel('ky');
    ax = gca; axpos = ax.Position; cpos = c.Position;cpos(3) = 0.5*cpos(3); c.Position = cpos; ax.Position = axpos;
    h=subplot(2,3,5);
    imagesc(kxnum,kynum,wlog), axis equal tight,c=colorbar();
    set(gca,'Ydir','Normal')
    xlabel('kx');
    ylabel('ky');
    title('log10(enstrophy power)');    
    %ax = gca; axpos = ax.Position; cpos = c.Position;cpos(1) = 0.5*cpos(1); c.Position = cpos; ax.Position = axpos;

    h=subplot(2,3,6);
    imagesc(kxnum,kynum,wlog), axis equal tight,c=colorbar();
    set(gca,'Ydir','Normal')
    xlabel('kx');
    ylabel('ky');
    title('log10(temperature power)');    
        
     %ax = gca; axpos = ax.Position; cpos = c.Position;cpos(3) = 0.5*cpos(3); c.Position = cpos; ax.Position = axpos;

    if(save_plots)
      % saveas(gcf,sprintf('plots/fig_%d.png',k));
      figs=gcf;
      figs.PaperUnits = 'inches';
      figs.PaperPosition=[0 0 16 8];
      print(sprintf('plots/fig_%d.png',k),'-dpng');
    end
    drawnow
    set(0,'CurrentFigure',fig2);
    kx_spec=energylog(NY_real/2,:);
    ky_spec=energylog(:,NX_real/2);
    kx_spec(NX_real/2)=0.5*(kx_spec(NX_real/2-1)+kx_spec(NX_real/2+1)); 
    ky_spec(NY_real/2)=0.5*(ky_spec(NY_real/2-1)+ky_spec(NY_real/2+1)); 
    [~,Ix] = max(kx_spec(:));
    [~,Iy] = max(ky_spec(:));
    if(diagnostics)
      disp(sprintf('%.10e %.10e',energy_tot,enstrophy_tot));
      disp(sprintf('Maximum DW mode:%.01f    Maximum ZF mode:%.01f',kynum(Iy), kxnum(Ix)));
    end
    plot(kxnum,kx_spec,kynum,ky_spec);
    legend('Zonal (ky=0)','DW (kx=0)');
    title(sprintf('DW/Zonal energy spectrum t=%.02f',t));
    xlabel('k');
    if(save_plots)
      saveas(gcf,sprintf('plots/k_spec_%d.png',k));
    end
    drawnow
    k=k+1; 
  end

  function dump(filename)
    fileID=fopen(filename,'w');
    fwrite(fileID,t,'double');
    fwrite(fileID,dt,'double');
    fwrite(fileID,i,'int');
    fwrite(fileID,k,'int');
    fwrite(fileID,dt1,'double');
    fwrite(fileID,dt2,'double');
    fwrite(fileID,dt3,'double');
    fwrite(fileID,real(cn1),'double'); fwrite(fileID,imag(cn1),'double');
    fwrite(fileID,real(cn2),'double'); fwrite(fileID,imag(cn2),'double');
    fwrite(fileID,real(cn3),'double'); fwrite(fileID,imag(cn3),'double');
    fwrite(fileID,real(cT1),'double'); fwrite(fileID,imag(cT1),'double');
    fwrite(fileID,real(cT2),'double'); fwrite(fileID,imag(cT2),'double');
    fwrite(fileID,real(cT3),'double'); fwrite(fileID,imag(cT3),'double');
    fwrite(fileID,real(n_hat),'double');fwrite(fileID,imag(n_hat),'double');
    fwrite(fileID,real(T_hat),'double');fwrite(fileID,imag(T_hat),'double');
    fclose(fileID);
  end

  function readdump(filename)
    fileID=fopen(filename,'r');
    t=fread(fileID,1,'double');
    dt=fread(fileID,1,'double');
    i=fread(fileID,1,'int');
    k=fread(fileID,1,'int');
    dt1=fread(fileID,1,'double');
    dt2=fread(fileID,1,'double');
    dt3=fread(fileID,1,'double');
      
    cn1=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    cn2=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    cn3=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    cT1=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    cT2=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    cT3=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    n_hat=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    T_hat=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    fclose(fileID);  
  end
end
