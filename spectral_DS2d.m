function spectral_DS2d(HM_in,TH_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1-field toy model of the Dimits shift. Similar to the Hasegawa-Mima           %
% equation or the Kuramoto-Sivashinsky equation.                                %
%                                                                               %
% Based off an MIT code originally made by Jean-Christophe Nave.                %
% Modified by Denis St-Onge                                                     %
%                                                                               %
% Laplacian(phi) = w                                                            %
% u = phi_y                                                                     %
% v =-phi_x                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear dto lg M r Q1 Q2 f1 f2 f3 isreal; 
clearvars -except HM_in TH_in;
clear ETDRK4; clear ETDRK3; clear ETDRK2;
basename='beta_scale/delt1.5_visc0.01/';
if(nargin > 0)
   basename=['TH-HW-',num2str(HM_in)]; %basename_in; 
end  
if(nargin > 1)
   basename=['d-',num2str(TH_in),'-h-',num2str(HM_in)]; %basename_in; 
end  

%c_map=parula;
cm_redblue =flip(cbrewer('div', 'RdBu',129));
%cm_magma=magma();
cm_inferno=inferno();
%cm_plasma=plasma();
%cm_viridis=viridis();
%cm_redblue=redblue();
c_map=cm_inferno;
c_maprb=cm_redblue;
%if(7==exist(basename,'dir'))
%  return;
%end
mkdir(basename);
mkdir([basename 'plots']); mkdir([basename 'data']);
delete([basename 'log.txt']);
scale=1;                 % quick scale of the linear terms
mu   =scale*[0];         % friction
nu   =scale*[0.01,0.0];      % viscosity
muZF =scale*[0.0];       % zonal friction
nuZF =scale*[0.0,0];         % zonal viscosity
l    =scale*0.0;           % Landau-like damping
gamma=scale*0.0;         % linear drive
rboost=1;
HM=scale*1;           % HM-type wave
delta_0=scale*1.5;    % Terry-Horton i delta
tau=0;                % FLR effects
neo=1;                % Neoclassical polarization factor 
neo1=-0.55;  % -0.55
neo2=4.0;    % 4.0
forcing=0.0;              % forcing magnitude
LX=2*pi*10;          % X scale
LY=2*pi*10;         % Y scale
NX_real=128;         % resolution in x
NY_real=128;         % resolution in y
dt=5e-2;            % time step. Should start small as CFL updated can pick up the pace
pert_size=1e-3;     % size of perturbation
TF=15000.0;           % final time
iF=100000000;  % final iteration, whichever occurs first
iRST=20000;    % write restart dump
i_report=100;
en_print=500;
trans_print=20000;
TSCREEN=10000;    % sreen update interval time (NOTE: plotting is usually slow)
initial_condition='random';   %'simple vortices' 'vortices' 'random' or 'random w' 
AB_order=-3; % Adams Bashforth order 1,2,3, or 4 (3 more accurate, 2 possibly more stable) -1 = RK3
linear_term='exact';  % CN, BE, FE, or exact
simulation_type='NL'; % NL, QL, or L (nonlinear, quasilinear and linear respectively)
padding = true;       % 3/2 padding,spspec otherwise 2/3 truncation (latter doesn't work yet)
with_plotting = true; % save plots to file
save_plots = true;   % save plots to file
system_type='MHM'; % NS, HM, MHM
holland=false;      % Holland-type model (remove adiabateic ExB nonlinearity)
waltz  =false;       % Waltz idelta (remove nonadiabtic ExB nonlinearity after linear terms calculated)
cfl_cadence=1;
cfl=0.5;
max_dt=5e-1; 
safety=0.8;
diagnostics=false;
rng(707296708);
%rng('shuffle');
s=rng;

if(nargin > 0)
  HM = HM_in; 
end
if(nargin > 1)
  delta_0 = TH_in; 
end

% print log file.
diary([basename 'log.txt']) 
diary on;

fig1=0;
fig2=0;
if(with_plotting)
  fig1=figure(1);
  fig2=figure(2);
end

% ensure parameters get printed  to log.
fprintf('mu: ');
for i = 1:length(mu), fprintf('  %.05e',mu(i)); end; fprintf('\n');
fprintf('nu: ');
for i = 1:length(nu), fprintf('  %.05e',nu(i)); end; fprintf('\n');
fprintf('muZF: ');
for i = 1:length(muZF), fprintf('  %.05e',muZF(i)); end; fprintf('\n');
fprintf('nuZF: ');
for i = 1:length(nuZF), fprintf('  %.05e',nuZF(i)); end; fprintf('\n');
fprintf('gamma:%.05e l: %.05e  HM:%.05e TH:%.05e\n',gamma,l,HM,delta_0);
fprintf('LX:%.02f LY:%.02f NX:%d NY:%d\n',LX, LY, NX_real, NY_real);
fprintf('scale:%d Tf:%.01f iF:%d\n', scale, TF, iF);
fprintf('muZF:%e nuZF:%e neo:%e tau:%e\n',muZF,nuZF,neo,tau);
fprintf('Nonlinear:%s padding:%d System:%s LinSolv:%s\n',simulation_type, padding,system_type,linear_term);
fprintf('random seed:%d AB order:%d CFL step:%d\n',s.Seed, AB_order, cfl_cadence);
fprintf('safety:%f perturbation size:%.05e\n',safety, pert_size);

energyFile=0; 
if(strcmp(initial_condition,'restart'))
  energyFile = fopen([basename 'energy.dat'],'a');
else
  energyFile = fopen([basename 'energy.dat'],'w');
end

fprintf(energyFile,'# [1] t  [2] dt  [3] energy  [4] enstrophy  [5] ZF    [6] DW    [7] flux\n');

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
minkx_p= -(NX/2 )*dkx;
maxkx_p=  (NX/2 - 1)*dkx;
minky_p= -(NY/2 )*dky;
maxky_p=  (NY/2 - 1)*dky;
kxnum_p= minkx_p:dkx:maxkx_p;
kynum_p= minky_p:dky:maxky_p;


t=0.;
i=0;
enk=0;
rsk=0;
rek=0;
trk=0;
c1=0;
c2=0;
c3=0;
dt1=0;
dt2=0;
dt3=0;
w_hat=0; %this is our field.
I=sqrt(-1);

kx=dkx*ones(1,NY)'*(mod((1:NX)-ceil(NX/2+1),NX)-floor(NX/2)); % matrix of wavenumbers in x direction 
ky=dky*(mod((1:NY)'-ceil(NY/2+1),NY)-floor(NY/2))*ones(1,NX); % matrix of wavenumbers in y direction 

% Cutting of frequencies using the 2/3 rule. 
% Discard asymmetric N/2 term. Otherwise reality is not enforced.
dealias=abs(kx/dkx) <1/3*NX & abs(ky/dky) <1/3*NY;
%dealias=dealias & kx.*ky >= 0
dealias(1,1)=0;

ksquare=kx.^2+ky.^2;                                % Laplacian in Fourier space
ksquare(1,1)=1; 
kmu = build_friction(mu);               % Friction
knu= build_viscosity(nu);              % Viscosity

G0 = exp(-tau*ksquare).*besseli(0,tau*ksquare);
G1 = sqrt(G0);

%TH= I*delta_0*ky.*ksquare./(1.0 + ksquare); % Terry-Horton term
TH= I*delta_0*ky;

ksquare_poisson=-(ksquare + ones(NY, NX) - TH);    % Poisson equation in Fourier space
if(tau> 0)
  ksquare_poisson=-((1-G0)/tau + ones(NY, NX) - TH)./G1; 
end
if(strcmpi(system_type,'NS'))
  ksquare_poisson=-ksquare;    % Poisson equation in Fourier space
end
ksquare_poisson(1,1)=-1;  % fixed Laplacian in Fourier space for Poisson's equation

gd = (-gamma * ky.^2./ksquare_poisson + HM*I*ky./ksquare_poisson).*G1;
%gd = (-gamma * abs(ky)./ksquare_poisson + HM*I*ky./ksquare_poisson).*G1;
%gd = (gamma * abs(ky)./ksquare_poisson.^2 + HM*I*ky./ksquare_poisson).*G1;
%gd = (-gamma./ksquare_poisson.*(ky.^2./ksquare) + HM*I*ky./ksquare_poisson).*G1;


%ksquare_poisson=ksquare - ones(NY, NX);  
%ksquare_poisson(1,1)=-1;

zonal_part = zeros(NY,NX);
zonal_part(1,:) = zonal_part(1,:) + ones(1,NX);
fluct_part = ones(NY,NX) - zonal_part;

lin_growth = kmu + knu + gd - l*abs(ky); % this is the linear growth rate used in computations


% No damping on zonal modes. Modified Poisson equation (proper adiabatic electron response)
if(strcmpi(system_type,'MHM'))
  kmuZF=build_friction(muZF);
  knuZF=build_viscosity(nuZF);
  lin_growth(1,:) = zeros(1,NX) + kmuZF(1,:) + knuZF(1,:);
  
  ksquare_poisson(1,:)=-neo*ksquare(1,:);    % Poisson equation in Fourier space
  if(tau> 0)
    G0_neo = exp(-neo*tau*ksquare.*(1.0+neo1*tanh(neo2.*ksquare))).*besseli(0,neo*tau*ksquare.*(1.0+neo1*tanh(neo2.*ksquare)));
    ksquare_poisson(1,:)=-((1-G0_neo(1,:))/tau )./G1(1,:); 
  end 
end

numel(find(lin_growth(:)>0))
ksquare_poisson(1,1)=-1;

lin_trunc = dealias.*lin_growth;
[max_growth,indg] = max(real(lin_trunc(:)))
sum(lin_trunc(:))
max_freq = max(imag(lin_trunc(:)))
max_rate = max(abs(lin_trunc(:)))
if(~strcmpi(linear_term,'exact'))
  max_dt = min(safety * cfl /max_rate,max_dt);
  if(dt > max_dt)
    dt = max_dt;
  end
end

if max_growth == 0
% cd('..');
% return;
end

kpo=-(ksquare + ones(NY, NX));    % No idelta exb
kpo(1,:)=-neo*ksquare(1,:); 
if(waltz)
  kpo = ksquare_poisson;
  ksquare_poisson=-(ksquare + ones(NY, NX));    % No idelta exb
  ksquare_poisson(1,:)=-neo*ksquare(1,:);    % Poisson equation in Fourier space
end

% forcing stuff?
kf_min = 16*dkx;
dkf = 0.5*dkx;
forcing_base = abs(ksquare) <= (kf_min+dkf)^2 & abs(ksquare) >= kf_min^2;
nforce = 1.0/sum(forcing_base(:));

% Define initial vorticity distribution
switch lower(initial_condition)
  case {'vortices'}
    [i,j]=meshgrid(1:NX,(1:NY));
    phi=exp(-((i*dkx-pi).^2+(j*dky-pi+pi/4).^2)/(0.2))+exp(-((i*dx-pi).^2+(j*dky-pi-pi/4).^2)/(0.2))-0.5*exp(-((i*dkx-pi-pi/4).^2+(j*dky-pi-pi/4).^2)/(0.4));
    case {'simple vortices'}
      [i,j]=meshgrid((1:NX)*(2*pi/NX),(1:NY)*(2*pi/NY));
      phi=1*sin(i).*cos(j).^2;
    case {'zero'}
      w_hat = zeros(NY,NX);
    case {'random phi'}
      
      phi=pert_size*(2*rand(NX,NY) - 1);%normally 1e-3
      w_hat=ksquare_poisson.*fft2(phi);
    case {'4mt'}
      w=pert_size*(2*rand(NY,NX)-1);%normally 5e-2
      w_hat=fft2(w);
      w_hat(2,1)=real(w_hat(2,1));
      w_hat(NY,1)=real(w_hat(NY,1));
      w_hat(1,2)=real(w_hat(1,2));
      w_hat(1,NX)=real(w_hat(1,NX));
      %dealias.*w_hat
      %make is sine
      w_hat(1,1)=0;

  %    dealias.*w_hat
 %     w_hat = dealias.*w_hat.*ksquare_poisson;
      
      %w_hat(1,:)=zeros(1,NX);
      %w_hat(1,2) =abs(14.3986/ksquare_poisson(1,2));
      %w_hat(1,NY) = w_hat(1,2);
      %w_hat=zonal_part.*w_hat;
    case {'random'}
      %w=pert_size*(2*rand(NY,NX)-1);%normally 5e-2
      w=pert_size*randn(NY,NX);%normally 5e-2
      
      w_hat=fft2(w);
      %w_hat(NX/2 + 1,:) = 0;
      w_hat(1,:)=zeros(1,NX);
      %w_hat=zonal_part.*w_hat;
    case {'waves'}
      [i,j]=meshgrid((1:NX)*(2*pi/NX),(1:NY)*(2*pi/NY));
      %w=20*sin(32*j) + 1e-1*cos(13*i);
      w=2e1*sin(32*j) + 5e-2*(2*rand(NY,NX)-1);
      w_hat=fft2(w);
      w_hat(1,:)=fluct_part.*w_hat;
    case {'restart'}
      fileID=fopen('final_state.bin','r');
      t=fread(fileID,1,'double');
      dt=fread(fileID,1,'double');
      i=fread(fileID,1,'int');
      enk=fread(fileID,1,'int');
      trk=fread(fileID,1,'int');
      rsk=fread(fileID,1,'int');
      rek=fread(fileID,1,'int');
      dt1=fread(fileID,1,'double');
      dt2=fread(fileID,1,'double');
      dt3=fread(fileID,1,'double');
      %c1=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');%fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
      %c2=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');%fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
      %c3=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');%fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');    
      c1=fread(fileID,1,'double')+I*fread(fileID,1,'double');%fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
      c2=fread(fileID,1,'double')+I*fread(fileID,1,'double');%fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
      c3=fread(fileID,1,'double')+I*fread(fileID,1,'double');%fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');    
  
      w_hat=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');            
     
      %c1 = c1.*ksquare_poisson./kpo;
      %c2 = c2.*ksquare_poisson./kpo;
      %c3 = c3.*ksquare_poisson./kpo;
      %w_hat = w_hat.*ksquare_poisson./kpo;
      fclose(fileID);        
    otherwise
      disp('Unknown initial conditions !!!');
      return
end

if(padding) % ensure there's no energy in the padded region of k-space.
  w_hat=dealias.*w_hat;
end
w_hat=enforceReality(w_hat);
w_hat(1,1)=0; % Gauge condition. Should be redundant.
w_hat0=w_hat; % Keep initial conditions for possible future diagnostics.

u=0;
v=0;
U_ZF=0;
force=zeros(NY,NX);

if(with_plotting && save_plots) %plot growth rate contours
  plotgrowth()
end


nextScreen = 0;
tic
while t<TF && i<iF
    dwen=1;
    phi_hat = w_hat./ksquare_poisson;  % Solve Poisson's Equation
    if(any(isnan(w_hat(:)))) % Something really bad has happened.
        disp(sprintf('Divergence at iteration: %d',i));
        return;
    end
    if(mod(i,i_report)==0)
       disp(sprintf('iteration: %d    dt:%.02e     t:%.03e     step time:%.1d s',i,dt,t,toc)); 
       tic;
        diary off; %flush diary
    diary on;
    end
    if (mod(i,en_print)== 0) 
        dwen=outputEnergy();
    end
    
   if(mod(i,trans_print)==0)
     %calculate_transfer_map(0,indg);
   end
    
   if (mod(i,TSCREEN)== 0 && with_plotting) 
    % if(t > nextScreen && with_plotting)
        plotfunc();
        nextScreen = t + TSCREEN;
    end
    
    if (mod(i,iRST) == 0)
      dump(sprintf('restart_%d.bin',rek));
      rek=rek+1;
    end
    

    if(forcing ~= 0)
       force=calculate_forcing(); 
    end
    
    w_hat_new=0;
    conv_hat =0;
    % Compute the non-linear term
    if (AB_order < 0 ) % RK
        if(AB_order == -1)
        w_hat_new = RK3(w_hat);
        elseif(AB_order == -4)
        w_hat_new = ETDRK4(w_hat);
        elseif(AB_order == -3)
        w_hat_new = ETDRK3(w_hat);
        elseif(AB_order == -2)
        w_hat_new = ETDRK2(w_hat);
        else
          exit();      
        end
      if(mod(i+1,cfl_cadence)==0 && i > 3) % compute new timestep from CFL condition.
          abs_u=abs(u);
          abs_v=abs(v);
          abs_ZF=abs(U_ZF);
          maxV= max(abs_u(:)+abs_ZF(:))/dx + max(abs_v(:))/dy;
          if(maxV>0)
              new_dt=1/maxV;
          else
              new_dt=inf;
            end
            target_dt=min(cfl*new_dt,min(4.0*dt/safety,max_dt/safety));
            if(target_dt < dt)
                 disp('WARNING: New dt fell below safety.')
            end
            if target_dt < dt/safety || target_dt > 3.2*dt
                dt = max_dt/safety;
                while dt > target_dt 
                    dt=dt/2.0;
                end
                dt=safety*dt;
             end
      end
    else 
        conv_hat = calc_Nonlinear(w_hat,simulation_type);
         
     
    if(mod(i+1,cfl_cadence)==0 && i > 3) % compute new timestep from CFL condition.
        abs_u=abs(u);
        abs_v=abs(v);
        abs_ZF=abs(U_ZF);
        maxV= max(abs_u(:)+abs_ZF(:))/dx + max(abs_v(:))/dy;
        if(maxV>0)
            new_dt=1/maxV;
        else
            new_dt=inf;
        end
        target_dt=min(max(0.5*dt/safety,min(cfl*new_dt,1.1*dt/safety)),max_dt/safety);
        if(target_dt < dt)
            disp('WARNING: New dt fell below safety.')
        end
        dt=safety*target_dt;
    end
    
    conv_hat = dealias.*(conv_hat + nforce*forcing*force/sqrt(dt));
    
    L1=0;
    L2=0;
    w_prev=w_hat;
    
    A=1.0; B=0; C=0; D=0;
    w_hat_new = 0;
    if (i < 1 || AB_order == 1) %Forward-Euler to generate history. Should run a small time-step at first. 
       %do nothing w_hat_new = w_hat_new - L2.*conv_hat;
    elseif (i < 2 || AB_order == 2)
      w1=dt1/dt;
      A=(1.0 + 0.5*w1);
      B=0.5*w1;
     elseif (i < 3 || AB_order == 3)
      w1=dt1/dt;
      w2=dt2/dt;
      A=(2.0 + 3.0*w2 + 6.0*w1*(1+w1+w2))/(6.0*w1*(w1+w2));
      B=(2.0 + 3.0*w1 + 3.0*w2)/(6.0*w1*w2);
      C=(2.0 + 3.0*w1)/(6.0*w2*(w1+w2));
     elseif (i < 4 || AB_order == 4)
        w1=dt1/dt;
        w2=dt2/dt;
        w3=dt3/dt;
        A=(3.0 + 8.0*w2 + 4.0*w3 + 6.0*(2*w1^3 + w2*(w2+w3) +2.0*w1*(1.0+w2)*(1.0+w2+w3) + w1^2*(3.0 + 4.0*w2+2*w3)));
        A=A/(12.0*w1*(w1+w2)*(w1+w2+w3));
        B=(3.0 + 6.0*(w1+w2)*(w1+w2+w3) + 4.0*(2.0*(w1+w2)+w3))/(12.0*w1*w2*(w2+w3));
        C=(3.0 + 6.0*w1*(w1+w2+w3)+4.0*(2.0*w1 + w2+w3))/(12.0*w2*(w1+w2)*w3);
        D=(3.0 + 6.0*w1*(w1+w2)+4.0*(2.0*w1 + w2))/(12.0*w3*(w2+w3)*(w1+w2+w3));
    end
  
    
    
    switch upper(linear_term)
        case {'CN'} %crank nicolson
            L1=(1 + 0.5*dt*lin_growth)./(1 - 0.5*dt*lin_growth);
            L2=dt./(1 - 0.5*dt*lin_growth);    
        case {'BE'} %backward euler
            L1=1.0./(1 - dt*lin_growth);
            L2=dt./(1 - dt*lin_growth);    
        case {'EXACT'}%exact
            L1=exp(dt*lin_growth);
            L2=dt*L1;   
            B=B.*exp((dt1)*lin_growth);
            C=C.*exp((dt1+dt2)*lin_growth);
            D=D.*exp((dt1+dt2+dt3)*lin_growth);
        otherwise
            disp('Unknown linear handling type.');
            return;
    end
    
    % Compute Solution at the next step
  % Implicitly solve the linear term with 2nd order Crank-Nicholson
    w_hat_new = L1.*w_prev;
    w_hat_new = w_hat_new - L2.*(A.*conv_hat - B.*c1 + C.*c2 - D.*c3);
  
    end
    t=t+dt;
    i=i+1;

    w_hat=w_hat_new;
    
    w_hat=enforceReality(w_hat);
    w_hat(1,1)=0;
    
    c3=c2;
    c2=c1;
    c1=conv_hat;
    
    dt3=dt2;
    dt2=dt1;
    dt1=dt;
  %  if dwen < 1e-15
  %      break
  %  end
end
fclose(energyFile);

%write final state 
fprintf('Simulation ended at time %.03e and iteration %d\n',t,i);
dump('final_state.bin');

if(strcmpi(simulation_type,'L') && diagnostics)
   w0f=w_hat0.*exp(t*lin_growth); 
   comp=abs(w_hat-w0f)./abs(w0f);
   fprintf('Maximum linear error:%e\n',max(comp(:)));
end

% SIMULATION ENDS HERE

%%----------------------------------%%
%% Helper Function Definitions      %% 
%%----------------------------------%%

    function v=build_friction(amp)
        v = zeros(NY,NX);
        for i1=1:length(amp)
           v = v - amp(i1)*ones(NY,NX) ./ ksquare.^(i1-1); 
        end
    end

    function v=build_viscosity(amp)
        v = zeros(NY,NX);
        for i1=1:length(amp)
           v = v - amp(i1)*ones(NY,NX) .* (rboost*kx.^2+ky.^2).^(i1); 
          %v = v - amp(i1)*ones(NY,NX) .* (ky.^2).^(i1); 
        end
    end

  function f=calculate_forcing()
    fint=zeros(NY,NX);
    %for j1 = 1:NX
    %  for i1 = 1:NY
    %    if(forcing_base(i1,j1) ~=0 )
    %      fint(i1,j1) = randn(1)+I*randn(1);
    %    end    
    %  end
    %end
    %f=enforceReality(fint)*NX*NY;
    f=enforceReality(forcing_base.*(randn(NY,NX) + I*randn(NY,NX)))*NX*NY;
  end

    function x=RK3(w_h)
        a =[8./15.,5./12.,0.75];
        b =[0,-17.0/60.0,-5.0/12.0];
        Fa =0.0;
        u_new=w_h;
        for j1 = 1:3 
            Fb = Fa;
            Fa = -calc_Nonlinear(u_new,simulation_type);
            u_new = (1 + 0.5*(a(j1)+b(j1))*dt*lin_growth)./(1 - 0.5*(a(j1)+b(j1))*dt*lin_growth).*u_new;
            u_new = u_new + dt*(a(j1)*Fa +b(j1)*Fb)./(1-0.5*(a(j1)+b(j1))*dt*lin_growth);
        end
        x=u_new;
    end   

    function x=ETDRK4(w_h)
        u_new=w_h;
        persistent dto lg M r Q1 f1 f2 f3 isreal;
        if(isempty(dto))
           dto=-1;
           lg=lin_growth(:);
       if(any(imag(lg(:))))
         M=128;
         r=exp(2*I*pi*((1:M)-0.5)/M);
               isreal=0;
       else
         M=64;
         r=exp(I*pi*((1:M)-0.5)/M);
               isreal=1;
       end
        end
        
        st=simulation_type;
        eh = exp(0.5*lin_growth*dt);
        e1 = exp(lin_growth*dt);

        if(dto ~= dt)
            LR=dt*lg(:,ones(M,1)) + r(ones(NY*NX,1),:);
            elr=exp(LR);
            if(isreal)
               Q1 =dt*real(mean(       (exp(LR/2)-1)./LR              ,2));
               f1 =dt*real(mean(     (-4   -LR        +elr.*(4-3*LR+LR.^2))./LR.^3  ,2));
               f2 =dt*real(mean(    ( 2  +LR        +elr.*(-2+LR))./LR.^3      ,2));
               f3 =dt*real(mean(     (-4 -3*LR -LR.^2 +elr.*(4-LR))./LR.^3      ,2));
            else
               Q1 =dt*(mean(       (exp(LR/2)-1)./LR              ,2));
               f1 =dt*(mean(     (-4   -LR        +elr.*(4-3*LR+LR.^2))./LR.^3  ,2));
               f2 =dt*(mean(    ( 2  +LR        +elr.*(-2+LR))./LR.^3      ,2));
               f3 =dt*(mean(     (-4 -3*LR -LR.^2 +elr.*(4-LR))./LR.^3      ,2));
            end
            Q1=reshape(Q1,NY,NX);
            f1=reshape(f1,NY,NX);
            f2=reshape(f2,NY,NX);
            f3=reshape(f3,NY,NX);
            dto=dt;
        end
        A1 = -calc_Nonlinear(u_new,st);
        
        an = u_new.*eh + Q1.*A1;

        A2 = -calc_Nonlinear(an,st);
        
        bn = u_new.*eh + Q1.*A2;
      
        A3 = -calc_Nonlinear(bn,st);
        
        cn = an.*eh + Q1.*(2.0*A3-A1);
        
        A4 = -calc_Nonlinear(cn,st);
        
        dn = u_new.*e1 + (A1.*f1 + 2*(A2+A3).*f2 + A4.*f3); 
        
        x = dealias.*dn;   

    end   

    function x=ETDRK3(w_h)
        u_new=w_h;
        persistent dto lg M r Q1 Q2 f1 f2 f3 isreal;
        if(isempty(dto))
           dto=-1;
           lg=lin_growth(:);
       if(any(imag(lg(:))))
         M=128;
         r=exp(2*I*pi*((1:M)-0.5)/M);
               isreal=0;
       else
         M=64;
         r=exp(I*pi*((1:M)-0.5)/M);
               isreal=1;
       end
        end
        
        st=simulation_type;
        eh = exp(0.5*lin_growth*dt);
        e1 = exp(lin_growth*dt);

        if(dto ~= dt)
            LR=dt*lg(:,ones(M,1)) + r(ones(NY*NX,1),:);
            elr=exp(LR);
           if(isreal)
            Q1 =dt*real(mean(       (exp(LR/2)-1)./LR              ,2));
            Q2 =dt*real(mean(       (elr-1)./LR                                ,2));
            f1 =dt*real(mean(     (-4   -LR        +elr.*(4-3*LR+LR.^2))./LR.^3  ,2));
            f2 =dt*real(mean(    4*( 2  +LR        +elr.*(-2+LR))./LR.^3      ,2));
            f3 =dt*real(mean(     (-4 -3*LR -LR.^2 +elr.*(4-LR))./LR.^3      ,2));
           else
             Q1 =dt*(mean(       (exp(LR/2)-1)./LR              ,2));
            Q2 =dt*(mean(       (elr-1)./LR                                ,2));
            f1 =dt*(mean(     (-4   -LR        +elr.*(4-3*LR+LR.^2))./LR.^3  ,2));
            f2 =dt*(mean(    4*( 2  +LR        +elr.*(-2+LR))./LR.^3      ,2));
            f3 =dt*(mean(     (-4 -3*LR -LR.^2 +elr.*(4-LR))./LR.^3      ,2));     
           end
            Q1=reshape(Q1,NY,NX);
            Q2=reshape(Q2,NY,NX);
            f1=reshape(f1,NY,NX);
            f2=reshape(f2,NY,NX);
            f3=reshape(f3,NY,NX);
            dto=dt;
        end
        A1 = -calc_Nonlinear(u_new,st);
        
        an =     u_new.*eh + Q1.*A1;
      
        A2 = -calc_Nonlinear(an,st);
        
        bn = u_new.*e1 + Q2.*(2.0*A2-A1);
        
        A3 = -calc_Nonlinear(bn,st);
        
        cn = u_new.*e1 + (A1.*f1 + A2.*f2 + A3.*f3); 
        
        x = dealias.*cn;   

    end   

    function x=ETDRK2(w_h)
        u_new=w_h;
        persistent dto lg M r Q1 f1 isreal;
        if(isempty(dto))
           dto=-1;
           lg=lin_growth(:);
       if(any(imag(lg(:))))
         M=128;
         r=exp(2*I*pi*((1:M)-0.5)/M);
               isreal=0;
       else
         M=64;
         r=exp(I*pi*((1:M)-0.5)/M);
               isreal=1;
       end
        end
        
        st=simulation_type;
        e1 = exp(lin_growth*dt);

        if(dto ~= dt)
            LR=dt*lg(:,ones(M,1)) + r(ones(NY*NX,1),:);
            elr=exp(LR);
            if(isreal)
              Q1 =dt*real(mean(       (elr-1)./LR           ,2));
              f1 =dt*real(mean(     (elr - 1 - LR)./LR.^2 ,2));
            else
              Q1 =dt*(mean(       (elr-1)./LR           ,2));
              f1 =dt*(mean(     (elr - 1 - LR)./LR.^2 ,2));    
            end
            Q1=reshape(Q1,NY,NX);
            f1=reshape(f1,NY,NX);
            dto=dt;
        end
        
        A1 = -calc_Nonlinear(u_new,st);
        
        an =     u_new.*e1 + Q1.*A1;
      
        A2 = -calc_Nonlinear(an,st);
        
        bn = an + f1.*(A2-A1);
        
        x = dealias.*bn;   

    end   

    function y=enforceReality(x) % x and y are in Fourier space
        x_HT = conj(circshift(rot90(x,2),[1,1]));% Hermitian conjugate transpose
        y = 0.5*(x + x_HT);
    end
%{
    function infCheck(inver,name)
       if(any(isinf(inver(:))))
              disp(sprintf('%s infinity at iteration: %d',name, i));
       end
    end

    function nanCheck(inver,name)
       if(any(isnan(inver(:))))
              disp(sprintf('%s NaN at iteration: %d',name, i));
       end
    end
%}

  function y=calc_Nonlinear(w_h,type)
    c_hat = 0;
    phi_h = G1.*w_h./ksquare_poisson;  % Solve Poisson's Equation
    w_curr= w_h;
    if(holland)
        w_curr(1,:) =-phi_h(1,:).*(1+ksquare(1,:));
    end
    
    switch upper(type)
          case {'NL'} %full non-linear
              uhat = -I*ky.*phi_h;
              vhat =  I*kx.*phi_h;
              w_xhat = I*kx.*w_curr;
              w_yhat = I*ky.*w_curr;
            
              % dealiasing here truncates if not padded, other it has no effect
              u  =real(ifft2(dealias.*uhat));      % Compute  y derivative of stream function ==> u
              v  =real(ifft2(dealias.*vhat));      % Compute -x derivative of stream function ==> v
              w_x=real(ifft2(dealias.*w_xhat));      % Compute  x derivative of vorticity
              w_y=real(ifft2(dealias.*w_yhat));
              conv     = u.*w_x + v.*w_y;         % evaluate the convective derivative (u,v).grad(w)   
              c_hat = fft2(conv);              % go back to Fourier space
      
          case {'QL'} %quasi-linear
              %for zonal part
              uphat = -I*ky.*fluct_part.*phi_h;
              vphat =  I*kx.*fluct_part.*phi_h;
        
              %for fluctuation part
              U_hat = I*kx.*zonal_part.*phi_h;
              U_xxhat = -(kx.^2).*U_hat;
              w_xhat = I*kx.*w_curr;
              w_yhat = I*ky.*fluct_part.*w_curr;       
        
              % dealiasing here truncates if not padded, other it has no effect
              u    =real(ifft2(dealias.*uphat));      % Compute  y derivative of stream function ==> u
              v    =real(ifft2(dealias.*vphat));      % Compute -x derivative of stream function ==> v
              U_ZF =real(ifft2(dealias.*U_hat));      % Compute zonal velocity
              U_xx =real(ifft2(dealias.*U_xxhat));      % Compute zonal velocity
              w_x  =real(ifft2(dealias.*w_xhat));      % Compute  x derivative of vorticity
              w_y  =real(ifft2(dealias.*w_yhat));

              conv_fluct = U_ZF.*w_y + u.*U_xx;         % evaluate the convective derivative (u,v).grad(w)   
              conv_fluct_hat = fft2(conv_fluct);               % go back to Fourier space

              conv_zonal = u.*w_x + v.*w_y;         % evaluate the convective derivative (u,v).grad(w)   
              conv_zonal_hat = fft2(conv_zonal);

              c_hat = fluct_part.*conv_fluct_hat + zonal_part.*conv_zonal_hat;
          case {'L'}%full linear 
              %do nothing
            otherwise
              disp('Unknown simulation type.');
              return 
    end
    if(padding)
        c_hat = dealias.*c_hat;
    end
    y=c_hat;
  end

  function plotgrowth()
      subplot(1,1,1);
      plotg=0;
      f=randn(NY,NX).*forcing_base;
      if(padding)
           plotg=circshift(dealias.*lin_growth,[NY_real/2,NX_real/2]); 
           plotg=plotg(2:NY_real,2:NX_real);
      else
           plotg=circshift(dealias.*lin_growth,[NY_real/2,NX_real/2]); 
           plotg=plotg(2:NY_real,2:NX_real);
      end
      %force=real(ifft2(f));
      imagesc(kxnum,kynum,real(plotg)), axis equal tight, colorbar
      %imagesc(kxnum,kynum,f), axis equal tight, colorbar
      set(gca,'Ydir','Normal')
      colormap(c_map)
      title('growth rates');    
      xlabel('kx');
      ylabel('ky');
      if(save_plots)
          saveas(gcf,[basename 'growth.png']);
     end
     drawnow
    end

    function calculate_transfer_map(i0,j0)
      phi_curr=w_hat./ksquare_poisson;
      
      i0c = i0;
      j0c = j0;

      if(i0c < 0);  i0c  = i0c + NX ; end;
      if(j0c < 0);  j0c  = j0c + NY ; end;
      
      
      maxktr=4;
      nx=ceil(maxktr/dkx);
      ny=ceil(maxktr/dky);
      l1=kxnum_p < maxktr & kxnum_p > -maxktr;
      l2=kynum_p < maxktr & kynum_p > -maxktr;

      transfer_phi = zeros(NY,NX);
      transfer_en = zeros(NY,NX);
      transfer_ens = zeros(NY,NX);
%  zonal=true;
      %norm= 1.0/(NX*NY*sqrt(mean(abs(phi_curr(j0c+1,:)).^4)));
      norm= 1.0/(NX*NY*sum(abs(phi_curr(j0c+1,:)).^2));
      
      for il = -nx:nx
        for i0o = -nx:nx

          mxc = il;    % matrix location
          myc = i0o;

          if(mxc < 0);  mxc  = mxc + NX ; end;
          if(myc < 0);  myc  = myc + NX ; end;

          ip = il+i0o; %sideband x
          jp = j0; %sideband y

          iz = (i0+i0o)-(il+i0o);          
          jz = 0;

          izc = iz;
          jzc = jz;
          ipc = ip;
          jpc = jp;

          i0c = i0 + i0o;
          j0c = j0;

          if(i0c < 0);  i0c  = i0c + NX ; end;
          if(j0c < 0);  j0c  = j0c + NX ; end;

          if(izc < 0); izc = izc + NX ; end;
          if(jzc < 0); jzc = jzc + NX ; end;
          if(ipc < 0); ipc = ipc + NX ; end;
          if(jpc < 0); jpc = jpc + NX ; end;            

          ww=ksquare_poisson;
          if(holland)              
            ww(1,:)=-(1+ksquare(1,:));
          end

          if(abs(ip) < NX/2 && abs(jp) < NY/2)
            transfer_phi(myc+1,mxc+1) = 0*transfer_phi(myc+1,mxc+1) + ...
               real(dkx*dky*(iz*jp - jz*ip)*(ww(jpc+1,ipc+1)-ww(jzc+1,izc+1)) * ...
              phi_curr(jpc+1,ipc+1)*phi_curr(jzc+1,izc+1)*conj(phi_curr(j0c+1,i0c+1))/ww(j0c+1,i0c+1))*norm;

            transfer_en(myc+1,mxc+1) = 0*transfer_en(myc+1,mxc+1) + ...
              real(dkx*dky*(iz*jp - jz*ip)*(ww(jpc+1,ipc+1)-ww(jzc+1,izc+1)) * ...
              phi_curr(jpc+1,ipc+1)*phi_curr(jzc+1,izc+1)*conj(phi_curr(j0c+1,i0c+1)))*norm;

            transfer_ens(myc+1,mxc+1)= 0*transfer_ens(myc+1,mxc+1) + ...
              real(dkx*dky*(iz*jp - jz*ip)*(ww(jpc+1,ipc+1)-ww(jzc+1,izc+1)) * ...
              phi_curr(jpc+1,ipc+1)*phi_curr(jzc+1,izc+1)*conj(phi_curr(j0c+1,i0c+1))*conj(ww(j0c+1,i0c+1)))*norm;
          end
        end
      end
      
      %fs_phi =fftshift(transfer_phi./abs(phi_curr(j0+1,i0+1)).^2)/(NX*NY);
      %fs_phi =fftshift(transfer_phi./mean(abs(fluct_part(:).*phi_curr(:).*w_hat(:))));%/(NX*NY);
      fs_phi1 = fftshift(transfer_phi);%./sqrt(mean(transfer_phi(:).^2));
      fs_en1  = fftshift(transfer_en) ;%./sqrt(mean(transfer_en(:).^2));
      fs_ens1 = fftshift(transfer_ens);%./sqrt(mean(transfer_ens(:).^2));

      save_binary_matrix(sprintf('%splots/transfer_zf_phi_%d.bin',basename,trk),kxnum_p(l1),kynum_p(l2),fs_phi1(l2,l1));
      save_binary_matrix(sprintf('%splots/transfer_zf_en_%d.bin',basename,trk), kxnum_p(l1),kynum_p(l2),fs_en1(l2,l1));     
      save_binary_matrix(sprintf('%splots/transfer_zf_ens_%d.bin',basename,trk),kxnum_p(l1),kynum_p(l2),fs_ens1(l2,l1));
     
      %norm= 1.0/(NX*NY*abs(phi_curr(j0c+1,i0c+1))^2);
      norm= 1.0/(NX*NY*sum(abs(phi_curr(j0c+1,:)).^2));
      for il = -nx:nx
        for jl = -ny:ny
          ip = i0-il;
          jp = j0-jl;

          ilc = il;
          jlc = jl;
          ipc = ip;
          jpc = jp;

          if(ilc < 0); ilc  = ilc + NX ; end;
          if(jlc < 0); jlc  = jlc + NY ; end;
          if(ipc < 0); ipc = ipc + NX ; end;
          if(jpc < 0); jpc = jpc + NY ; end;

          ww=ksquare_poisson;
          if(holland)
            ww(1,:)=-(1+ksquare(1,:));
          end

          if(abs(ip) < NX/2 && abs(jp) < NY/2)
            transfer_phi(jlc+1,ilc+1) = 0*transfer_phi(jlc+1,ilc+1) + ...
              real(dkx*dky*(il*jp - jl*ip)*(ww(jpc+1,ipc+1)-ww(jlc+1,ilc+1)) * ...
              phi_curr(jpc + 1,ipc + 1)*phi_curr(jlc+1, ilc+1)*conj(phi_curr(j0+1,i0+1))/ww(j0+1,i0+1))*norm;
            
            transfer_en(jlc+1,ilc+1) =  0*transfer_en(jlc+1,ilc+1) + ...
              real(dkx*dky*(il*jp - jl*ip)*(ww(jpc+1,ipc+1)-ww(jlc+1,ilc+1)) * ...
              phi_curr(jpc + 1,ipc + 1)*phi_curr(jlc+1, ilc+1)*conj(phi_curr(j0+1,i0+1)))*norm;
            
            transfer_ens(jlc+1,ilc+1) = 0*transfer_ens(jlc+1,ilc+1) + ...
              real(dkx*dky*(il*jp - jl*ip)*(ww(jpc+1,ipc+1)-ww(jlc+1,ilc+1)) * ...
              phi_curr(jpc + 1,ipc + 1)*phi_curr(jlc+1, ilc+1)*conj(phi_curr(j0+1,i0+1))*conj(ww(j0+1,i0+1)))*norm;
          end
        end
      end      
 
      fs_phi1 = fftshift(transfer_phi);%./sqrt(mean(transfer_phi(:).^2));
      fs_en1  = fftshift(transfer_en); %./sqrt(mean(transfer_en(:).^2));
      fs_ens1 = fftshift(transfer_ens);%./sqrt(mean(transfer_ens(:).^2));

      save_binary_matrix(sprintf('%splots/transfer_phi_%d.bin',basename,trk),kxnum_p(l1),kynum_p(l2),fs_phi1(l2,l1));
      save_binary_matrix(sprintf('%splots/transfer_en_%d.bin',basename,trk),kxnum_p(l1),kynum_p(l2),fs_en1(l2,l1));     
      save_binary_matrix(sprintf('%splots/transfer_ens_%d.bin',basename,trk),kxnum_p(l1),kynum_p(l2),fs_ens1(l2,l1));
     
      trk=trk+1;
      
      if(false)
        transfer_phi = transfer_phi/trans_ave;
        transfer_en  = transfer_en/trans_ave;
        transfer_ens = transfer_ens/trans_ave;

        rrr=2;
        m_tt = max(max(abs(transfer_phi(:))),10^(-100))/(NX*NY)^2;
        m_tt1 = max(max(abs(transfer_en(:))),10^(-100))/(NX*NY)^2;
        m_tt2 = max(max(abs(transfer_ens(:))),10^(-100))/(NX*NY)^2;
        l1=kxnum_p < rrr & kxnum_p > -rrr;
        l2=kynum_p < rrr & kynum_p > -rrr;
        set(0,'CurrentFigure',fig2);
        subplot(1,1,1)
        fs_phi = fftshift(transfer_phi)/(NX*NY)^2;
        fs_en  = fftshift(transfer_en) /(NX*NY)^2;
        fs_ens = fftshift(transfer_ens)/(NX*NY)^2;
        fs_phi1 = zeros(NY,NX);
        fs_en1  = zeros(NY,NX);
        fs_ens1 = zeros(NY,NX);
        for jj=1:NY
          for ii=1:NX
            ic=min(max(1,ii-(NY/2-jj+1)),NX);
            %ic=ii;
            %ic=min(max(1,ii-10),NY);
            fs_phi1(jj,ii) = fs_phi(jj,ic);
            fs_en1(jj,ii)  = fs_en(jj,ic);
            fs_ens1(jj,ii) = fs_ens(jj,ic);
          end
        end
      %mm=max(log(abs(fftt(:))));
      %imagesc(kxnum_p(l1),kynum_p(l2),max(log(abs(fftt(l2,l1))),mm - 2)), axis equal tight, colorbar

        subplot(1,3,1)
        imagesc(kxnum_p(l1),kynum_p(l2),real(fs_phi1(l2,l1)),[-m_tt m_tt]), axis equal tight, colorbar
        title(sprintf('spectral transfer of p = (%.1f,%.1f) at t = %.02f',dkx*i0,dky*j0,t));
        set(fig2.CurrentAxes,'Ydir','Normal')
        colormap(c_maprb)
        xlabel('qx')
        ylabel('py');
        subplot(1,3,2)
        imagesc(kxnum_p(l1),kynum_p(l2),real(fs_en1(l2,l1)),[-m_tt1 m_tt1]), axis equal tight, colorbar
        title(sprintf('spectral transfer of p = (%.1f,%.1f) at t = %.02f',dkx*i0,dky*j0,t));
        set(fig2.CurrentAxes,'Ydir','Normal')
        colormap(c_maprb)
        xlabel('qx')
        ylabel('py');
        subplot(1,3,3)
        imagesc(kxnum_p(l1),kynum_p(l2),real(fs_ens1(l2,l1)),[-m_tt2 m_tt2]), axis equal tight, colorbar
        title(sprintf('spectral transfer of p = (%.1f,%.1f) at t = %.02f',dkx*i0,dky*j0,t));
        set(fig2.CurrentAxes,'Ydir','Normal')
        colormap(c_maprb)
        xlabel('qx')
        ylabel('py');
        if(zonal)
          ylabel('px');
        end
        if(save_plots)  
          %saveas(gcf,sprintf('plots/fig_%d.ps',k),'psc');
  %         fp = fopen(sprintf('plots/ascii_transfer_%d.dat',enk),'w');
  %         for jj=1:NY
  %           for ii=1:NX
  %             fprintf(fp,'%e %e %e %e %e\n',kxnum_p(ii),kynum_p(jj),...
  %                        real(fs_phi1(jj,ii)),real(fs_en1(jj,ii)),real(fs_ens1(jj,ii)));
  %           end
  %           fprintf(fp,'\n');
  %         end
  %         fclose(fp);
          %saveas(fig2,sprintf('plots/transfer_%d.png',enk));
          figs=fig2;
          figs.PaperUnits = 'inches';
          figs.PaperPosition=[0 0 16 8];
          print(fig2,sprintf('%splots/transfer_%d.png',basename,enk),'-dpng','-r0');
        end
        drawnow
      end
    end

    function y=outputEnergy()
        diary off; %flush diary
    diary on;
        
          
        w_curr=w_hat;
        phi_curr=phi_hat;
        phi_y = I*phi_curr.*ky;
        
        enstrophy = 0.5*w_curr.*conj(w_curr)/(NX*NY)^2;
        energy = 0.5*real(-conj(phi_curr).*w_curr)/(NX*NY)^2;
        flux = 0.5*real(conj(phi_y).*w_curr)/(NX*NY)^2;
             
        enstrophy_tot = sum(enstrophy(:));
        energy_tot    = sum(energy(:)); 
        ZF_energy = sum(zonal_part(:).*energy(:));
        DW_energy = sum(fluct_part(:).*energy(:));
        
        
        flux_tot = sum(flux(:));
        
        %for zonal part
        uphat = -ky.*fluct_part.*phi_curr;
        vphat =  kx.*fluct_part.*phi_curr;
        
        %for fluctuation part
        U_hat = I*kx.*zonal_part.*phi_curr;
        Up_hat = I*kx.^2.*zonal_part.*phi_curr;
 
        u    =real(ifft2(dealias.*uphat));      % Compute  y derivative of stream function ==> u
        v    =real(ifft2(dealias.*vphat));      % Compute -x derivative of stream function ==> v
            
      
        flux1 = Up_hat.*u.*v;
        flux2 = TH.*phi_curr.*U_hat.*u;
        
        flux1_tot = mean(flux1(:));
        flux2_tot = mean(flux2(:));
        fprintf(energyFile,'%e %e %.15e %e %e %e %e %e %e \n',t,dt,energy_tot,enstrophy_tot,ZF_energy,DW_energy,flux_tot,flux1_tot, flux2_tot);
        
        iso_spectrum(energy,enstrophy);
        y= DW_energy;

    end
    function plotfunc()
              % Go back in real space omega in real space for plotting
      %calculate_transfer_map(0,indg,1)
      diary off; %flush diary
      diary on;
      w_curr=w_hat;
      phi_curr=phi_hat;
      phi=real(ifft2(phi_curr));
      w=real(ifft2(w_curr)); 

      enstrophy = 0.5*w_curr.*conj(w_curr)/(NX*NY)^2;
      energy = 0.5*real(-conj(phi_curr).*w_curr)/(NX*NY)^2;                  


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
      m_phi = max(abs(phi(:)));
      m_w = max(abs(w(:)));
      set(0,'CurrentFigure',fig1);
      clf(fig1,'reset')
      cf=subplot(2,2,1);
      imagesc(LXnum,LYnum,phi, [-m_phi m_phi]), axis equal tight, colorbar
      set(fig1.CurrentAxes,'Ydir','Normal')
      set(fig1.CurrentAxes,'Xdir','Normal')
      colormap(cf,c_maprb)
      title(sprintf('potential t=%.02f',t));
      xlabel('x');
      ylabel('y');
      cf=subplot(2,2,2);
      imagesc(LXnum,LYnum,w,[-m_w m_w]), axis equal tight, colorbar
      colormap(cf,c_maprb)
      title(sprintf('vorticity t=%.02f',t));
      set(fig1.CurrentAxes,'Ydir','Normal')
      set(fig1.CurrentAxes,'Xdir','Normal')
      xlabel('x');
      ylabel('y');
      cf=subplot(2,2,3);
      imagesc(kxnum,kynum, energylog), axis equal tight, colorbar
      colormap(cf,c_map)
      title('log10(Energy power spectrum)');    
      set(fig1.CurrentAxes,'Ydir','Normal')
      xlabel('kx');
      ylabel('ky');
      cf=subplot(2,2,4);
      imagesc(kxnum,kynum,wlog), axis equal tight, colorbar
      colormap(cf,c_map)
      set(fig1.CurrentAxes,'Ydir','Normal')
      xlabel('kx');
      ylabel('ky');
      title('log10(vorticity/Enstrophy power spectrum)');    
      if(save_plots)
        save_binary_matrix(sprintf('%splots/phi_%d.bin',basename,enk),((1:NX)-0.5)*dx,((1:NY)-0.5)*dy,phi);
        save_binary_matrix(sprintf('%splots/w_%d.bin',basename,enk),((1:NX)-0.5)*dx,((1:NY)-0.5)*dy,w);
        %saveas(gcf,sprintf('plots/fig_%d.ps',k),'psc');
        saveas(fig1,sprintf('%splots/fig_%d.png',basename,enk));
%         fp = fopen(sprintf('plots/ascii_%d.dat',enk),'w');
%         for ii=1:NX
%           for jj=1:NY
%             fprintf(fp,'%e %e %e\n',dx*ii,dy*jj,w(jj,ii));
%           end
%           fprintf(fp,'\n');
%         end
%         fclose(fp);
      end
      drawnow
      enk=enk+1;
    end

    function iso_spectrum(energy,enstrophy)
        dw_energy= fluct_part.*energy;
        dw_enstrophy = fluct_part.*enstrophy;
        raden_arr=cat(2,sqrt(abs(ksquare(:))),energy(:),enstrophy(:),dw_energy(:),dw_enstrophy(:));
        raden_arr=sortrows(raden_arr);
        radenergy=zeros(max(NX,NY),1);
        radenstrophy=zeros(max(NX,NY),1);
        raddwenergy=zeros(max(NX,NY),1);
        raddwenstrophy=zeros(max(NX,NY),1);
        ik_old=1;
        nk=0;
        counts=zeros(max(NX,NY),1);
        for n = 1:length(ksquare(:))
          ik=round(raden_arr(n,1)/min(dkx,dky)) + 1;
          radenergy(ik) = radenergy(ik) + raden_arr(n,2);
          radenstrophy(ik) = radenstrophy(ik) + raden_arr(n,3);
          raddwenergy(ik) = raddwenergy(ik) + raden_arr(n,4);
          raddwenstrophy(ik) = raddwenstrophy(ik) + raden_arr(n,5);   
          if(ik~=ik_old)
            counts(ik_old) = nk;
            nk=0;
            ik_old=ik;
          end
          nk=nk+1;
        end    
        
        range =0:min(dkx,dky):sqrt(2)*max(maxkx,maxky);
        l=length(range);
        
        surface = 2*pi/min(dkx,dky);
        comp=surface*range./(counts(1:l).');
        
        kxpos=0:min(dkx,dky):max(maxkx,maxky);
        specFile = fopen(sprintf('%sdata/rad_spec_%04d.dat',basename,rsk),'w');
        
        fprintf(specFile,'# Spectra at time t = %f',t);
        fprintf(specFile,'# [1] k  [2] energy  [3] enstrophy [4] dw en [5] dw ens [6] ZF en [7] ZF ens [8] kx=0 en [9] ky=0 ens\n');
        radenergy_comp = radenergy(1:l).*comp.';
        radenstrophy_comp = radenstrophy(1:l).*comp.';
        raddwenergy_comp = raddwenergy(1:l).*comp.';
        raddwenstrophy_comp = raddwenstrophy(1:l).*comp.';
        
        
        ZF_en = zeros(l);
        ZF_ens = zeros(l);
        ZF_en(1:(NX/2)) = energy(1,1:(NX/2));     
        ZF_ens(1:(NX/2)) = enstrophy(1,1:(NX/2));

        DW_en = zeros(l);
        DW_ens = zeros(l);
        DW_en(1:(NY/2)) = energy(1:(NY/2),1);     
        DW_ens(1:(NY/2)) = enstrophy(1:(NY/2),1);
        [~,maxQ] = max(ZF_en(:));
        %range(maxQ)

        
        for n = 1:l
            fprintf(specFile,'%e %e %e %e %e %e %e %e %e\n',range(n),radenergy_comp(n),radenstrophy_comp(n), ...
                                                    raddwenergy_comp(n),raddwenstrophy_comp(n),   ...                                             
                                                     ZF_en(n),ZF_ens(n), DW_en(n),DW_ens(n));       
        end
        
        fclose(specFile); 
        rsk=rsk+1;
    end


    function dump(filename)
        fileID=fopen(filename,'w');
        fwrite(fileID,t,'double');
        fwrite(fileID,dt,'double');
        fwrite(fileID,i,'int');
        fwrite(fileID,enk,'int');
        fwrite(fileID,trk,'int');
        fwrite(fileID,rsk,'int');
        fwrite(fileID,rek,'int');
        fwrite(fileID,dt1,'double');
        fwrite(fileID,dt2,'double');
        fwrite(fileID,dt3,'double');
        fwrite(fileID,real(w_hat),'double');fwrite(fileID,imag(w_hat),'double');
        fwrite(fileID,real(c1),'double'); fwrite(fileID,imag(c1),'double');
        fwrite(fileID,real(c2),'double'); fwrite(fileID,imag(c2),'double');
        fwrite(fileID,real(c3),'double'); fwrite(fileID,imag(c3),'double');
        fclose(fileID);
    end
end
