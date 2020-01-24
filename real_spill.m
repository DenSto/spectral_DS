function real_spill(HM_in,TH_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-field two-scale model of the Terry-Horton/Hasegawa-Mima equation            %
%                                                                               %
% Based off an MIT code originally made by Jean-Christophe Nave.                %
% Modified by Denis St-Onge                                                     %
%                                                                               %
% Should be OK!                                                                 %
%                                                                               %
% Laplacian(phi) = w                                                            %
% u = phi_y                                                                     %
% v =-phi_x                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear dto lg M r Q1 Q2 f1 f2 f3 isreal;
clearvars -except HM_in TH_in;
basename='spreading_big';
if(nargin > 0)
  basename=['TH-HW-',num2str(HM_in)]; %basename_in;
end
if(nargin > 1)
  basename=['d-',num2str(TH_in),'-h-',num2str(HM_in)]; %basename_in;
end
cm_redblue =flip(cbrewer('div', 'RdBu',129));
cm_inferno=inferno();
c_map=cm_inferno;
c_maprb=cm_redblue;
%if(7==exist(basename,'dir'))
%  return;
%end
mkdir(basename);cd(basename);
mkdir('plots'); mkdir('data');
delete('log.txt');


scale=1;                 % quick scale of the linear terms
mu_l     = scale*[1];         % friction
nu_l     = scale*[0.01,0,1e-5];      %viscosity
hm_l     = scale*6.5;            % HM-type wave
delta0_l = scale*1.5;        % Terry-Horton i delta 
l_scale  = 1;


mu_r     = mu_l;         % friction
nu_r     = nu_l;      %viscosity
hm_r     = scale*3.5;            % HM-type wave
delta0_r = scale*1.5;        % Terry-Horton i delta
r_scale  = 1;

mu_c     = mu_l;         % friction
nu_c     = nu_l;      %viscosity
hm_c     = 0.5*(hm_l + hm_r);    % HM-type wave
delta0_c = scale*1.5;        % Terry-Horton i delta
hm_sec   = 0;

shear=0;

forcing=0.0;              % forcing magnitude

LX   = 2*pi*10;         % X scale
LY   = 2*pi*10;         % Y scale

NX_real=256;        % resolution in x
NY_real=256;        % resolution in y
NXc_real=1024;     % resolution in x (center)

NXneigh=48;     % number of extra cells
LX_c = 2*pi*40;       % X scale

hm_p     = (hm_r - hm_l)/LX_c;

dt=5e-4;            % time step. Should start small as CFL updated can pick up the pace
pert_size=1e-3;     % size of perturbation
TF=2000.0;           % final time
iF=100000;  % final iteration, whichever occurs first
iRST=20000;    % write restart dump
i_report=100;
en_print=500;

t_begin=-1;
TSCREEN=500;    % sreen update interval time (NOTE: plotting is usually slow)
initial_condition='random w';   %'simple vortices' 'vortices' 'random' or 'random w'
AB_order=3; % Adams Bashforth order 1,2,3, or 4 (3 more accurate, 2 possibly more stable) -1 = RK3
linear_term='exact';  % CN, BE, FE, or exact
simulation_type='NL'; % NL, QL, or L (nonlinear, quasilinear and linear respectively)
with_plotting = true; % save plots to file
save_plots = true;   % saveplots to file
zonal=true;
cfl_cadence=1;
cfl=0.2;
max_dt=5e-2;
safety=0.5;
rng(707296708);
%rng('shuffle');
s=rng;

if(nargin > 0)
  hm = HM_in;
end
if(nargin > 1)
  delta_0 = TH_in;
end

% print log file.
diary('log.txt')
diary on;

fig1=0;
%fig2=0;
if(with_plotting)
  fig1=figure(1);
  %  fig2=figure(2);
end

% ensure parameters get printed  to log.
fprintf('mu_l: ');
for i = 1:length(mu_l), fprintf('  %.05e',mu_l(i)); end; fprintf('\n');
fprintf('nu_l: ');
for i = 1:length(nu_l), fprintf('  %.05e',nu_l(i)); end; fprintf('\n');
fprintf('HM_l:%.05e TH_l:%.05e\n',hm_l,delta0_l);
fprintf('mu_c: ');
for i = 1:length(mu_c), fprintf('  %.05e',mu_c(i)); end; fprintf('\n');
fprintf('nu_c: ');
for i = 1:length(nu_c), fprintf('  %.05e',nu_c(i)); end; fprintf('\n');
fprintf('HM_c:%.05e TH_c:%.05e\n',hm_c,delta0_c);
fprintf('mu_r: ');
for i = 1:length(mu_r), fprintf('  %.05e',mu_r(i)); end; fprintf('\n');
fprintf('nu_r: ');
for i = 1:length(nu_r), fprintf('  %.05e',nu_r(i)); end; fprintf('\n');
fprintf('HM_r:%.05e TH_r:%.05e\n',hm_r,delta0_r);

fprintf('LX:%.02f LY:%.02f NX:%d NY:%d\n',LX, LY, NX_real, NY_real);
fprintf('LXc:%.02f NXc:%d neigh:%d\n',LX_c,  NXc_real, NXneigh);
fprintf('scale:%d Tf:%.01f iF:%d\n', scale, TF, iF);
fprintf('Nonlinear:%s\n',simulation_type);
fprintf('random seed:%d AB order:%d CFL step:%d\n',s.Seed, AB_order, cfl_cadence);
fprintf('safety:%f perturbation size:%.05e\n',safety, pert_size);

energyFile=0;
if(strcmp(initial_condition,'restart'))
  energyFile = fopen('energy.dat','a');
else
  energyFile = fopen('energy.dat','w');
end

fprintf(energyFile,'# [1] t  [2] dt  [3] energy  [4] enstrophy  [5] ZF    [6] DW    [7] flux\n');

NX  = NX_real;
NY  = NY_real;
NX_c= NXc_real;
NX  = 3*NX/2;
NY  = 3*NY/2;
NX_c= 3*NX_c/2;

dx=LX/NX;
dy=LY/NY;
dkx=2*pi/LX;
dky=2*pi/LY;
LXnum=0:dx:LX;
LYnum=0:dy:LY;

dx_c  = LX_c/NX_c;
dkx_c = 2*pi/LX_c;
LXcnum = 0:dx_c:LX_c;

minkx= -(NX_real/2 - 1)*dkx;
maxkx=  (NX_real/2 - 1)*dkx;
minky= -(NY_real/2 - 1)*dky;
maxky=  (NY_real/2 - 1)*dky;
kxnum= minkx:dkx:maxkx;
kynum= minky:dky:maxky;

minkx_c = -(NXc_real/2 - 1)*dkx_c;
maxkx_c =  (NXc_real/2 - 1)*dkx_c;
kxnum_c =    minkx_c:dkx_c:maxkx_c;


%[x,y] = meshgrid((0:(NX-1))*dx,(0:(NY-1))*dy);
%x= x - LX/2;
%y= y - LY/2;
%xsq=x.^2;


t=0.;
i=0;
ic=0;
enk=0;
rsk=0;
rek=0;
cl_1=0;
cl_2=0;
cl_3=0;
cr_1=0;
cr_2=0;
cr_3=0;
cc_1=0;
cc_2=0;
cc_3=0;
dt1=0;
dt2=0;
dt3=0;
dtc=dt;

%these are our fields.
phi_l_hat=0; %short-scale
phi_r_hat=0; %long_scale
phi_c_hat=0;
I=sqrt(-1);

kx  = dkx*ones(1,NY)'*(mod((1:NX)-ceil(NX/2+1),NX)-floor(NX/2)); % matrix of wavenumbers in x direction
ky  = dky*(mod((1:NY)'-ceil(NY/2+1),NY)-floor(NY/2))*ones(1,NX); % matrix of wavenumbers in y direction
kx_c= dkx_c*ones(1,NY)'*(mod((1:NX_c)-ceil(NX_c/2+1),NX_c)-floor(NX_c/2)); % matrix of wavenumbers in x direction
ky_c= dky*(mod((1:NY)'-ceil(NY/2+1),NY)-floor(NY/2))*ones(1,NX_c); % matrix of wavenumbers in y direction

kx_cl = kx_c(1,:);
x   = (0:(NX_c-1))*dx_c; %- LX_c/2;
x_d = (0:(NXc_real-1))*LX_c/NXc_real; %- LX_c/2;
hsec = sech(6*(x-LX_c/2)/LX_c).^2;



% Cutting of frequencies using the 2/3 rule.
% Discard asymmetric N/2 term. Otherwise reality is not enforced.
dealias=abs(kx/dkx) <1/3*NX & abs(ky/dky) <1/3*NY;
dealias(1,1)=0;
dealias_c=abs(kx_c/dkx_c) <1/3*NX_c & abs(ky_c/dky) <1/3*NY;
dealias_c(1,1)=0;


dealias_t=abs(kx_c/dkx_c) <1/3*(NX_c - 3*NXneigh) & abs(ky_c/dky) <1/3*NY;
dealias_t(1,1)=0;

zonal_part   = zeros(NY, NX);
fluct_part   =  ones(NY, NX);

zonal_part_c   = zeros(NY, NX_c);
fluct_part_c   = ones(NY, NX_c);

if(zonal)
  zonal_part(1,:) = zonal_part(1,:) + ones(1,NX);
  fluct_part = fluct_part - zonal_part;
end

if(zonal)
  zonal_part_c(1,:) = zonal_part_c(1,:) + ones(1,NX_c);
  fluct_part_c = fluct_part_c - zonal_part_c;
end

ksquare = kx.^2 + ky.^2;                                % Laplacian in Fourier space
TH_l = I*delta0_l*ky;
kmu_l = build_friction(mu_l,NX,ksquare);               % Friction
knu_l = build_viscosity(nu_l,NX,ksquare);              % Viscosity
ksquare_poisson_l  =-(ksquare + fluct_part - TH_l);    % Poisson equation in Fourier space
ksquare_poisson_l(1,1)  = -1;
gd_l = hm_l*I*ky./ksquare_poisson_l;
lin_growth_l   = kmu_l   + knu_l   + gd_l; % this is the linear growth rate used in computations

ksquare = kx.^2 + ky.^2;                                % Laplacian in Fourier space
TH_r = I*delta0_r*ky;
kmu_r = build_friction(mu_r,NX,ksquare);               % Friction
knu_r = build_viscosity(nu_r,NX,ksquare);              % Viscosity
ksquare_poisson_r  =-(ksquare + fluct_part - TH_r);    % Poisson equation in Fourier space
ksquare_poisson_r(1,1)  = -1;
gd_r = hm_r*I*ky./ksquare_poisson_r;
lin_growth_r   = kmu_r   + knu_r   + gd_r; % this is the linear growth rate used in computations


ksquare_c = kx_c.^2 + ky_c.^2;                                % Laplacian in Fourier space
TH_c= I*delta0_c*ky_c;
kmu_c = build_friction(mu_c,NX_c,ksquare_c);               % Friction
knu_c = build_viscosity(nu_c,NX_c,ksquare_c);              % Viscosity
ksquare_poisson_c = -(ksquare_c + fluct_part_c - TH_c);  
ksquare_poisson_c(1,1)= -1;
gd_c = hm_c*I*ky_c./ksquare_poisson_c;
lin_growth_c = kmu_c + knu_c + gd_c; % this is the linear growth rate used in computations


ccount=NX_c - 2*NXneigh;

shear_C = zeros(1,NX_c);
shear_C(1:NXneigh) = -shear;
shear_C((NX_c-NXneigh+1):NX_c) = shear;
shear_C((NXneigh+1):(NX_c-NXneigh)) = 2*shear*((1:ccount) -ccount/2)/ccount; 


% forcing stuff
kf_min = 16*dkx;
dkf = 0.5*dkx;
forcing_base = abs(ksquare) <= (kf_min+dkf)^2 & abs(ksquare) >= kf_min^2;
nforce = 1.0/sum(forcing_base(:));

% Define initial vorticity distribution
switch lower(initial_condition)
  case {'zero'}
    phi_l_hat = zeros(NY,NX);
    phi_r_hat = zeros(NY,NX);
  case {'vortices'}
    [i,j]=meshgrid(1:NX,(1:NY));
    phi_l_hat=exp(-((i*dkx-pi).^2+(j*dky-pi+pi/4).^2)/(0.2))+exp(-((i*dx-pi).^2+(j*dky-pi-pi/4).^2)/(0.2))-0.5*exp(-((i*dkx-pi-pi/4).^2+(j*dky-pi-pi/4).^2)/(0.4));
    phi_r_hat=0;
  case {'simple vortices'}
    phi0=pert_size*sin(dkx.*x).*cos(dky.*y).^2;
    phi_l_hat=fft2(phi0);
    
    phi_r=1e-1*phi0;
    phi_r_hat=fft2(phi_r);
  case {'random'}
    phi_l = pert_size*(2*rand(NY,NX) - 1);%normally 1e-3
    phi_r = pert_size*(2*rand(NY,NX) - 1);%normally 1e-3
    phi_l_hat=fluct_part.*fft2(phi_l);
    phi_r_hat=fluct_part.*fft2(phi_r);
  case {'random w'}
    w_l = pert_size*(2*rand(NY,NX) - 1);%normally 1e-3
    w_r = pert_size*(2*rand(NY,NX) - 1);%normally 1e-3
    w_c = pert_size*(2*rand(NY,NX_c) - 1);%normally 1e-3
    %w_r = w_l;
    %w_c = w_l;
    
    phi_l_hat=fluct_part.*fft2(w_l)./ksquare_poisson_l;
    phi_r_hat=fluct_part.*fft2(w_r)./ksquare_poisson_r;
    phi_c_hat=fluct_part_c.*fft2(w_c)./ksquare_poisson_c;
  case {'restart'}
    fileID=fopen('res.bin','r');
    t=fread(fileID,1,'double');
    dt=fread(fileID,1,'double');
    i=fread(fileID,1,'int');
    enk=fread(fileID,1,'int');
    rsk=fread(fileID,1,'int');
    rek=fread(fileID,1,'int');
    dt1=fread(fileID,1,'double');
    dt2=fread(fileID,1,'double');
    dt3=fread(fileID,1,'double');
    cl_1=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');%fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    cl_2=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');%fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    cl_3=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');%fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    cr_1=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');%fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    cr_2=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');%fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    cr_3=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');%fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    cc_1=fread(fileID,[NY NX_c],'double')+I*fread(fileID,[NY NX_c],'double');%fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    cc_2=fread(fileID,[NY NX_c],'double')+I*fread(fileID,[NY NX_c],'double');%fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    cc_3=fread(fileID,[NY NX_c],'double')+I*fread(fileID,[NY NX_c],'double');%fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    phi_l_hat=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    phi_r_hat=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
    phi_c_hat=fread(fileID,[NY NX_c],'double')+I*fread(fileID,[NY NX_c],'double');
    fclose(fileID);
  otherwise
    disp('Unknown initial conditions !!!');
    return
end

phi_l_hat=dealias.*phi_l_hat;
phi_r_hat=dealias.*phi_r_hat;
phi_c_hat=dealias_c.*phi_c_hat;
phi_l_hat=enforceReality(phi_l_hat);
phi_r_hat=enforceReality(phi_r_hat);
phi_c_hat=enforceReality(phi_c_hat);

phi_l_hat(1,1) = 0;
phi_r_hat(1,1) = 0;
phi_c_hat(1,1) = 0;

u=0;
v=0;
force=zeros(NY,NX);
%if(with_plotting && save_plots) %plot growth rate contours
%  plotgrowth()
%end

[~,u_l,v_l] = calc_Nonlinear(phi_l_hat,ksquare_poisson_l.*phi_l_hat,-1);
[~,u_r,v_r] = calc_Nonlinear(phi_r_hat,ksquare_poisson_r.*phi_r_hat, 1);


nextScreen = 0;
tic
while t<TF && i<iF
  dwen=1;
  
  phi_c_hat = enforce_boundaries(phi_l_hat,phi_r_hat,phi_c_hat);
  
  if(any(isnan(phi_l_hat(:)))) % Something really bad has happened.
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
    outputEnergy();
  end
  
  
  if (mod(i,TSCREEN)== 0 && with_plotting)
    % if(t > nextScreen && with_plotting)
    plotfunc();
    %qqq=3;
    %drawnow
    nextScreen = t + TSCREEN;
  end
  
  if (mod(i,iRST) == 0)
    dump(sprintf('restart_%d.bin',rek));
    rek=rek+1;
  end
  
  
  if(forcing ~= 0)
    force=calculate_forcing();
  end
    
  if(t > t_begin)
    phi_c_hat = phi_c_hat.*dealias_c;
    
    %phi_c1=real(ifft2(phi_c_hat));
    
    phi_c_hat = enforce_boundaries(phi_l_hat,phi_r_hat,phi_c_hat);
    
    %phi_l =real(ifft2(  I.*kx.*phi_l_hat));
    %phi_c2=real(ifft2(I.*kx_c.*phi_c_hat));
    %phi_r =real(ifft2(  I.*kx.*phi_r_hat));
    
    
    [conv_c_hat,u_c,v_c] = calc_Nonlinear_C(phi_c_hat);
    
    if(mod(ic+1,cfl_cadence)==0) % compute new timestep from CFL condition.
      abs_u = abs(u_c);
      abs_v = abs(v_c);
      maxV_c = max(abs_u(:))/dx_c + max(abs_v(:))/dy;
      abs_u = abs(u_l);
      abs_v = abs(v_l);
      maxV_l = max(abs_u(:))/dx + max(abs_v(:))/dy;
      abs_u = abs(u_r);
      abs_v = abs(v_r);
      maxV_r = max(abs_u(:))/dx + max(abs_v(:))/dy;
      maxV=max(max(maxV_l,maxV_c),maxV_r);
      
      if(maxV>0)
        new_dt=1/maxV;
      else
        new_dt=inf;
      end
      target_dt=min(max(0.5*dtc/safety,min(cfl*new_dt,1.1*dtc/safety)),max_dt/safety);
      if(target_dt < dt)
        disp('WARNING: New dt fell below safety.')
      end
      dtc=safety*target_dt;
      dt=dtc;
    end
    
    AB1=1.0; AB2=0; AB3=0; AB4=0;
    w_c_prev = ksquare_poisson_c.*phi_c_hat;
    
    if (ic < 1 || AB_order == 1) %Forward-Euler to generate history. Should run a small time-step at first.
      %do nothing
    elseif (ic < 2 || AB_order == 2)
      w1=dt1/dtc;
      AB1=(1.0 + 0.5*w1);
      AB2=0.5*w1;
    elseif (ic < 3 || AB_order == 3)
      w1=dt1/dtc;
      w2=dt2/dtc;
      AB1=(2.0 + 3.0*w2 + 6.0*w1*(1+w1+w2))/(6.0*w1*(w1+w2));
      AB2=(2.0 + 3.0*w1 + 3.0*w2)/(6.0*w1*w2);
      AB3=(2.0 + 3.0*w1)/(6.0*w2*(w1+w2));
    elseif (ic < 4 || AB_order == 4)
      w1=dt1/dtc;
      w2=dt2/dtc;
      w3=dt3/dtc;
      AB1=(3.0 + 8.0*w2 + 4.0*w3 + 6.0*(2*w1^3 + w2*(w2+w3) +2.0*w1*(1.0+w2)*(1.0+w2+w3) + w1^2*(3.0 + 4.0*w2+2*w3)));
      AB1=AB1/(12.0*w1*(w1+w2)*(w1+w2+w3));
      AB2=(3.0 + 6.0*(w1+w2)*(w1+w2+w3) + 4.0*(2.0*(w1+w2)+w3))/(12.0*w1*w2*(w2+w3));
      AB3=(3.0 + 6.0*w1*(w1+w2+w3)+4.0*(2.0*w1 + w2+w3))/(12.0*w2*(w1+w2)*w3);
      AB4=(3.0 + 6.0*w1*(w1+w2)+4.0*(2.0*w1 + w2))/(12.0*w3*(w2+w3)*(w1+w2+w3));
    end
    
    % Compute Solution at the next step
    % Implicitly solve the linear term with 2nd order Crank-Nicholson
    
    switch upper(linear_term)
      case {'CN'} %crank nicolson
        L1=(1 + 0.5*dtc*lin_growth_c)./(1 - 0.5*dtc*lin_growth_c);
        L2=dtc./(1 - 0.5*dtc*lin_growth_c);
      case {'BE'} %backward euler
        L1=1.0./(1 - dtc*lin_growth_c);
        L2=dtc./(1 - dtc*lin_growth_c);
      case {'EXACT'}%exact
        L1=exp(dtc*lin_growth_c);
        L2=dtc*L1;
        AB2=AB2.*exp((dt1)*lin_growth_c);
        AB3=AB3.*exp((dt1+dt2)*lin_growth_c);
        AB4=AB4.*exp((dt1+dt2+dt3)*lin_growth_c);
      otherwise
        disp('Unknown linear handling type.');
        return;
    end
    
    % Compute Solution at the next step
    % Implicitly solve the linear term with 2nd order Crank-Nicholson
    w_c_hat_new = L1.*w_c_prev;
    w_c_hat_new = w_c_hat_new - L2.*(AB1.*conv_c_hat - AB2.*cc_1 + AB3.*cc_2 - AB4.*cc_3);
    
    ic=ic+1;
    
    phi_c_hat = w_c_hat_new./ksquare_poisson_c;
    phi_c_hat = dealias_c.*enforceReality(phi_c_hat); 
    
    %shouldn't be using these!
    %phi_c_hat(1,1) = 0;
    
    cc_3=cc_2;
    cc_2=cc_1;
    cc_1=conv_c_hat;
  end 
  
  
  conv_l_hat =0;
  conv_r_hat =0;
  % Compute the non-linear term
  if (AB_order < 0 ) % RK
    if(AB_order == -1)
      [phi0_hat_new, phil_hat_new] = RK3(phi_l_hat, phi_r_hat);
    else
      exit();
    end
    if(mod(i+1,cfl_cadence)==0 && i > 3) % compute new timestep from CFL condition.
      abs_u = abs(u_c);
      abs_v = abs(v_c);
      maxV_c = max(abs_u(:))/dx_c + max(abs_v(:))/dy;
      abs_u = abs(u_l);
      abs_v = abs(v_l);
      maxV_l = max(abs_u(:))/dx + max(abs_v(:))/dy;
      abs_u = abs(u_r);
      abs_v = abs(v_r);
      maxV_r = max(abs_u(:))/dx + max(abs_v(:))/dy;
      maxV=max(max(maxV_l,maxV_c),maxV_r);
      
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
        dt = min(max_dt/safety,dtc);
        while dt > target_dt
          dt=dt/2.0;
        end
        dt=safety*dt;
      end
    end
  else
    
    w_l_prev = ksquare_poisson_l.*phi_l_hat;
    w_r_prev = ksquare_poisson_r.*phi_r_hat;
    
    [conv_l_hat,u_l,v_l] = calc_Nonlinear(phi_l_hat,w_l_prev,-1);
    [conv_r_hat,u_r,v_r] = calc_Nonlinear(phi_r_hat,w_r_prev, 1);
      
    conv_l_hat = dealias.*(conv_l_hat + nforce*forcing*force/sqrt(dt));
    conv_r_hat = dealias.*(conv_r_hat + nforce*forcing*force/sqrt(dt));
    
    AB1=1.0; AB2=0; AB3=0; AB4=0;

    if (i < 1 || AB_order == 1) %Forward-Euler to generate history. Should run a small time-step at first.
      %do nothing
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
    
    % Compute Solution at the next step
    % Implicitly solve the linear term with 2nd order Crank-Nicholson
    
    switch upper(linear_term)
      case {'CN'} %crank nicolson
        L1=(1 + 0.5*dt*lin_growth_l)./(1 - 0.5*dt*lin_growth_l);
        L2=dt./(1 - 0.5*dt*lin_growth_l);
      case {'BE'} %backward euler
        L1=1.0./(1 - dt*lin_growth_l);
        L2=dt./(1 - dt*lin_growth_l);
      case {'EXACT'}%exact
        L1=exp(dt*lin_growth_l);
        L2=dt*L1;
        AB2=AB2.*exp((dt1)*lin_growth_l);
        AB3=AB3.*exp((dt1+dt2)*lin_growth_l);
        AB4=AB4.*exp((dt1+dt2+dt3)*lin_growth_l);
      otherwise
        disp('Unknown linear handling type.');
        return;
    end
    
    % Compute Solution at the next step
    % Implicitly solve the linear term with 2nd order Crank-Nicholson
    w_l_hat_new = L1.*w_l_prev;
    w_l_hat_new = w_l_hat_new - L2.*(AB1.*conv_l_hat - AB2.*cl_1 + AB3.*cl_2 - AB4.*cl_3);
    
    switch upper(linear_term)
      case {'CN'} %crank nicolson
        L1=(1 + 0.5*dt*lin_growth_r)./(1 - 0.5*dt*lin_growth_r);
        L2=dt./(1 - 0.5*dt*lin_growth_r);
      case {'BE'} %backward euler
        L1=1.0./(1 - dt*lin_growth_r);
        L2=dt./(1 - dt*lin_growth_r);
      case {'EXACT'}%exact
        L1=exp(dt*lin_growth_r);
        L2=dt*L1;
        AB2=AB2.*exp((dt1)*lin_growth_r);
        AB3=AB3.*exp((dt1+dt2)*lin_growth_r);
        AB4=AB4.*exp((dt1+dt2+dt3)*lin_growth_r);
      otherwise
        disp('Unknown linear handling type.');
        return;
    end
    
    % Compute Solution at the next step
    % Implicitly solve the linear term with 2nd order Crank-Nicholson
    
    w_r_hat_new = L1.*w_r_prev;
    w_r_hat_new = w_r_hat_new - L2.*(AB1.*conv_r_hat - AB2.*cr_1 + AB3.*cr_2 - AB4.*cr_3);
    
  end
  t=t+dt;
  i=i+1;
  
  phi_l_hat = w_l_hat_new./ksquare_poisson_l;
  phi_r_hat = w_r_hat_new./ksquare_poisson_r;
  
  phi_l_hat = dealias.*enforceReality(phi_l_hat);
  phi_r_hat = dealias.*enforceReality(phi_r_hat);
  
  
  %shouldn't be using these!
  %phi0_hat(1,1) = 0;
  %phil_hat(1,1) = 0;
  
  cl_3=cl_2;
  cl_2=cl_1;
  cl_1=conv_l_hat;
  
  cr_3=cr_2;
  cr_2=cr_1;
  cr_1=conv_r_hat;
  
  dt3=dt2;
  dt2=dt1;
  dt1=dt;
  
end
fclose(energyFile);

%write final state
fprintf('Simulation ended at time %.03e and iteration %d\n',t,i);
dump('final_state.bin');

cd('..');

% SIMULATION ENDS HERE

%%----------------------------------%%
%% Helper Function Definitions      %%
%%----------------------------------%%

  function ret=build_friction(amp,nnx,ks)
    ret = zeros(NY,nnx);
    for i1=1:length(amp)
      ret = ret - amp(i1)*ones(NY,nnx) ./ ks.^(i1-1);
    end
  end

  function ret=build_viscosity(amp,nnx,ks)
    ret = zeros(NY,nnx);
    for i1=1:length(amp)
      ret = ret - amp(i1)*ones(NY,nnx).*ks.^(i1);
    end
  end

  function f=calculate_forcing()
    %fint=zeros(NY,NX);
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

  function x = RK3(p0_h)
    a =[8./15.,5./12.,0.75];
    b =[0,-17.0/60.0,-5.0/12.0];
    Fa =0.0;
    u_new=p0_h;
    for j1 = 1:3
      Fb = Fa;
      Fa = -calc_Nonlinear(u0_new,u1_new);
      u_new = (1 + 0.5*(a(j1)+b(j1))*dt*lin_growth)./(1 - 0.5*(a(j1)+b(j1))*dt*lin_growth).*u_new;
      u_new = u_new + dt*(a(j1)*Fa +b(j1)*Fb)./(1-0.5*(a(j1)+b(j1))*dt*lin_growth);
    end
    x=u_new;
  end

  function y=enforceReality(x) % x and y are in Fourier space
    x_HT = conj(circshift(rot90(x,2),[1,1]));% Hermitian conjugate transpose
    y = 0.5*(x + x_HT);
    %y=fftn(real(ifftn(x)));
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

  function ret=enforce_boundaries(ph1,ph2,phc)
    
    ph1y=ifft(ph1.*ksquare_poisson_l,[],2);
    ph2y=ifft(ph2.*ksquare_poisson_r,[],2);
    ff  =ifft(phc.*ksquare_poisson_c,[],2);
     
    for ii=1:NXneigh
      ff(:,ii) = ph1y(:,NX - NXneigh +  ii);
      ff(:,NX_c - NXneigh +  ii) = ph2y(:,ii);
    end
    ret=dealias_c.*fft(ff,[],2)./ksquare_poisson_c; 
    
    
%     ph1y=ifft(ph1,[],2);
%     ph2y=ifft(ph2,[],2);
%     ff  =ifft(phc,[],2);
%      
%     for ii=1:NXneigh
%       ff(:,ii) = ph1y(:,NX - NXneigh +  ii);
%       ff(:,NX_c - NXneigh +  ii) = ph2y(:,ii);
%     end
%     ret=dealias_c.*fft(ff,[],2); 
    
  end

  function [y1,u,v] = calc_Nonlinear_C(phi_h)
    c_hat = 0; u = 0; v = 0;
    type='NL';
    switch upper(type)
      case {'NL'} %full non-linear
        uhat = -I*ky_c.*phi_h;
        vhat =  I*kx_c.*phi_h;
        w_xhat = I*kx_c.*ksquare_poisson_c.*phi_h;
        w_yhat = I*ky_c.*ksquare_poisson_c.*phi_h;
        
        
        % dealiasing here truncates if not padded, other it has no effect
        u  =real(ifft2(dealias_c.*uhat));      % Compute  y derivative of stream function ==> u
        v  =real(ifft2(dealias_c.*vhat));      % Compute -x derivative of stream function ==> v
        w_x=real(ifft2(dealias_c.*w_xhat));      % Compute  x derivative of vorticity
        w_y=real(ifft2(dealias_c.*w_yhat));
        conv     = u.*w_x + (v + shear_C).*w_y + hm_p.*u.*(x - LX_c/2) + hm_sec.*u.*hsec;         % evaluate the convective derivative (u,v).grad(w)
        c_hat = fft2(conv);              % go back to Fourier space
        
      case {'QL'} %quasi-linear
        % nope not here!
      case {'L'}%full linear
        %do nothing
      otherwise
        disp('Unknown simulation type.');
        return
    end
    y1 = dealias_t.*c_hat;
  end

  function [y1,u,v] = calc_Nonlinear(phi_h,w_hat,pm)
    c_hat = 0; u = 0; v = 0;
    type='NL';
    switch upper(type)
      case {'NL'} %full non-linear
        uhat = -I*ky.*phi_h;
        vhat =  I*kx.*phi_h;
        w_xhat = I*kx.*w_hat;
        w_yhat = I*ky.*w_hat;
        
        % dealiasing here truncates if not padded, other it has no effect
        u  =real(ifft2(dealias.*uhat));      % Compute  y derivative of stream function ==> u
        v  =real(ifft2(dealias.*vhat));      % Compute -x derivative of stream function ==> v
        w_x=real(ifft2(dealias.*w_xhat));      % Compute  x derivative of vorticity
        w_y=real(ifft2(dealias.*w_yhat));
        conv     = u.*w_x + (v  + pm*shear).*w_y;         % evaluate the convective derivative (u,v).grad(w)
        c_hat = fft2(conv);              % go back to Fourier space
        
      case {'QL'} %quasi-linear
        % nope not here!
      case {'L'}%full linear
        %do nothing
      otherwise
        disp('Unknown simulation type.');
        return
    end
    y1 = dealias.*c_hat;
  end

  function time_integrate()
    AB1=1.0; AB2=0; AB3=0; AB4=0;

    if (i < 1 || AB_order == 1) %Forward-Euler to generate history. Should run a small time-step at first.
      %do nothing
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
    
    % Compute Solution at the next step
    % Implicitly solve the linear term with 2nd order Crank-Nicholson
    
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
        AB2=AB2.*exp((dt1)*lin_growth);
        AB3=AB3.*exp((dt1+dt2)*lin_growth);
        AB4=AB4.*exp((dt1+dt2+dt3)*lin_growth);
      otherwise
        disp('Unknown linear handling type.');
        return;
    end
    
    % Compute Solution at the next step
    % Implicitly solve the linear term with 2nd order Crank-Nicholson
    w_new = L1.*w_prev;
    w_new = w_new - L2.*(AB1.*conv_hat - AB2.*cl_1 + AB3.*cl_2 - AB4.*cl_3);
  end


  function y=outputEnergy()
    diary off; %flush diary
    diary on;
    
    w_l_curr=phi_l_hat.*ksquare_poisson_l;
    w_c_curr=phi_c_hat.*ksquare_poisson_c;
    w_r_curr=phi_r_hat.*ksquare_poisson_r;
    
    energy_l = 0.5*real(-conj(phi_l_hat).*w_l_curr)/(NX*NY)^2;
    energy_c = 0.5*real(-conj(phi_c_hat).*w_c_curr)/(NX_c*NY)^2;
    energy_r = 0.5*real(-conj(phi_r_hat).*w_r_curr)/(NX*NY)^2;
    
    enstrophy_l = 0.5*w_l_curr.*conj(w_l_curr)/(NX*NY)^2;
    enstrophy_c = 0.5*w_c_curr.*conj(w_c_curr)/(NX_c*NY)^2;
    enstrophy_r = 0.5*w_r_curr.*conj(w_r_curr)/(NX*NY)^2;
    
    energy_l_tot    = sum(energy_l(:));
    energy_c_tot    = sum(energy_c(:));
    energy_r_tot    = sum(energy_r(:));
    enstrophy_l_tot = sum(enstrophy_l(:));
    enstrophy_c_tot = sum(enstrophy_c(:));
    enstrophy_r_tot = sum(enstrophy_r(:));
    
    %flux = 0.5*real(conj(phi0_y).*w_curr)/(NX*NY)^2;
    %ZF_energy = sum(zonal_part(:).*energy(:));
    %DW_energy = sum(fluct_part(:).*energy(:));
    
    %flux_tot = sum(flux(:));
    
    %for zonal part
    %     uphat = -ky.*fluct_part.*phi_curr;
    %     vphat =  kx.*fluct_part.*phi_curr;
    %
    %     %for fluctuation part
    %     U_hat = kx.*zonal_part.*phi_curr;
    %     Up_hat = kx.^2.*zonal_part.*phi_curr;
    %
    %     u    =real(ifft2(dealias.*uphat));      % Compute  y derivative of stream function ==> u
    %     v    =real(ifft2(dealias.*vphat));      % Compute -x derivative of stream function ==> v
    %
    %     flux1 = Up_hat.*u.*v;
    %     flux2 = TH.*phi_curr.*U_hat.*u;
    
    %     flux1_tot = mean(flux1(:));
    %     flux2_tot = mean(flux2(:));
    fprintf(energyFile,'%e ', t);
    fprintf(energyFile,'%e ', dt);
    fprintf(energyFile,'%e ', energy_l_tot);
    fprintf(energyFile,'%e ', energy_c_tot);
    fprintf(energyFile,'%e ', energy_r_tot);
    fprintf(energyFile,'%e ', enstrophy_l_tot);
    fprintf(energyFile,'%e ', enstrophy_c_tot);
    fprintf(energyFile,'%e ', enstrophy_r_tot);
    fprintf(energyFile,'\n');
    %iso_spectrum(energy,enstrophy);
    %y= DW_energy;
    
  end


  function plotgrowth()
    return;
  end

  function plotfunc()
    % Go back in real space omega in real space for plotting
    diary off; %flush diary
    diary on;
    wl_hat   = phi_l_hat.*ksquare_poisson_l;
    
    phil = real(ifft2(phi_l_hat));
    wl   = real(ifft2(wl_hat));
   
    phic = real(ifft2(phi_c_hat));
    wc   = real(ifft2(phi_c_hat.*ksquare_poisson_c));
    
    philrms = phil/sqrt(mean(phil(:).^2));
    phicrms = phic/sqrt(mean(phic(:).^2));
    wlrms   = wl/sqrt(mean(wl(:).^2));
    wcrms   = wc/sqrt(mean(wc(:).^2));
    
    enstrophy = 0.5*wl_hat.*conj(wl_hat)/(NX*NY)^2;
    energy = 0.5*real(-conj(phi_l_hat).*wl_hat)/(NX*NY)^2;
    
    %phic=phic(:,nslosh_ind_f);
    
    enstrophy=circshift(enstrophy,[NY_real/2,NX_real/2]);
    enstrophy=enstrophy(2:NY_real,2:NX_real);
    energy=circshift(energy,[NY_real/2,NX_real/2]);
    energy=energy(2:NY_real,2:NX_real);
    
    wlog=max(log10(enstrophy),-10);
    energylog=max(log10(energy),-10);
    m_phil = max(max(abs(philrms(:))),1e-60);
    m_phic = max(max(abs(phicrms(:))),1e-60);
    m_wl   = max(max(abs(wlrms(:))),1e-60);
    m_wc   = max(max(abs(wcrms(:))),1e-60);
    m_phil  = 3;
    m_phic = 4;
    m_wl  = 3;
    m_wc = 4;
    set(0,'CurrentFigure',fig1);
    clf(fig1,'reset')
    cf=subplot(2,3,1);
    imagesc(LXnum,LYnum,philrms, [-m_phil m_phil]), axis equal tight, colorbar
    set(fig1.CurrentAxes,'Ydir','Normal')
    set(fig1.CurrentAxes,'Xdir','Normal')
    colormap(cf,c_maprb)
    title(sprintf('potential t=%.02f',t));
    xlabel('x');
    ylabel('y');
    cf=subplot(2,3,2);
    imagesc(LXnum,LYnum,wlrms,[-m_wl m_wl]), axis equal tight, colorbar
    colormap(cf,c_maprb)
    title(sprintf('vorticity t=%.02f',t));
    set(fig1.CurrentAxes,'Ydir','Normal')
    set(fig1.CurrentAxes,'Xdir','Normal')
    xlabel('x');
    ylabel('y');
    cf=subplot(2,3,3);
    imagesc(LXcnum,LYnum,phicrms,[-m_phic m_phic]), axis equal tight, colorbar
    colormap(cf,c_maprb)
    title(sprintf('center potential t=%.02f',t));
    set(fig1.CurrentAxes,'Ydir','Normal')
    set(fig1.CurrentAxes,'Xdir','Normal')
    xlabel('x');
    ylabel('y');
    cf=subplot(2,3,4);
    imagesc(kxnum,kynum, energylog), axis equal tight, colorbar
    colormap(cf,c_map)
    title('log10(Energy power spectrum)');
    set(fig1.CurrentAxes,'Ydir','Normal')
    xlabel('kx');
    ylabel('ky');
    cf=subplot(2,3,5);
    imagesc(kxnum,kynum,wlog), axis equal tight, colorbar
    colormap(cf,c_map)
    set(fig1.CurrentAxes,'Ydir','Normal')
    xlabel('kx');
    ylabel('ky');
    title('log10(vorticity/Enstrophy power spectrum)');
    cf=subplot(2,3,6);
    imagesc(LXcnum,LYnum,wcrms, [-m_wc m_wc]), axis equal tight, colorbar
    colormap(cf,c_maprb)
    set(fig1.CurrentAxes,'Ydir','Normal')
    xlabel('x');
    ylabel('y');
    title('log10(vorticity/Enstrophy power spectrum)');
    drawnow
    if(save_plots)  
      wr_hat   = phi_r_hat.*ksquare_poisson_r;  
      phir = real(ifft2(phi_r_hat));
      wr   = real(ifft2(wr_hat));
      
      save_binary_matrix(sprintf('plots/phi_l_%d.bin',enk),((1:NX)-0.5)*dx,((1:NY)-0.5)*dy,phil);
      save_binary_matrix(sprintf('plots/w_l_%d.bin',enk)  ,((1:NX)-0.5)*dx,((1:NY)-0.5)*dy,wl);

      save_binary_matrix(sprintf('plots/phi_r_%d.bin',enk),((1:NX)-0.5)*dx,((1:NY)-0.5)*dy,phir);
      save_binary_matrix(sprintf('plots/w_r_%d.bin',enk)  ,((1:NX)-0.5)*dx,((1:NY)-0.5)*dy,wr);
      
      save_binary_matrix(sprintf('plots/phi_c_%d.bin',enk),((1:NX_c)-0.5)*dx_c,((1:NY)-0.5)*dy,phic);
      save_binary_matrix(sprintf('plots/w_c_%d.bin',enk)  ,((1:NX_c)-0.5)*dx_c,((1:NY)-0.5)*dy,wc);
      
      %saveas(gcf,sprintf('plots/fig_%d.ps',k),'psc');
      saveas(fig1,sprintf('plots/fig_%d.png',enk));
      %         fp = fopen(sprintf('plots/ascii_%d.dat',enk),'w');
      %         for ii=1:NX
      %           for jj=1:NY
      %             fprintf(fp,'%e %e %e\n',dx*ii,dy*jj,w(jj,ii));
      %           end
      %           fprintf(fp,'\n');
      %         end
      %         fclose(fp);
    end
    
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
    specFile = fopen(sprintf('data/rad_spec_%04d.dat',rsk),'w');
    
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
    fwrite(fileID,rsk,'int');
    fwrite(fileID,rek,'int');
    fwrite(fileID,dt1,'double');
    fwrite(fileID,dt2,'double');
    fwrite(fileID,dt3,'double');
    fwrite(fileID,real(cl_1),'double'); fwrite(fileID,imag(cl_1),'double');
    fwrite(fileID,real(cl_2),'double'); fwrite(fileID,imag(cl_2),'double');
    fwrite(fileID,real(cl_3),'double'); fwrite(fileID,imag(cl_3),'double');
    fwrite(fileID,real(cr_1),'double'); fwrite(fileID,imag(cr_1),'double');
    fwrite(fileID,real(cr_2),'double'); fwrite(fileID,imag(cr_2),'double');
    fwrite(fileID,real(cr_3),'double'); fwrite(fileID,imag(cr_3),'double');
    fwrite(fileID,real(cc_1),'double'); fwrite(fileID,imag(cc_1),'double');
    fwrite(fileID,real(cc_2),'double'); fwrite(fileID,imag(cc_2),'double');
    fwrite(fileID,real(cc_3),'double'); fwrite(fileID,imag(cc_3),'double');
    fwrite(fileID,real(phi_l_hat),'double');fwrite(fileID,imag(phi_l_hat),'double');
    fwrite(fileID,real(phi_r_hat),'double');fwrite(fileID,imag(phi_r_hat),'double');
    fwrite(fileID,real(phi_c_hat),'double');fwrite(fileID,imag(phi_c_hat),'double');
    fclose(fileID);
  end
end
