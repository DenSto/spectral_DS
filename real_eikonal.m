function real_eikonal(HM_in,TH_in)
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
basename='test/';
if(nargin > 0)
  basename=['TH-HW-',num2str(HM_in)]; %basename_in;
end
if(nargin > 1)
  basename=['d-',num2str(TH_in),'-h-',num2str(HM_in)]; %basename_in;
end
maxNumCompThreads(4)
cm_redblue =flip(cbrewer('div', 'RdBu',129));
cm_inferno=inferno();
c_map=cm_inferno;
c_maprb=cm_redblue;
%if(7==exist(basename,'dir'))
%  return;
%end

xref=0;
LX_a=0.1;

scale=1;                   % quick scale of the linear terms
mu     = scale*[0];         % friction
nu     = scale*[0.01,0];      %viscosity
nz     = scale*[0.0];      %viscosity
%hm     = scale*6.5;            % HM-type wave
delta0 = scale*0.0;        % Terry-Horton i delta - should vary too?

positions = [0,0.25,0.5,0.75,1];

forcing=0.0;              % forcing magnitude
LX   = 2*pi*10;         % X scale
LY   = 2*pi*10;         % Y scale
NX_real=64;        % resolution in x
NY_real=64;        % resolution in y
dt=5e-4;            % time step. Should start small as CFL updated can pick up the pace
pert_size=1e-3;     % size of perturbation
TF=20000.0;           % final time
iF=10000000;  % final iteration, whichever occurs first
iRST=10000;    % write restart dump
i_report=100;
en_print=100;

TSCREEN=1000;    % sreen update interval time (NOTE: plotting is usually slow)
initial_condition='random w';   %'simple vortices' 'vortices' 'random' or 'random w'
AB_order=3; % Adams Bashforth order 1,2,3, or 4 (3 more accurate, 2 possibly more stable) -1 = RK3
linear_term='exact';  % CN, BE, FE, or exact
simulation_type='NL'; % NL, QL, or L (nonlinear, quasilinear and linear respectively)
with_plotting = true; % save plots to file
save_plots = true;   % saveplots to file
zonal=true;
zonal_damping=false;
cfl_cadence=1;
cfl=0.2;
max_dt=5e-2;
safety=0.5;
rng(707296708);
%rng('shuffle');
s=rng;

nbox = length(positions);
rhostar = LX_a/LX;

if(nargin > 0)
  hm = HM_in;
end
if(nargin > 1)
  delta0 = TH_in;
end

mkdir(basename);
mkdir([basename 'plots']); mkdir([basename 'data']);
delete([basename 'log.txt']);

% print log file.
diary([basename 'log.txt'])
diary on;

fig1=0;
%fig2=0;
if(with_plotting)
  fig1=figure(1);
  %  fig2=figure(2);
end

[~,dens,hm,d2ndr2]=get_profiles(positions,xref,0.3,0.2);

% ensure parameters get printed  to log.
fprintf('mu: ');
for i = 1:length(mu), fprintf('  %.05e',mu(i)); end; fprintf('\n');
fprintf('nu: ');
for i = 1:length(nu), fprintf('  %.05e',nu(i)); end; fprintf('\n');
fprintf('HM:');
for i = 1:nbox, fprintf('  %.05e',hm(i)); end; fprintf('\n');
fprintf('TH:%.05e\n',delta0);

fprintf('LX:%.02f LY:%.02f NX:%d NY:%d\n',LX, LY, NX_real, NY_real);
fprintf('scale:%d Tf:%.01f iF:%d\n', scale, TF, iF);
fprintf('Nonlinear:%s\n',simulation_type);
fprintf('random seed:%d AB order:%d CFL step:%d\n',s.Seed, AB_order, cfl_cadence);
fprintf('safety:%f perturbation size:%.05e\n',safety, pert_size);

if(~strcmp(initial_condition,[basename 'restart']))
  for ib=1:nbox
    delete(sprintf('%senergy_%d.dat',basename,ib));
  end
end



NX  = NX_real;
NY  = NY_real;
NX  = 3*NX/2;
NY  = 3*NY/2;

dx=LX/NX;
dy=LY/NY;
dkx=2*pi/LX;
dky=2*pi/LY;
LXnum=0:dx:LX;
LYnum=0:dy:LY;

dX=0;
if(nbox>1)
  dX=(positions(2)-positions(1))/rhostar;
end

minkx= -(NX_real/2 - 1)*dkx;
maxkx=  (NX_real/2 - 1)*dkx;
minky= -(NY_real/2 - 1)*dky;
maxky=  (NY_real/2 - 1)*dky;
kxnum= minkx:dkx:maxkx;
kynum= minky:dky:maxky;


enk=0;
rsk=0;
rek=0;
c_1=zeros(NY,NX,nbox);
c_2=zeros(NY,NX,nbox);
c_3=zeros(NY,NX,nbox);

dt1=0;
dt2=0;
dt3=0;

I=sqrt(-1);

kx  = dkx*ones(1,NY)'*(mod((1:NX)-ceil(NX/2+1),NX)-floor(NX/2)); % matrix of wavenumbers in x direction
ky  = dky*(mod((1:NY)'-ceil(NY/2+1),NY)-floor(NY/2))*ones(1,NX); % matrix of wavenumbers in y direction

% Cutting of frequencies using the 2/3 rule.
% Discard asymmetric N/2 term. Otherwise reality is not enforced.
dealias=abs(kx/dkx) <1/3*NX & abs(ky/dky) <1/3*NY;

zonal_part   = zeros(NY, NX);
fluct_part   =  ones(NY, NX);

if(zonal)
  zonal_part(1,:) = zonal_part(1,:) + ones(1,NX);
  fluct_part = fluct_part - zonal_part;
end

lin_growth = zeros(NY,NX,nbox);

ksquare = kx.^2 + ky.^2;                                % Laplacian in Fourier space
TH = I*delta0*ky;
kmu = build_friction(mu,NX,ksquare);               % Friction
knu = build_viscosity(nu,NX,ksquare);              % Viscosity
ksquare_poisson  =-(ksquare + fluct_part - TH);    % Poisson equation in Fourier space
iksq_poisson = 1.0./ksquare_poisson;
if(ksquare_poisson(1,1) == 0)
  iksq_poisson(1,1) = 0.0;
end

if(~zonal_damping)
  kmu = kmu.*fluct_part;
  knu = knu.*fluct_part + zonal_part.*build_viscosity(nz,NX,ksquare);
end

for ib=1:nbox
  gd = hm(ib)*I*ky.*iksq_poisson/dens(ib);
  lin_growth(:,:,ib)   = kmu   + knu   + gd; % this is the linear growth rate used in computations
end
max(real(lin_growth(:)))

% forcing stuff
kf_min = 16*dkx;
dkf = 0.5*dkx;
forcing_base = abs(ksquare) <= (kf_min+dkf)^2 & abs(ksquare) >= kf_min^2;
nforce = 1.0/sum(forcing_base(:));

% Define initial vorticity distribution
phi_hat = zeros(NY,NX,nbox);
w_hat   = zeros(NY,NX,nbox);
switch lower(initial_condition)
  case {'zero'}
    phi_hat = zeros(NY,NX,nbox);
  case {'vortices'}
    [i,j]=meshgrid(1:NX,(1:NY));
    phi_hat(:,:,1)= exp(-((i*dkx-pi).^2+(j*dky-pi+pi/4).^2)/(0.2))+exp(-((i*dx-pi).^2+(j*dky-pi-pi/4).^2)/(0.2))-0.5*exp(-((i*dkx-pi-pi/4).^2+(j*dky-pi-pi/4).^2)/(0.4));
  case {'simple vortices'}
    phi0=pert_size*sin(dkx.*x).*cos(dky.*y).^2;
    phi_hat=fft2(phi0);
    
    phi=1e-1*phi0;
    phi_hat=fft2(phi);
  case {'random'}
    for ib=1:nbox
      phi = pert_size*(2*rand(NY,NX) - 1);%normally 1e-3
      phi_hat(:,:,ib)=fluct_part.*fft2(phi);
    end
  case {'random w'}
    for ib=1:nbox
      w = pert_size*(2*rand(NY,NX) - 1);%normally 1e-3   
      phi_hat(:,:,ib)=fluct_part.*fft2(w).*iksq_poisson/dens(ib);
    end
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

for ib=1:nbox
  phi_hat(:,:,ib)=dealias.*phi_hat(:,:,ib);
  phi_hat(:,:,ib)=enforceReality(phi_hat(:,:,ib));
  phi_hat(1,1,ib) = 0;
end

mat=0;
if(zonal && nbox > 2)
  mat=zeros(nbox,nbox);
  for ib= 2:(nbox-1)
    mat(ib-1,ib) =  1;
    mat(ib  ,ib) = -2;
    mat(ib+1,ib) =  1;
  end
  mat(1,1) = -2; mat(2,1) = 2; % E(0) = 0 (zero electric field at center)
  mat(nbox-1,nbox) = 1; mat(nbox,nbox) = -2;         % phi(a+ep) = 0 (zero potential at edge)
  mat=mat'/dX^2;
 % inverse=inv(mat);
end  

w_zero = zeros(1,nbox);
force=zeros(NY,NX);
%if(with_plotting && save_plots) %plot growth rate contours
plotgrowth()
%end

%nextScreen = 0;
tic

t=0.;
i=0;
maxV=0;

while t<TF && i<iF

  for ib=1:nbox
    w_hat(:,:,ib) = dens(ib)*ksquare_poisson.*phi_hat(:,:,ib);
  end
  
  [dphidx,d2phidx2]    = get_derivatives(rhostar,positions,phi_hat);
  [dwdx,d2wdx2]        = get_derivatives(rhostar,positions,w_hat);
  
  if(any(isnan(phi_hat(:)))) % Something really bad has happened.
    fprintf('Divergence at iteration: %d\n',i);
    return;
  end
  if(mod(i,i_report)==0)
    fprintf('iteration: %d    dt:%.02e     t:%.03e     step time:%.1d s\n',i,dt,t,toc);
    tic;
    diary off; %flush diary
    diary on;
  end
  if (mod(i,en_print)== 0)
    outputEnergy();
  end
  
  
  if (mod(i,TSCREEN)== 0 && with_plotting)
    % if(t > nextScreen && with_plotting)
    %plotfunc();
    %qqq=3;
    %drawnow
    %nextScreen = t + TSCREEN;
  end
  
  if (mod(i,iRST) == 0)
    %dump(sprintf('restart_%d.bin',rek));
    rek=rek+1;
  end
  
  
  if(forcing ~= 0)
    force=calculate_forcing();
  end
  
  if(mod(i+1,cfl_cadence)==0 && i > 3) % compute new timestep from CFL condition.       
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
  maxV=0;

  for ib=1:nbox
    conv_hat =0;
    % Compute the non-linear term
    if (AB_order < 0 ) % RK
      if(AB_order == -1)
        [phi_hat(:,:,ib)] = RK3(phi_hat(:,:,ib));
      else
        exit();
      end
    else

      w_prev = dens(ib)*(ksquare_poisson.*phi_hat(:,:,ib)  ...
                         + 2*I*kx.*dphidx(:,:,ib) ...
                         + d2phidx2(:,:,ib));

      [conv_hat,u,v,dpdx_l] = calc_Nonlinear(phi_hat(:,:,ib),w_prev, ...
                                       dphidx(:,:,ib),dwdx(:,:,ib));

      abs_u = abs(u);
      abs_v = abs(v);
      abs_up= abs(dpdx_l);
      maxV_s = (max(abs_u(:)) + max(abs_up(:)))/dx + max(abs_v(:))/dy;
      if(maxV_s > maxV)
        maxV = maxV_s;
      end

      conv_hat = dealias.*(conv_hat + nforce*forcing*force/sqrt(dt)); %...
                 %          - 2*nu(1)*I*kx.*fluct_part.*dwdx(:,:,ib) ...
                 %          -   nu(1).*fluct_part.*d2wdx2(:,:,ib));
                         
      conv_hat(1,1) = conv_hat(1,1) + rhostar*nu(1)*fluct_part(1,1)*d2ndr2(ib);                
      

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
          L1=(1 + 0.5*dt*lin_growth(ib))./(1 - 0.5*dt*lin_growth(ib));
          L2=dt./(1 - 0.5*dt*lin_growth(ib));
        case {'BE'} %backward euler
          L1=1.0./(1 - dt*lin_growth(ib));
          L2=dt./(1 - dt*lin_growth(ib));
        case {'EXACT'}%exact
          L1=exp(dt*lin_growth(:,:,ib));
          L2=dt*L1;
          AB2=AB2.*exp((dt1)*lin_growth(ib));
          AB3=AB3.*exp((dt1+dt2)*lin_growth(ib));
          AB4=AB4.*exp((dt1+dt2+dt3)*lin_growth(ib));
        otherwise
          disp('Unknown linear handling type.');
          return;
      end

      % Compute Solution at the next step
      % Implicitly solve the linear term with 2nd order Crank-Nicholson
      w_hat_new = L1.*w_prev;
      w_hat_new = w_hat_new - L2.*(AB1.*conv_hat    - AB2.*c_1(:,:,ib) ... 
                                 + AB3.*c_2(:,:,ib) - AB4.*c_3(:,:,ib));

    end
  
    phi_hat(:,:,ib) = iksq_poisson.*(w_hat_new/dens(ib) ...
                                     -2*I*kx.*dphidx(:,:,ib) ...
                                     - d2phidx2(:,:,ib));
    phi_hat(:,:,ib) = dealias.*enforceReality(phi_hat(:,:,ib));

    w_zero(ib) = w_hat_new(1,1)/dens(ib);
    
    c_3(:,:,ib)=c_2(:,:,ib);
    c_2(:,:,ib)=c_1(:,:,ib);
    c_1(:,:,ib)=conv_hat; 
  end
  
  if(zonal)
    phi_hat(1,1,:) = mat\w_zero';
  end 
  
  t=t+dt;
  i=i+1; 
  
  dt3=dt2;
  dt2=dt1;
  dt1=dt;
  
end

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

  function [y1,u,v,dpdx_l] = calc_Nonlinear(phi_h,w_h, dpdx_h, dwdx_h)
    c_hat = 0; u = 0; v = 0;
    type='NL';
    switch upper(type)
      case {'NL'} %full non-linear
        uhat = -I*ky.*phi_h;
        vhat =  I*kx.*phi_h;
        w_xhat = I*kx.*w_h;
        w_yhat = I*ky.*w_h;
        
        
        % dealiasing here truncates if not padded, other it has no effect
        u  =real(ifft2(dealias.*uhat));      % Compute  y derivative of stream function ==> u
        v  =real(ifft2(dealias.*vhat));      % Compute -x derivative of stream function ==> v
        w_x=real(ifft2(dealias.*w_xhat));      % Compute  x derivative of vorticity
        w_y=real(ifft2(dealias.*w_yhat));
        dpdx_l = real(ifft2(dealias.*dpdx_h));
        dwdx_l = real(ifft2(dealias.*dwdx_h));
        conv     = u.*w_x + v.*w_y + dpdx_l.*w_y + u.*dwdx_l;         % evaluate the convective derivative (u,v).grad(w)
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

  function outputEnergy()
    diary off; %flush diary
    diary on;
    
    for ibb=1:nbox
      energyFile = fopen(sprintf('%senergy_%d.dat',basename,ibb),'a');
      
      if((~strcmp(initial_condition,[basename 'restart'])) && (i == 0))
        fprintf(energyFile,'# [1] t  [2] dt  [3] energy  [4] enstrophy  [5] ZF    [6] DW    [7] flux [8] box\n');
      end

      
      w_curr=dens(ibb)*(ksquare_poisson.*phi_hat(:,:,ibb)  ...
                         + 2*I*kx.*dphidx(:,:,ibb) ...
                         + d2phidx2(:,:,ibb));
      phi_y = I*phi_hat(:,:,ibb).*ky;

      energy = 0.5*real(-conj(phi_hat(:,:,ibb)).*w_curr)/(NX*NY)^2;

      enstrophy = 0.5*w_curr.*conj(w_curr)/(NX*NY)^2;

      energy_tot    = sum(energy(:));
      enstrophy_tot = sum(enstrophy(:));

      flux = 0.5*real(conj(phi_y).*w_curr)/(NX*NY)^2;
      ZF_energy = sum(zonal_part(:).*energy(:));
      DW_energy = sum(fluct_part(:).*energy(:));
      zero_mode = -w_curr(1,1)/(NX*NY);
       
      flux_tot = sum(flux(:));

      fprintf(energyFile,'%e ', t);
      fprintf(energyFile,'%e ', dt);
      fprintf(energyFile,'%e ', energy_tot);
      fprintf(energyFile,'%e ', enstrophy_tot);
      fprintf(energyFile,'%e ', ZF_energy);
      fprintf(energyFile,'%e ', DW_energy);
      fprintf(energyFile,'%e ', flux_tot);     
      fprintf(energyFile,'%e ', zero_mode);
      fprintf(energyFile,'\n');
      
      fclose(energyFile);
      %iso_spectrum(energy,enstrophy);
      %y= DW_energy;
    end
    
  end


  function plotgrowth()
    for ibb=1:nbox
      subplot(1,1,1);
      
      plotg=circshift(dealias.*lin_growth(:,:,ibb),[NY_real/2,NX_real/2]); 
      plotg=plotg(2:NY_real,2:NX_real);
      %force=real(ifft2(f));
      imagesc(kxnum,kynum,real(plotg)); 
      axis equal tight; 
      colorbar;
      %imagesc(kxnum,kynum,f), axis equal tight, colorbar
      set(gca,'Ydir','Normal')
      colormap(c_map);
      title('growth rates');    
      xlabel('kx');
      ylabel('ky');
      if(save_plots)
      %    saveas(gcf,'growth.png');
      end
      drawnow;
    end
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
