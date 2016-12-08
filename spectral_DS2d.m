function spectral_DS2d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1-field toy model of the Dimits shift. Similar to the Hasegawa-Mima           %
% equation or the Kuramoto-Sivashinsky equation.                                %
%                                                                               %
% Based off an MIT code originally made by Jean-Christophe Nave.                %
% Modified by Denis St-Onge                                                     %
%                                                                               %
% Laplacian(psi) = w                                                            %
% u = psi_y                                                                     %
% v =-psi_x                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
mkdir('plots'); mkdir('data');
delete('log.txt');

scale=1;         % quick scale of the linear terms
mu   =scale*2;      % friction
nu   =scale*0.01; % viscosity
muZF =scale*0e-4; %scale*2*0; % zonal friction
nuZF =scale*0e-5; %scale*5.0e-4; % zonal viscosity
l    =scale*0;     % Landau-like damping
gamma=scale*2.5;   % linear drive 2.4203
HM=scale*0;		 % HM-type wave
TH=scale*0;      % Terry-Horton i delta
h=1;             % hyperviscosity factor
hZF=2;             % hyperviscosity factor
forcing=0; 		 % forcing magnitude
LX=2*pi*10;      % X scale
LY=2*pi*10;      % Y scale
NX_real=256;     % resolution in x
NY_real=256;     % resolution in y
dt=1e-5;    % time step. Should start small as CFL updated can pick up the pace
pert_size=1e-2; % size of perturbation
TF=1000.0;   % final time
iF=2000000;  % final iteration, whichever occurs first
iRST=10000; % write restart dump
i_report=100;
en_print=100;
TSCREEN=1000; % sreen update interval time (NOTE: plotting is usually slow)
initial_condition='random';   %'simple vortices' 'vortices' 'random' or 'random w' 
AB_order=-1; % Adams Bashforth order 1,2,3, or 4 (3 more accurate, 2 possibly more stable) -1 = RK3
linear_term='exact'; % CN, BE, FE, or exact
simulation_type='NL'; % NL, QL, or L (nonlinear, quasilinear and linear respectively)
padding = true; % 3/2 padding, otherwise 2/3 truncation.
save_plots = true; % save plots to file
system_type='MHM'; % NS, HM, MHM
cfl_cadence=5;
cfl=0.4
max_dt=1e-2;
safety=0.8;
diagnostics=false;

%rng(707296708);
rng('shuffle');
s=rng;


% print log file.
diary('log.txt') 
diary on;

fig1=figure(1);
fig2=figure(2);

% ensure parameters get printed  to log.
fprintf('nu: %.05e  mu:%.05e gamma:%.05e\n',nu,mu,gamma);
fprintf('l: %.05e  HM:%.05e TH:%.05e\n',l,HM,TH);
fprintf('LX:%.02f LY:%.02f NX:%d NY:%d\n',LX, LY, NX_real, NY_real);
fprintf('h:%d scale:%d Tf:%.01f iF:%d\n',h, scale, TF, iF);
fprintf('hZF:%d muZF:%e nuZF:%e\n',hZF,muZF,nuZF);
fprintf('Nonlinear:%s padding:%d System:%s\n',simulation_type, padding,system_type);
fprintf('random seed:%d AB order:%d CFL step:%d\n',s.Seed, AB_order, cfl_cadence);
fprintf('safety:%f\n',safety);

energyFile=0; 
if(strcmp(initial_condition,'restart'))
    energyFile = fopen('energy.dat','a');
else
    energyFile = fopen('energy.dat','w');
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

t=0.;
i=0;
k=0;
c1=0;
c2=0;
c3=0;
dt1=0;
dt2=0;
dt3=0;
w_hat=0; %this is our field.
I=sqrt(-1);

kx=dkx*I*ones(1,NY)'*(mod((1:NX)-ceil(NX/2+1),NX)-floor(NX/2)); % matrix of wavenumbers in x direction 
ky=dky*I*(mod((1:NY)'-ceil(NY/2+1),NY)-floor(NY/2))*ones(1,NX); % matrix of wavenumbers in y direction 

% Cutting of frequencies using the 2/3 rule. 
% Discard asymmetric N/2 term. Otherwise reality is not enforced.
dealias=abs(kx/dkx) <1/3*NX & abs(ky/dky) <1/3*NY;


ksquare=kx.^2+ky.^2;                                % Laplacian in Fourier space
ksquare(1,1)=-1; 
kmu = -(-1)^h*mu * ones(NY,NX)./(ksquare.^(h-1)) ;                            % drag force array
ksquare_viscous=-(-1)^h*nu*(ksquare.^h);            % Viscosity array

ksquare_poisson=ksquare - ones(NY, NX) + TH*ky;		% Poisson equation in Fourier space
if(strcmpi(system_type,'NS'))
    ksquare_poisson=ksquare;		% Poisson equation in Fourier space
end
ksquare_poisson(1,1)=-1;  % fixed Laplacian in Fourier space for Poisson's equation


gd = gamma * ky.^2./ksquare_poisson;

zonal_part = zeros(NY,NX);
zonal_part(1,:) = zonal_part(1,:) + ones(1,NX);
fluct_part = ones(NY,NX) - zonal_part;

lin_growth = -kmu + ksquare_viscous + gd - l*abs(ky); % this is the linear growth rate used in computations


% No damping on zonal modes. Modified Poisson equation (proper adiabatic electron response)
if(strcmpi(system_type,'MHM'))
    kxzf = kx(1,:);
    lin_growth(1,:) = zeros(1,NX) - muZF*ones(1,NX) -(-1)^hZF*nuZF*(kxzf.^(2*hZF)) ;
    ksquare_poisson(1,:) = ksquare_poisson(1,:) + ones(1,NX);
	ksquare_poisson(1,1)=-1;
end
numel(find(lin_growth(:)>0))

lin_trunc = dealias.*lin_growth;
max_growth = max(lin_trunc(:))
max_rate = max(abs(lin_trunc(:)))
0.25/max_rate

lg2 = lin_growth.^2;
lg3 = lin_growth.^3;

% forcing stuff?
kf_min = 16*dkx;
dkf = 0.5*dkx;
forcing_base = abs(ksquare) <= (kf_min+dkf)^2 & abs(ksquare) >= kf_min^2;

% Define initial vorticity distribution
switch lower(initial_condition)
   case {'vortices'}
      [i,j]=meshgrid(1:NX,(1:NY));
      psi=exp(-((i*dkx-pi).^2+(j*dky-pi+pi/4).^2)/(0.2))+exp(-((i*dx-pi).^2+(j*dky-pi-pi/4).^2)/(0.2))-0.5*exp(-((i*dkx-pi-pi/4).^2+(j*dky-pi-pi/4).^2)/(0.4));
    case {'simple vortices'}
        [i,j]=meshgrid((1:NX)*(2*pi/NX),(1:NY)*(2*pi/NY));
        psi=1*sin(i).*cos(j).^2;    
    case {'random psi'}
      psi=pert_size*(2*rand(NX,NY) - 1);%normally 1e-3
      w_hat=ksquare_poisson.*fft2(psi);
    case {'random'}
      w=pert_size*(2*rand(NY,NX)-1);%normally 5e-2
      w_hat=fft2(w);
      w_hat(1,:)=zeros(1,NX);
      %w_hat=zonal_part.*w_hat;
    case {'waves'}
      [i,j]=meshgrid((1:NX)*(2*pi/NX),(1:NY)*(2*pi/NY));
      %w=20*sin(32*j) + 1e-1*cos(13*i);
      w=2e1*sin(32*j) + 5e-2*(2*rand(NY,NX)-1);
      w_hat=fft2(w);
      w_hat(1,:)=fluct_part.*w_hat;
    case {'restart'}
        fileID=fopen('restart.bin','r');
        t=fread(fileID,1,'double');
        dt=fread(fileID,1,'double');
        i=fread(fileID,1,'int');
        k=fread(fileID,1,'int');
        dt1=fread(fileID,1,'double');
        dt2=fread(fileID,1,'double');
        dt3=fread(fileID,1,'double');
        c1=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
        c2=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
        c3=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
        w_hat=fread(fileID,[NY NX],'double')+I*fread(fileID,[NY NX],'double');
        fclose(fileID);        
    otherwise
      disp('Unknown initial conditions !!!');
      return
end

if(padding) % ensure there's no energy in the padded region of k-space.
   w_hat=dealias.*w_hat;
end
w_hat(1,1)=0; % Gauge condition. Should be redundant.
w_hat0=w_hat; % Keep initial conditions for possible future diagnostics.


u=0;
v=0;
U_ZF=0;
force=zeros(NY,NX);
w2=w_hat;

if(save_plots) %plot growth rate contours
	plotgrowth()
end

tic
while t<TF && i<iF
    psi_hat = w_hat./ksquare_poisson;  % Solve Poisson's Equation
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
        outputEnergy();
    end
    
    if (mod(i,TSCREEN)== 0) 
        plotfunc();
    end
    
    if (mod(i,iRST) == 0)
        dump('restart.bin');
    end
    

    if(forcing ~= 0)
       force=calculate_forcing(); 
    end
    
    w_hat_new=0;
    conv_hat =0;
    % Compute the non-linear term
    if (AB_order == -1 ) % RK3
	      w_hat_new = ETDRK3(w_hat);
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
            target_dt=min(cfl*new_dt,max_dt);
            if(target_dt < dt)
                 disp('WARNING: New dt fell below safety.')
            end
            dt=max_dt;        
            while dt > target_dt 
                dt=dt/2.0;
            end
            dt=safety*dt;
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
        target_dt=min(max(0.5*dt,min(CFL*new_dt,1.1*dt)),max_dt);
        if(target_dt < dt)
            disp('WARNING: New dt fell below safety.')
        end
        dt=safety*target_dt;
    end
    
    conv_hat = dealias.*(conv_hat + forcing*kf_min*force/sqrt(dt));
    
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
    
    c3=c2;
    c2=c1;
    c1=conv_hat;
    
    dt3=dt2;
    dt2=dt1;
    dt1=dt;
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

	function f=calculate_forcing()
        f=zeros(NY,NX);
        for i1 = 1:NY
           for j1 = 1:NX
                if(forcing_base(i1,j1) ~=0 )
                    f(i1,j1) = randn(1);
                end    
           end
        end
	end

    function x=RK3(w_h)
        a =[8./15.,5./12.,0.75];
        b =[0,-17.0/60.0,-5.0/12.0];
        Fa =0.0;
        Fb =0.0;
        u_new=w_h;
        for j = 1:3 
            Fb = Fa;
            Fa = calc_Nonlinear(u_new,simulation_type);
            u_new = (1 + 0.5*(a(j)+b(j))*dt*lin_growth)./(1 - 0.5*(a(j)+b(j))*dt*lin_growth).*u_new;
            u_new = u_new + a(j)*Fa +b(j)*Fb;
        end
        x=u_new;
    end   

    function x=ETDRK3(w_h)
        u_new=w_h;
        persistent dto;
        persistent lg;
        if(isempty(dto))
           dto=-1;
           lg=lin_growth(:);
        end
        persistent Q1;
        persistent Q2;
        persistent f1;
        persistent f2;
        persistent f3;
        
        lin_zeroes = (dt*abs(lin_growth)) < 1e-3;
        lin_NZ = ones(NY,NX) - lin_zeroes;
        
        lin_div = lin_growth + lin_zeroes.*ones(NY,NX);
        st=simulation_type;
        eh = exp(0.5*lin_growth*dt);
        e1 = exp(lin_growth*dt);

        M=16; 	
		r=exp(I*pi*((1:M)-0.5)/M);

        if(dto ~= dt)
            LR=dt*lg(:,ones(M,1)) + r(ones(NY*NX,1),:);
            elr=exp(LR);
            Q1 =dt*real(mean(			 (exp(LR/2)-1)./LR							,2));
            Q2 =dt*real(mean(			 (elr-1)./LR                                ,2));
            f1 =dt*real(mean(		 (-4   -LR        +elr.*(4-3*LR+LR.^2))./LR.^3	,2));
            f2 =dt*real(mean(		4*( 2  +LR        +elr.*(-2+LR))./LR.^3			,2));
            f3 =dt*real(mean(		 (-4 -3*LR -LR.^2 +elr.*(4-LR))./LR.^3			,2));
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
        
        lin_zeroes = dt*abs(lin_growth) == 0;
        lin_NZ = ones(NY,NX) - lin_zeroes; 
        
        lin_div = lin_growth + lin_zeroes.*ones(NY,NX);
        st=simulation_type;
        e1 = exp(lin_growth*dt);
        
        A1 = -calc_Nonlinear(u_new,st);
        
        an =     lin_NZ.*(u_new.*e1 + (e1 - ones(NY,NX)).*A1./lin_div);
        an = an + lin_zeroes.*(u_new.*e1 + dt*A1);
        
        A2 = -calc_Nonlinear(an,st);
        
        bn = lin_NZ.*(u_new.*e1 + (e1 - ones(NY,NX)).*(2.0*A2-A1)./lin_div);
        bn = bn + lin_zeroes.*(u_new.*e1 + dt*(2.0*A2-A1));
        
        cnNZ = an + (A2 - A1).*(e1 - ones(NY,NX) - dt*lin_growth)./(dt*lin_div.^2);
        cnZ = an + 0.5*dt*(A2-A1);
        x = dealias.*(lin_zeroes.*cnZ + lin_NZ.*cnNZ);
 end   

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

	function y=calc_Nonlinear(w_h,type)
    conv_hat = 0;
    psi_h = w_h./ksquare_poisson;  % Solve Poisson's Equation
	switch upper(type)
        	case {'NL'} %full non-linear
            	uhat = -ky.*psi_h;
            	vhat =  kx.*psi_h;
            	w_xhat = kx.*w_h;
            	w_yhat = ky.*w_h;
            
            	% dealiasing here truncates if not padded, other it has no effect
            	u  =real(ifft2(dealias.*uhat));      % Compute  y derivative of stream function ==> u
            	v  =real(ifft2(dealias.*vhat));      % Compute -x derivative of stream function ==> v
            	w_x=real(ifft2(dealias.*w_xhat));      % Compute  x derivative of vorticity
            	w_y=real(ifft2(dealias.*w_yhat));
            	conv     = u.*w_x + v.*w_y;         % evaluate the convective derivative (u,v).grad(w)   
            	conv_hat = fft2(conv);              % go back to Fourier space
      
        	case {'QL'} %quasi-linear
            	%for zonal part
            	uphat = -ky.*fluct_part.*psi_h;
            	vphat =  kx.*fluct_part.*psi_h;
        
            	%for fluctuation part
            	U_hat = kx.*zonal_part.*psi_h;
            	U_xxhat = (kx.^2).*U_hat;
            	w_yhat = ky.*fluct_part.*w_h;       
        
            	% dealiasing here truncates if not padded, other it has no effect
            	u    =real(ifft2(dealias.*uphat));      % Compute  y derivative of stream function ==> u
            	v    =real(ifft2(dealias.*vphat));      % Compute -x derivative of stream function ==> v
            	U_ZF =real(ifft2(dealias.*U_hat));      % Compute zonal velocity
            	U_xx =real(ifft2(dealias.*U_xxhat));      % Compute zonal velocity
            	w_y  =real(ifft2(dealias.*w_yhat));
            	conv_fluct = U_ZF.*w_y + u.*U_xx;         % evaluate the convective derivative (u,v).grad(w)   
            	conv_zonal = u.*v;         % evaluate the convective derivative (u,v).grad(w)       
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
        conv_hat = dealias.*conv_hat;
    end
    y=conv_hat;
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
    	force=real(ifft2(f));
    	imagesc(kxnum,kynum,real(plotg)), axis equal tight, colorbar
    	%imagesc(kxnum,kynum,f), axis equal tight, colorbar
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
    function outputEnergy()
        diary off; %flush diary
		diary on;
        
        w_curr=w_hat;
        psi_curr=psi_hat;
     
        enstrophy = 0.5*w_curr.*conj(w_curr);
        energy = 0.5*real(-conj(psi_curr).*w_curr);
     %   flux = -0.5*real(psi_y.*(psi_y.*psi_y + psi_x.*psi_x  + psi.*psi));
             
        enstrophy_tot = sum(enstrophy(:))/(NX*NY)^2;
        energy_tot    = sum(energy(:))/(NX*NY)^2; 
        ZF_energy = sum(sum(zonal_part.*energy))/(NX*NY)^2;
        DW_energy = sum(sum(fluct_part.*energy))/(NX*NY)^2;
       % flux_tot = sum(flux(:))/(NX*NY);
        
        fprintf(energyFile,'%e %e %.15e %e %e %e \n',t,dt,energy_tot,enstrophy_tot,ZF_energy,DW_energy);
        
        
    end
    function plotfunc()
              % Go back in real space omega in real space for plotting
		diary off; %flush diary
		diary on;
        w_curr=w_hat;
        psi_curr=psi_hat;
        psi=real(ifft2(psi_curr));
        w=real(ifft2(w_curr)); 
     
        enstrophy = 0.5*w_curr.*conj(w_curr);
        energy = 0.5*real(-conj(psi_curr).*w_curr);                  
        
		%iso_spectrum(energy,enstrophy);

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
        set(0,'CurrentFigure',fig1);
        subplot(2,2,1)
        imagesc(LXnum,LYnum,psi), axis equal tight, colorbar
        set(gca,'Ydir','Normal')
        colormap(jet)
        title(sprintf('psi t=%.02f',t));
        xlabel('x');
        ylabel('y');
        subplot(2,2,2)
        imagesc(LXnum,LYnum,w), axis equal tight, colorbar
        title(sprintf('vorticity t=%.02f',t));
        set(gca,'Ydir','Normal')
        xlabel('x');
        ylabel('y');
        subplot(2,2,3)
        imagesc(kxnum,kynum, energylog), axis equal tight, colorbar
        title('log10(Energy power spectrum)');    
        set(gca,'Ydir','Normal')
        xlabel('kx');
        ylabel('ky');
        subplot(2,2,4)
        imagesc(kxnum,kynum,wlog), axis equal tight, colorbar
        set(gca,'Ydir','Normal')
        xlabel('kx');
        ylabel('ky');
        title('log10(vorticity/Enstrophy power spectrum)');    
        if(save_plots)
            saveas(gcf,sprintf('plots/fig_%d.png',k));
        end
        drawnow
        set(0,'CurrentFigure',fig2);
        if(1==1)
       %     kx_spec=energylog(NY_real/2,:);
       %     ky_spec=energylog(:,NX_real/2);
       %     kx_spec(NX_real/2)=0.5*(kx_spec(NX_real/2-1)+kx_spec(NX_real/2+1)); 
       %     ky_spec(NY_real/2)=0.5*(ky_spec(NY_real/2-1)+ky_spec(NY_real/2+1)); 
       %     [~,Ix] = max(kx_spec(:));
       %     [~,Iy] = max(ky_spec(:));
       %     if(diagnostics)
       %        disp(sprintf('%.10e %.10e',energy_tot,enstrophy_tot));
       %         disp(sprintf('Maximum DW mode:%.01f    Maximum ZF mode:%.01f',abs(kynum(Iy)), abs(kxnum(Ix))));
       %     end
       %     plot(kxnum,kx_spec,kynum,ky_spec);
       %     legend('Zonal (ky=0)','DW (kx=0)');
       %     title(sprintf('DW/Zonal energy spectrum t=%.02f',t));
       %     xlabel('k');
        else
            kxpos=0:dkx:maxkx;
            loglog(kxpos,radenergy(1:NX_real/2),kxpos,radenstrophy(1:NX_real/2));
            title(sprintf('Radial energy/enstrophy spectrum t=%.02f',t));
            legend('Energy','Enstrophy')
            xlabel('k'); 
        end
        drawnow
        if(save_plots)
            saveas(gcf,sprintf('plots/k_spec_%d.png',k));
        end
        k=k+1;
    end

    function iso_spectrum(energy,enstrophy)
        raden_arr=cat(2,sqrt(abs(ksquare(:))),energy(:),enstrophy(:));
        raden_arr=sortrows(raden_arr);
        radenergy=zeros(max(NX,NY),1);
        radenstrophy=zeros(max(NX,NY),1);
        ik_old=1;
        nk=0;
        counts=zeros(max(NX,NY),1);
        for n = 1:length(ksquare(:))
            ik=round(raden_arr(n,1)/min(dkx,dky)) + 1;
            radenergy(ik) = radenergy(ik) + raden_arr(n,2);
            radenstrophy(ik) = radenstrophy(ik) + raden_arr(n,3);
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
        specFile = fopen(sprintf('data/rad_spec_%d.dat',k),'w');
        radenergy_comp = radenergy(1:l).*comp.';
        radenstrophy_comp = radenstrophy(1:l).*comp.';
        
		ZF_en = zeros(l);
        ZF_ens = zeros(l);
        ZF_en(1:(NX/2)) = energy(1,1:(NX/2));  	 
        ZF_ens(1:(NX/2)) = enstrophy(1,1:(NX/2));

		DW_en = zeros(l);
        DW_ens = zeros(l);
        DW_en(1:(NY/2)) = energy(1:(NY/2),1);  	 
        DW_ens(1:(NY/2)) = enstrophy(1:(NY/2),1);

        for n = 1:l
            fprintf(specFile,'%e %e %e %e %e %e %e\n',range(n),radenergy_comp(n),radenstrophy_comp(n),ZF_en(n),ZF_ens(n), DW_en(n),DW_ens(n));       
        end
        
        fclose(specFile); 
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
        fwrite(fileID,real(c1),'double'); fwrite(fileID,imag(c1),'double');
        fwrite(fileID,real(c2),'double'); fwrite(fileID,imag(c2),'double');
        fwrite(fileID,real(c3),'double'); fwrite(fileID,imag(c3),'double');
        fwrite(fileID,real(w_hat),'double');fwrite(fileID,imag(w_hat),'double');
        fclose(fileID);
    end
end
