function [ssx,ssz,wf,mov]=pspec2dpsv(Lamda,miu,ro,dx,dz,sxzm,gxzm,sw,dt,nt,freq,wfop,imop,movop)
% PSPEC2DPSV - Wave propagation and seismograms in 2-D heterogeneous media
% for elastic P-SV waves. Calculates the elastic wave equation in 2-D
% by using finite difference time integration, fft (pseudo-spectral method)
% for x and z deriviative  
% The formulation is Stress-Velocity formulation. See Kozlov et al, 1990.
%
%[ssx,ssz,wf,mov]=pspec2dpsv(Lamda,miu,ro,dx,dz,sxzm,gxzm,sw,dt,nt,freq,wfop,imop,movop)
%
% LAMDA: Lame' COF. grid of size NX by NZ,  ; 
% MIU:	 Shear modulus grid of size NX by NZ,  
%  RO:  density grid of size NX by NZ, 
%       
% DX:   x (horizontal) grid spacing
% DZ:   z (vertical) grid spacing
% SXZM: mask for source positions; matrix same size as Lamda with 1's at the
%       source positions and 0 everywhere else; plane waves simulated by
%       sources along a line.
% GXZM: receiver positions mask; matrix same size as Lamda, 1's at receiver
%       positions, 0 elsewhere
% SW:   source wavelet; give [] to specify default wavelet;
% DT:   time sampling; giving [] will set dt to default value of 
%       dt = min(dx,dz)/(3*max_velocity)
% FREQ:  Wavelet frequency
% NT:   number of time steps to run
% WFOP: =-1 to not save the wavefields at each step
%       =n (>=1) to save wavefield at every nth timestep
% IMOP: =1 to show wavefield at each time step (increases run time)
%       =-1 to not display wavefield at each time step
% MOVOP:=1 to create movie matrix MOV that can later be played back
%       =n (>1) captures every nth frame, saves memory
% SS:   output seismograms at receiver positions
% WF:   wavefield at every nth time step. The wavefield matrix at each
%       timestep is saved columnwise as a column of the wf matrix
%
%EXAMPLE for using 2D P-SV Matlab Code
% Dry Sand over water saturated sand
%Nx=150;Nz=250;Nl=40;%Nl is the layer where the properties change.
%Lamda=zeros(Nz,Nx)+300^2*1650-2*150^2*1650;
%miu=zeros(Nz,Nx)+150^2*1650;
%ro=zeros(Nz,Nx)+1650;
%Lamda(1:60,1:60)=4000^2*1650-2*150^2*1650;
%Lamda(Nl:Nz,:)=2050*1500^2- 2*150^2*20650;
%miu(Nl:Nz,:)=150^2*1650/2050;
%ro(Nl:Nz,:)= 2050;    
%dx=0.4;dz=0.2;
%sxzm=zeros(Nz,Nx);sxzm(1,75)=1;
%gxzm=sxzm*0;gxzm(3,1:Nx)=1;
%nt=3000;dt=[];sw=[];wfop=0;imop=1;movop=30;
%freq=200;
%[ssx,ssz,wf,mov]=pspec2dpsv(Lamda,miu,ro,dx,dz,sxzm,gxzm,...
%                               sw,dt,nt,freq,wfop,imop,movop);


% Written by Ran Bachrach 2000 after
% pspec2dsh, Modified Tapan Mukerji 1996
% 
%_____________________________________________________________________
% Vx and Vzare the velocities in the horizontal (x) and vertical (z)
% direction in the x-z plane.  
% WW is the absorbing matrix
%______________________________________________________________________

[Nz,Nx]=size(Lamda);
X=(1:Nx)*dx;Z=(1:Nz)*dz;
Vx=zeros(Nz,Nx); Vxold=Vx; Vxnew=Vx; ia=0+1i; tmpp=0;,tmpz=0;
Vz=zeros(Nz,Nx); Vzold=Vz; Vznew=Vz; %
Sxx=Vx;Sxxold=Vx;Sxxnew=Sxx;
Szz=Vx;Szzold=Vx;Szznew=Sxx;
Sxz=Vx;Sxzold=Vx;Sxznew=Sxx;
dVxdz=Vx;dVzdz=Vx;dSxxdz=Vx;dSzzdz=Vx;dSxzdz=Vx;
if isempty(sw), sw=sourcewvlt1; end;
if movop>0, mov=moviein(nt/movop+2); end; 
M=Lamda+2*miu;
alfa=sqrt(M./ro);beta=sqrt(miu./ro);
figure(1)
subplot(2,1,1),imagesc(alfa), colormap(copper), title('P-Velocity'), drawnow; 
subplot(2,1,2),imagesc(beta), colormap(copper), title('S-Velocity'), drawnow; 


if movop>0, mov(:,1)=getframe; end; frame=1; 


if isempty(dt), dt=min(dx,dz)/(3*max(max(alfa))); end;
t0=0.15/1000;  dtf=1/freq/50;sf=round(dtf/dt);

sf
%sf=1;
if sf>1;sw=resample(sw,sf,1);end

%dKx=2*pi/(Nx*dx); dKz=2*pi/(Nz*dz);
figure(2)
dA1=ro*0;dA2=dA1;c1=9/8;c2=-1/24;
for it=1:nt  			% starting to march in time

tt=(it-1)*dt ;
%Vertical stress at the source (dipole) or both vertical and Horizontal (monopole) 
	
	
	Szz=source(sw,tt,dt,Nx,Nz,sxzm)+Szz;% Stess source
	Sxx=source(sw,tt,dt,Nx,Nz,sxzm)+Sxx;


dA1=fastfurdx(Sxx,dx,2);
dA2=fastfurdx(Sxz,dz,1);
Vxnew=Vx+dt*(dA1+dA2)./ro;
dA1=fastfurdx(Sxz,dx,2);
dA2=fastfurdx(Szz,dz,1);

Vznew=Vz+dt*(dA1+dA2)./ro;

dA1=fastfurdx(Vxnew,dx,2);
dA2=fastfurdx(Vznew,dz,1);
Sxxnew=Sxx+dt*(M.*dA1+Lamda.*dA2);
Szznew=Szz+dt*(Lamda.*dA1+M.*dA2);


dA1=fastfurdx(Vznew,dx,2);
dA2=fastfurdx(Vxnew,dz,1);


Sxznew=Sxz+dt*miu.*(dA1+dA2);




if ( (movop>0) & (it==1 | rem(it,movop)==0) )  
figure(2)
subplot(2,1,1)
imagesc(X,Z,Vz),
title(['time step = ',num2str(it)]),drawnow,%pause
subplot(2,1,2)
tx=(1:it)*dt;
imagesc(X,tx,ssz),colormap(bone),drawnow,

;
end;

%if ( (movop>0) & (it==1 | rem(it,movop)==0) ) 
%mov(:,frame+1)=getframe; frame=frame+1;
%end;

%%%%%%%%%%%%%%%%%%%%%%%%5
% BOUNDARY CONDITIONS:  Weightx function

Vznew=weightx(Vznew);
Vxnew=weightx(Vxnew);
Szznew=weightx(Szznew);
Sxznew=weightx(Sxznew);
Sxxnew=weightx(Sxxnew);

%
%Vnew=Weight(Vnew);

if ( (wfop>0) & (it==1 | rem(it,wfop)==0) )
  wf=[wf, Vznew(:)];
end
%ss(it,:)=Vznew(find(gxzm));
tmp=Vznew(find(gxzm));
ssz(it,:)=tmp';
tmp=Vxnew(find(gxzm));
ssx(it,:)=tmp';
Vx=Vxnew; Vz=Vznew;
Sxx=Sxxnew;Szz=Szznew;Sxz=Sxznew;
it
end


function WW=weightx(W)
%function WW=weight(W)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function calculate a matrix of zero-padding and weight
%  For the absorbing boundary and the free surface boundary condition.
%  WW is a matrix of unity with decaying boundaries  
%  at the Z direction after absorbing boundaries.

% Written by Ran Bachrach, Jul. 14, 2000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WW=W;
[Nz,Nx]=size(W);

nw=18;

i=1:nw; 
WW(:,i)=exp(-(0.25*(nw-i(ones(Nz,1),:))).^2).*W(:,i);
i=Nx-nw:Nx;
WW(:,i)=exp(-(0.25*(i(ones(Nz,1),:)-Nx+nw)).^2).*W(:,i);
i=Nz-2*nw:Nz-nw;
WW(i,:)=exp(-(0.25*(i(ones(Nx,1),:)'-Nz+nw+nw)).^2).*W(i,:);
i=Nz-nw:Nz;
WW(i,:)=0*W(i,:);

function dAdx=fastfurdx(A,dx,dir)
%----------------------------------------------------------------
% function dAdx=fastfurdx(A,dx,dir)
% This function calculate the first deriviative of a 2-D field U(x,z)
% by the fourier transform {i.e. dU/(dx)} 
% A-2D matrix size(Nx,Nz)
% dx- spacing 
% dir optional parameter: direction of diff.
% dir= 2 for deriviative in the x direction, 1 for z dir. 

%Written by R. Bachrach
%-----------------------------------------------------------------
ia=0+i;
if nargin<3;dir=1;end
if dir==2;
[Nz,Nx]=size(A);dKx=2*pi/(Nx*dx);
fU=fft(A,[],2);
j=1:Nx;
Kx= dKx*(j-1).*(j<=(Nx/2+1)) -dKx*(Nx-j+1).*(j>(Nx/2+1));
%a1=Kx*0+1;
a1=ones(Nz,1);
KKx=(a1*Kx);
%KKx=a1*Kx;
fdAdx=-(i*KKx).*fU;
dAdx=real(ifft(fdAdx,[],2));
else
[Nz,Nx]=size(A);dKz=2*pi/(Nz*dx);
fU=fft(A,[],1);
j=(1:Nz)';
Kz= dKz*(j-1).*(j<=(Nz/2+1)) -dKz*(Nz-j+1).*(j>(Nz/2+1));
a1=ones(1,Nx);
%a1=Kx*0+1;
%KKx=Kx'*a1;
KKx=(Kz*a1);
fdAdx=-(i*KKx).*fU;
dAdx=real(ifft(fdAdx,[],1));
end
