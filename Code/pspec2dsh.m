function [ss,wf,mov]=pspec2dsh(alfa,dx,dz,sxzm,gxzm,sw,dt,nt,wfop,imop,movop)
% PSPEC2DSH - Wave propagation and seismograms in 2-D heterogeneous media
% for SH/acoustic waves. Calculates the elastic (SH) wave equation in 2-D
% by using finite difference time integration and spectral
% method for the spatial deriviatives
%
%[SS,WF,MOV]=PSPEC2DSH(ALFA,DX,DZ,SXZM,GXZM,SW,DT,NT,WFOP,IMOP,MOVOP)
%
%
% ALFA: velocity grid of size NX by NZ, NX and NZ should be even; 
%       and NX=NZ
% DX:   x (horizontal) grid spacing
% DZ:   z (vertical) grid spacing
% SXZM: mask for source positions; matrix same size as alfa with 1's at the
%       source positions and 0 everywhere else; plane waves simulated by
%       sources along a line.
% GXZM: receiver positions mask; matrix same size as alfa, 1's at receiver
%       positions, 0 elsewhere
% SW:   source wavelet; give [] to specify default wavelet;
% DT:   time sampling; giving [] will set dt to default value of 
%       dt = min(dx,dz)/(3*max_velocity)
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

% Initial version written by Ran Bachrach 1996
% Modified Tapan Mukerji 1996
%_____________________________________________________________________
% W and U are the displacement in the horizontal (U) and vertical (W)
% direction in the x-z plane.  V(x,z,t) is the displacement in the y direction.
% WW is the absorbing matrix
%______________________________________________________________________

[Nx,Nz]=size(alfa);
V=zeros(Nx,Nz); Vold=V; Vnew=V; ia=0+1i; tmpp=0;,tmpz=0;
if isempty(sw), sw=sourcewvlt; end;

if movop>0, mov=moviein(nt/movop+2); end; 
imagesc(alfa), colormap(copper), title('Velocity'), drawnow; 
if movop>0, mov(:,1)=getframe; end; frame=1; 

alfa2=alfa.^2;
beta2=alfa2;

if isempty(dt), dt=min(dx,dz)/(3*max(max(alfa))); end;
t0=0.15/1000;  

dKx=2*pi/(Nx*dx); dKz=2*pi/(Nz*dz);

for it=1:nt  			% starting to march in time
tt=(it-1)*dt ;
	V=source(sw,tt,dt,Nx,Nz,sxzm)+V;

dVdx2dz2=real(fastfurdx2dz2(V,dx,dz)); % second derivatives in fourier space

tmpp=(beta2).*dVdx2dz2;

Vnew=tmpp*dt*dt+2*V-Vold;
if imop>0, imagesc(Vnew), colormap(jet);
title(['time step = ',num2str(it)]), drawnow; end;
if ( (movop>0) & (it==1 | rem(it,movop)==0) ) 
mov(:,frame+1)=getframe; frame=frame+1;
end;

Vnew=Weight(Vnew);

if ( (wfop>0) & (it==1 | rem(it,wfop)==0) )
  wf=[wf, Vnew(:)];
end
%ss(it,:)=Vnew(gxzm)';
tmp=Vnew(find(gxzm));
ss(it,:)=tmp';
Vold=V; V=Vnew;
end
