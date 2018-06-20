function U=source(sw,t,dt,Nx,Nz,sxzm)
%function U=source(sw,t,dt,Nx,Nz,sxzm)
%function U=source(sw,t,dt,Nx,Nz,j0,i0)
% U=source(t,Nx,Nz,i0,j0)
%t=time, Nx,Nz size of velocity matrix, i0,j0=source point.
%----------------------------------------------------
%
% This file is the source function for the pspec2dsh program

% Written by R. Bachrach
% Modified T. Mukerji
%
%-------------------------------------------------------

%U=zeros(N);

%dx=1;dz=1;

%j0=N/2;,i0=N/2; t0=0.15/1000;
t0=0.15/1e3; tau=dt/100;
%sw=sourcewvlt;

%for i=1:N
%for j=1:N
%Gxz=exp(-0.5*((i-i0)^2+(j-j0)^2));
%ft=exp(-1000*(t-t0)^2)*(t-t0);
%
%U(i,j)=Gxz*ft;
%
%end,end

[i,j]=meshgrid(1:Nx,1:Nz);
%U=exp(-0.5*((i-i0).^2+(j-j0).^2))*exp(-1000*(t-t0)^2)*(t-t0);
%U=exp(-0.5*((i-i0).^2+(j-j0).^2))*exp(-tau*(t-t0)^2)*(t-t0);
%U=exp(-0.5*((i-i0).^2+(j-j0).^2))*sw(t/dt+1);

%Gxz=zeros(Nx,Nz); Gxz(j0,i0)=1;
Gxz=sxzm;
U=Gxz*sw(round(t/dt)+1);
