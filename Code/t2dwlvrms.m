function [sd,d,vf]=t2dwlvrms(st,t0,dt,wl,dx)
% T2DWLVRMS - Convert time section (seismogram) to depth section (seismogram) 
%with specified sonic well log input.
% The well log velocity (interval) is converted to RMS velocity
%
% [SD,D,VF]=T2DWLVRMS(ST,T0,DT,WL,DX)
%
%ST=time section; T0=begining time on time section; 
%DT=time sampling; WL=[depth, velocity], 2 column sonic log 
%SD=depth section; D=depth axis;
%VF=[2-way_traveltime, Vrms], two column matrix with rms velocity
%DX=desired depth sampling (optional)

%Written by T. Mukerji
% Modified by R. Bachrach 

if (nargin<5), dx=[]; end;

t=cumsum(diff([0;wl(:,1)])./wl(:,2));	%1WTT
Vrms=sqrt(cumsum(wl(:,2).*(diff([0;wl(:,1)]))./t));

vf=[2*t,Vrms];				%2WTT

plot(wl(:,2),wl(:,1)),set(gca,'ydir','reverse');
title('velocity log'), drawnow, pause(2);
plot(vf(:,1),vf(:,2)),title('2-waytime-RMS vel'), drawnow, pause(2);
[sd,d]=t2dvf(st,t0,dt,vf,dx);
