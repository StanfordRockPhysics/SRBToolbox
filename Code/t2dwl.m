function [sd,d]=t2dwl(st,t0,dt,wl,dx)
% T2DWL - Convert time section (seismogram) to depth section (seismogram) with
% specified sonic well log input.
%
% [SD,D]=T2DWL(ST,T0,DT,WL,DX)
%
%ST=time section; T0=begining time on time section; 
%DT=time sampling; WL=[depth, velocity], 2 column sonic log 
%SD=depth section; D=depth axis;
%DX=desired depth sampling (optional)

%Written by T. Mukerji

if (nargin<5), dx=[]; end;

t=cumsum(diff(2*[0;wl(:,1)])./wl(:,2));
vf=[t,wl(:,2)];
%plot(wl(:,2),wl(:,1)),set(gca,'ydir','reverse');
%title('velocity log'), drawnow, pause(2);
%plot(vf(:,1),vf(:,2)),title('2-waytime-vel'), drawnow, pause;
[sd,d]=t2dvf(st,t0,dt,vf,dx);

