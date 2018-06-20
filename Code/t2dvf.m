function [sd,d]=t2dvf(st,t0,dt,vf,dx)
% T2DVF -  Convert time section (seismogram) to depth section (seismogram) with
% specified velocity function input. Uses simple linear interpolation of 
% the velocity function.
%
% [SD,D]=T2DVF(ST,T0,DT,VF,DX)
% 
%ST=time section; T0=begining time on time section; 
%DT=time sampling; VF=[2way-time, velocity], 2 column velocity function
%SD=depth section; D=depth axis;
%DX=desired depth sampling (optional); default=min_velocity*dt/2

%Written by T. Mukerji

if ((nargin<5)|isempty(dx)), dx=min(vf(:,2))*dt/2; end;

t=[t0:dt:t0+(size(st,1)-1)*dt]';
if min(t)<vf(1,1), vf=[min(t), vf(1,2);vf]; end;
if max(t)>vf(end,1), vf=[vf;max(t), vf(end,2)]; end;
plot(vf(:,1),vf(:,2)),title('velocity function'),drawnow,pause(2);
vis=interp1(vf(:,1),vf(:,2),t);
%plot(t,vis), title('t-vis'), drawnow, pause;
dnu=[vis(1)*t0/2;  cumsum(vis(2:length(vis)).*(0.5*dt))];
%plot(dnu,st(:,1)),title('dnu-st'), drawnow, pause;
d=[0:dx:max(dnu)]';
sd=interp1(dnu,st,d,'spline');
plot(d,sd(:,1)), drawnow;
