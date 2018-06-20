function [wf,wi,ws] = born2d(wvlt,v,c0,gx,dx,dz,dt)
%BORN APPROXIMATION
%function [wf,wi,ws] = born2d(wvlt,v,c0,gx,dx,dz,dt)
%BORN2D calculates acoustic 2-d plane-wave propagation 
%in heterogeneous media using Born Approximation with 
%the 2-d Green's function (Hankel function).
%The plane wave source is located at the left (column 1)
%of the heterogeneous medium.
%WVLT is the source wavelet
%V is the input velocity field (2-d matrix)
%C0 is the constant background velocity
%GX is the x co-ordinate (in column nos.) of the receiver line.
%The forward scattered field can be obtained by specifying
%GX to be the total no. of columns in V.
%The back scattered field can be obtained by specifying
%GX = 1. (receivers located at the source position).
%DX, DZ are the grid spacings along x and z axes.
%DT is the time sampling of the input wavelet.
%WF is the total wavefield at the receivers.
%WI is the incident field at the receivers.
%WS is the scattered field at the receivers.
%
% e.g.  w=sourcewvlt; w=w(1:1024);
%       vel=2000*ones(50,50);
%       vel(10:15,10:15)=1900; vel(40:45,10:15)=1900;
%       vel(21:30,21:30)=2200;
%       [wf,wi,ws]=born2d(w,vel,2000,50,100,50,1e-3); 
%       subplot(1,2,1), imagesc(vel');
%       subplot(1,2,2), imagesc(ws');

%Written by M. Sengupta, Tapan Mukerji 1996

[nz,nx] = size(v);
nt = length(wvlt);
fw = fft(wvlt);
freq = [1:(nt/2-1)]./(nt*dt);
fw = fw(2:nt/16);            %reduce 16 to cover more of the higher freqs.
fw = fw(:).'; 
f = c0^2./v.^2 -1;
mask=f~=0;
fmask=f(mask);
k0 = 2*pi*freq./c0;
u = zeros(nz,length(k0));
ui = zeros(nz,length(k0));
us = zeros(nz,length(k0));
k0=k0(1:length(fw));
[x,z] = meshgrid(0:dx:(nx-1)*dx,0:dz:(nz-1)*dz);
rs = x;
rs=rs(mask);
k0m=k0(ones(length(fmask),1),:);
fwm = fw(ones(length(fmask),1),:);
uinc = fw.*exp(-i*k0*(gx-1)*dx);

hwaitbar=waitbar(0,'Computing Born 2d...');
for gz = 1:nz
rr = sqrt(((gx-1)*dx-x).^2+((gz-1)*dz-z).^2);
rr=rr(mask);
rrm = rr; rrm = rrm(:,ones(1,length(k0)));
%gfunr=(-i/4)*(besselj(0,k0m.*rrm)-i*bessely(0,k0m.*(rrm+eps*(rrm==0)))).*fwm;
gfunr=(-i/4)*(besselh(0,2,k0m.*(rrm+eps*(rrm==0)))).*fwm;
              %changed 1999 to use besselh = besselj-i*bessely
              %single besselh faster than besselj+bessely
rsm = rs; rsm = rsm(:,ones(1,length(k0)));
gfuns = exp(-i*k0m.*rsm);
ff=fmask;
fm = ff(:,ones(1,length(k0)));
usc = (-k0m.^2.*fm.*gfunr.*gfuns);
nuone = isnan(usc);
index = find(nuone==1);
usc(index)=zeros(size(index));
usct=sum(usc);
u(gz,1:length(k0)) = uinc+usct;
ui(gz,1:length(k0)) = uinc;
us(gz,1:length(k0)) = usct;
waitbar(gz/nz,hwaitbar);
end

u =[zeros(size(u,1),1),u,zeros(size(u,1),1),conj(fliplr(u))];
ui =[zeros(size(ui,1),1),ui,zeros(size(ui,1),1),conj(fliplr(ui))];
us =[zeros(size(us,1),1),us,zeros(size(us,1),1),conj(fliplr(us))];
wf = ifft(u.').';
wf = real(wf);
wi = ifft(ui.').';
wi = real(wi);
ws = ifft(us.').';
ws = real(ws);

delete(hwaitbar);
