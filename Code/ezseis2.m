function [ss1,ss2,ttu]=ezseis2(vel1,dens1,vel2,dens2,dx,dz,top,freq,snr);
%[ss1,ss2,ttu]=ezseis2(vel1,dens1,vel2,dens2,dx,dz,top,freq,snr);
%
% Quick and easy normal incidence pair of synthetic seismic section. 
% Calculates and displays wiggle trace and color plot of the two sections.
% Inputs (reguired):  VEL1, DENS1, VEL2, DENS2, velocities and densities.
% These can be vectors (e.g. well logs) or 2-D matrices.
% The 1 and 2 could correspond to two different situations e.g. before
% and after fluid substitution.
% Optional input parameters: DX, DZ, horizontal and vertical grid spacing,
% TOP, depth to top of velocity image (or log), FREQ, seismic frequency,
% SNR, Noise to Signal energy ratio (usually < 1) [Note: not S/N].
% When not specified on input, the default values are: DX, DZ =1, 
% TOP = 0, FREQ = 25 (Hz), SNR = 0 (no random noise added).
% Outputs (optional): SS1, SS2, synthetic seismograms or seismic sections,
% with common time axis, TTU. The amplitudes are normalized by the same 
% value for both sections, and a single colormap is used,
% so they can be directly compared.
% If VEL, DENS are vectors, the single seismogram is repeated 25 times 
% before plotting the section.
%                 
% The algorithm consists of calculation of reflectivity from impedance,
% depth to time conversion, low pass filtering the reflectivity sequence
% based on the value of FREQ, and horizontal averaging over a Fresnel zone.
%
% See also: EZSEIS, KENNET, KENFRTT, PSPEC2DSH, BORN2D, BORNFILT
% For plotting seismic sections see SEISPLOT, SEISRWB

%Written by T. Mukerji, 1998

if nargin~=9
switch nargin
case 4, dx=1; dz=1; top=0; freq=25; snr=0;
case 5, dz=1; top=0; freq=25; snr=0;
case 6, top=0; freq=25; snr=0;
case 7, freq=25; snr=0;
case 8, snr=0;
otherwise, error('Incorrect number of input arguments'); end;
end;

if size(vel1,1)==1, vel1=vel1(:); vel2=vel2(:); end;
if size(dens1,1)==1, dens1=dens1(:); dens2=dens2(:); end;

imped1=vel1.*dens1; rz1=0.5*diff(log(imped1)); rz1=[zeros(1,size(rz1,2));rz1];
imped2=vel2.*dens2; rz2=0.5*diff(log(imped2)); rz2=[zeros(1,size(rz2,2));rz2];
ttop1=2*top./(mean(vel1(1,:))); tt1=cumsum(2*dz./vel1)+ttop1*ones(size(vel1));
ttop2=2*top./(mean(vel2(1,:))); tt2=cumsum(2*dz./vel2)+ttop2*ones(size(vel2));

ttu0=min([tt1(1,:) tt2(1,:)]); ttuend=max([tt1(end,:) tt2(end,:)]);
dttu=0.5*min( min(min(diff(tt1))), min(min(diff(tt2))) );
ttu=[ttu0:dttu:ttuend].';

nc1=size(tt1,2); nc2=size(tt2,2);
tt1=[(min(ttu)-dttu)*ones(1,nc1);tt1;(max(ttu)+dttu)*ones(1,nc1)];
tt2=[(min(ttu)-dttu)*ones(1,nc2);tt2;(max(ttu)+dttu)*ones(1,nc2)];
rz1=[zeros(1,size(rz1,2));rz1;zeros(1,size(rz1,2))];
rz2=[zeros(1,size(rz2,2));rz2;zeros(1,size(rz2,2))];

rz2tt1=zeros(length(ttu),size(rz1,2)); rz2tt2=zeros(length(ttu),size(rz2,2));

for k=1:size(rz2tt1,2), rz2tt1(:,k)=interp1q(tt1(:,k),rz1(:,k),ttu); end;
for k=1:size(rz2tt2,2), rz2tt2(:,k)=interp1q(tt2(:,k),rz2(:,k),ttu); end;

fn=0.5/(ttu(2)-ttu(1)); fc=min(freq/fn,0.99); b=fir1(9,fc);
ss1=filtfilt(b,1,rz2tt1); ss2=filtfilt(b,1,rz2tt2);

dst=top+dz*0.5*size(vel1,1); lmda=mean2(vel1)/freq; frnlz=sqrt(dst*lmda);
boxn=max(2,floor(frnlz/dx)); bb=boxcar(boxn); bb=bb./sum(bb);
 
if size(ss1,2)>1, 
ss1=conv2(1,bb,ss1,'same');
end;
if size(ss2,2)>1, 
ss2=conv2(1,bb,ss2,'same');
end;
if size(ss1,2)==1, ss1=ss1(:,ones(25,1)); end;
if size(ss2,2)==1, ss2=ss2(:,ones(25,1)); end;

%if size(ss1,2)>1, ss1=filtfilt(b,1,ss1.').'; end;
%if size(ss2,2)>1, ss2=filtfilt(b,1,ss2.').'; end;

if snr~=0, 
ss1=ss1+(sqrt(snr))*std(ss1(:))*randn(size(ss1));
ss2=ss2+(sqrt(snr))*std(ss2(:))*randn(size(ss2));
end;

caxs=[ min([ss1(:);ss2(:)]),max([ss1(:);ss2(:)]) ];
if nargout==0
subplot(211)
%imagesc(0:dx:dx*size(ss1,2)-1,ttu,ss1),colormap(rwb),caxis(caxs),hold on;
seisrwb(ss1,min(ttu),ttu(2)-ttu(1),0,dx,0,ss1);
%set(gca,'xticklabel',num2str(dx*(str2num(get(gca,'xticklabel')))) );
title('Seismograms before'); xlabel('distance'); ylabel('time (s)');
subplot(212)
%imagesc(0:dx:dx*size(ss2,2)-1,ttu,ss2),colormap(rwb),caxis(caxs),hold on;
seisrwb(ss2,min(ttu),ttu(2)-ttu(1),0,dx,0,ss1);
%set(gca,'xticklabel',num2str(dx*(str2num(get(gca,'xticklabel')))) );
title('Seismograms after'); xlabel('distance'); ylabel('time (s)');
end;

