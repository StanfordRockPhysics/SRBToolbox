function [ss,ttu]=ezseis(vel,dens,varargin);
%[ss,ttu]=ezseis(vel,dens,dx,dz,top,freq,snr);
%
% Quick and easy normal incidence synthetic seismic section.
% Calculates and displays wiggle trace and color plot of seismograms. 
% Inputs (reguired):  VEL, DENS, velocity and density.
% These can be vectors (e.g. well logs) or 2-D matrices of the same size.  
% Optional input parameters: DX, DZ, horizontal and vertical grid spacing,
% TOP, depth to top of velocity image (or log), FREQ, seismic frequency,
% SNR, Noise to Signal energy ratio (usually < 1) [Note: not S/N].
% When not specified on input, the default values are: DX, DZ =1,
% TOP = 0, FREQ = 25 (Hz), SNR = 0 (no random noise added).
% Outputs (optional): SS synthetic seismogram or seismic section, TTU time axis.
% If VEL, DENS are vectors, the single seismogram is repeated 25 times 
% before plotting the section.
%                 
% The algorithm consists of calculation of reflectivity from impedance,
% depth to time conversion, low pass filtering the reflectivity sequence
% based on the value of FREQ, and horizontal averaging over a Fresnel zone.
%
% See also: KENNET, KENFRTT, PSPEC2DSH, BORN2D, BORNFILT, EZSEIS2
% For plotting seismic sections see SEISPLOT, SEISRWB

switch length(varargin)
case 0, dx=1; dz=1; top=0; freq=25; snr=0;
case 1, dx=varargin{1}; dz=1; top=0; freq=25; snr=0;
case 2, dx=varargin{1};dz=varargin{2};top=0; freq=25; snr=0;
case 3, dx=varargin{1};dz=varargin{2};top=varargin{3};freq=25; snr=0;
case 4, dx=varargin{1};dz=varargin{2};top=varargin{3};freq=varargin{4};snr=0;
case 5, dx=varargin{1};dz=varargin{2};top=varargin{3};freq=varargin{4};
        snr=varargin{5};
otherwise, error('Incorrect number of input arguments'); end;


if size(vel,1)==1, vel=vel(:); end;
if size(dens,1)==1, dens=dens(:); end;

imped=vel.*dens; 

rz=0.5*diff(log(imped)); rz=[zeros(1,size(rz,2));rz];
ttop=2*top./(vel(1,:)); tt=cumsum(2*dz./vel)+ttop(ones(size(vel,1),1),:);

ttu=[max(tt(1,:)):0.5*min(min(diff(tt))):min(tt(end,:))].';
rz2tt=zeros(length(ttu),size(rz,2));

for k=1:size(rz2tt,2), rz2tt(:,k)=interp1q(tt(:,k),rz(:,k),ttu); end;

fn=0.5/(ttu(2)-ttu(1)); fc=min(freq/fn,0.99); b=fir1(9,fc);
ss=filtfilt(b,1,rz2tt);

dst=top+dz*0.5*size(vel,1); lmda=mean2(vel)/freq; frnlz=sqrt(dst*lmda);
boxn=max(2,floor(frnlz/dx)); bb=boxcar(boxn); bb=bb./sum(bb);

if size(ss,2)>1, ss=conv2(1,bb,ss,'same'); end;
%ss=filtfilt(bb,1,ss.').'; 
if size(ss,2)==1, ss=ss(:,ones(25,1)); end;

if snr~=0, ss=ss+(sqrt(snr))*std(ss(:))*randn(size(ss)); end;

imagesc(0:dx:dx*size(ss,2)-1,ttu,ss),colormap(rwb),hold on;
seisplot(ss,min(ttu),ttu(2)-ttu(1),0,dx,0);
%set(gca,'xticklabel',num2str(dx*(str2num(get(gca,'xticklabel')))) );



