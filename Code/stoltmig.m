function img = stoltmig(pad,dt,dx,vel,data);
%STOLT MIGRATION
%IMG = STOLTMIG(PAD,DT,DX,VEL,DATA)
% PAD      : 0=no padding, 1=zero padding on time axis
%		In order to avoid spatial aliasing, spatial-padding needed.   
% DT,DX    : sampling interval of time & space
% VEL      : velocity (scalar). Stolt migration can handle only constant VEL.
% DATA/IMG : input/output. DATA is the stack section, IMG is the migrated image.
%
% See also:  STOLTMOD, GAZADJ, KIRCHMIG, KIRCHADJ 

% Due to the sign convention of Claerbout's book (Basic Earth Imaging),
% Fourier transform of time axis == Inverse transform in MATLAB.
% Need ft1axis.m & ft2axis.m
%
% 1999. 12.  Youngseuk Keehm
%

% declare tmp. variables and parameters
sig1=1; sig2=-1;
[nt,nx]=size(data);
if(pad==1) data=[data; zeros(nt,nx)]; nt=2*nt; end

nz=nt; nw=nt; nkx=nx; nkz=nz; dz=dt*vel/2; vhalf=vel/2;
w0=-pi/dt; kz0=-pi/dz; kx0=-pi/dx;
dw=2*pi/(nt*dt); dkz=2*pi/(nz*dz); dkx=2*pi/(nx*dx); 
img = complex(zeros(nt,nx),0);

data=ft2axis(0,sig2,data);
data=ft1axis(0,sig1,data);
for ikx=2:nkx
   	kx=kx0+(ikx-1)*dkx;
   for ikz=2:nkz
      kz=kz0+(ikz-1)*dkz;
      w=-sign(kz)*vhalf*sqrt(kx*kx+kz*kz);
      ckzkx=cinterp1(w,nw,w0,dw,data(:,ikx));
      if(w~=0) 	img(ikz,ikx)=ckzkx*abs(kz/w);
      else 		img(ikz,ikx)=0;
      end
   end
end
img=ft1axis(1,-sig1,img);
img=ft2axis(1,sig2,img);
img=real(img);
if(pad==1) img=img(1:nt/2,:); end


function cbx=cinterp1(x,nx,x0,dx,cb);
%LINEAR INTERPOLATION

xc =(x-x0)/dx;
ixc=fix(xc);
fraction = xc - ixc;
ix = 1+ixc;
if(ix<1)			cbx=cb(1);
elseif(ix+1>nx) 	cbx=cb(nx);
else				cbx=(1-fraction)*cb(ix) + fraction*cb(ix+1);
end

