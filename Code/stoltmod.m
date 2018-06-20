function data = stoltmod(pad,dz,dx,vel,modl);
%STOLT MODELING
%IMG = STOLTMIG(PAD,DZ,DX,VEL,DATA)
% PAD       : 0=no padding, 1=zero padding on time axis
%		 In order to avoid spatial aliasing, spatial-padding needed.
% DZ,DX     : sampling interval of space
% VEL       : velocity (scalar). Can handle only constant velocity.
% MODL/DATA : input/output. MODL is the reflectivity section, DATA is the
%                           stack section.
%
% See also:   STOLTMIG, GAZADJ, KIRCHMIG, KIRCHADJ


% Due to the sign convention of Claerbout's book (Basic Earth Imaging),
% Fourier transform of time axis == Inverse transform in MATLAB.
% Need ft1axis.m & ft2axis.m
%
% 1999. 12.  Youngseuk Keehm
%

% declare tmp. variables and parameters
sig1=-1; sig2=-1;
[nz,nx]=size(modl);
if(pad==1) modl=[modl; zeros(nz,nx)]; nz=2*nz; end

nt=nz; nw=nt; nkx=nx; nkz=nz; dt=dz*2/vel; vhalf2=vel*vel/4;
w0=-pi/dt; kz0=-pi/dz; kx0=-pi/dx;
dw=2*pi/(nt*dt); dkz=2*pi/(nz*dz); dkx=2*pi/(nx*dx); 
data = complex(zeros(nt,nx),0);

modl=ft2axis(0,sig2,modl);
modl=ft1axis(0,sig1,modl);
for ikx=2:nkx
   kx=kx0+(ikx-1)*dkx;
   for iw=2:nw
      w=w0+(iw-1)*dw;
      kz=w*w/vhalf2-kx*kx;
      if(kz>=0)	kz=-sign(w)*sqrt(kz);
      else			kz=0;
      end
      cwkx=cinterp1(kz,nkz,kz0,dkz,modl(:,ikx));
      if(kz~=0) 	data(iw,ikx)=cwkx;
      else 			data(iw,ikx)=0;
      end
   end
end
data=ft1axis(1,-sig1,data);
data=ft2axis(1,sig2,data);
data=real(data);
if(pad==1) data=data(1:nz/2,:); end



function cbx=cinterp1(x,nx,x0,dx,cb);
%LINEAR INTERPOLATION

xc =(x-x0)/dx;
ixc=fix(xc);
fraction = xc - ixc;
ix = 1+ixc;
if(ix<1)			cbx=cb(1);
elseif(ix+1>nx) 	cbx=cb(nx);
else 				cbx=(1-fraction)*cb(ix) + fraction*cb(ix+1);
end

  
