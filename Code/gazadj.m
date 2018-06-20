function out=gazadj(adj, dt, dx, v, in)
%PHASE_SHIFT GAZDAG MODELING/MIGRATION
%OUT = GAZADJ(ADJ,DT,DX,V,IN)
% ADJ    : 0=forward modeling, 1=migration
% DT,DX  : time sampling, horizontal spacing
% V      : velocity (scalar is OK, if const. velocity)
%          velocity can vary with depth/time, but is constant laterally.
%          V(z) should be a vector with same length as rows of IN.
% IN,OUT : input/output data. In forward modeling, IN is the reflectivity
%	   image, OUT is the stack section. In the migration mode, IN is
%          the stack section, and out is the migrated image.
%
%	   e.g.
%		>> data=zeros(50,50);
%		>> data(24:26,24:26)=1;
%		>> model=gazadj(0,5e-3,20,2300,data);
%		>> imagesc(model)
%		>> migr=gazadj(1,5e-3,20,2300,model); 
%		>> imagesc(migr) 
%
% See also:  STOLTMIG, STOLTMOD, KIRCHMIG, KIRCHADJ


% Due to the sign convention of Claerbout's book (Basic Earth Imaging),
% Fourier transform of time axis == Inverse transform in MATLAB 
% Need ft1axis.m & ft2axis.m
%
% 1999. 12.  Youngseuk Keehm
%

[nt, nx]=size(in);
if(length(v)==1) 			v=v*ones(nt,1); 		end
if(adj==1) 	data=in; 	modl=zeros(nt,nx); 
else			modl=in;	data=zeros(nt,nx);
end

nz=nt;			qi=.5/(nt*dt);
w0 =-pi/dt; 	dw =2*pi/(nt*dt);
kx0=-pi/dx; 	dkx=2*pi/(nx*dx); 

if(adj==0)	modl=ft2axis(0,-1,modl);
else		data=ft2axis(0,-1,data);	data=ft1axis(0,1,data);
end

for ikx=2:nx
   kx=kx0+(ikx-1)*dkx;
   for iw=2:1+fix(nt/2)
      w=w0+(iw-1)*dw;
      if(adj == 0)
         data(iw,ikx)=modl(nz,ikx);
         for iz=nz-1:-1:1
            data(iw,ikx)=data(iw,ikx)*eiktau(dt,w,v(iz)*kx,qi)+modl(iz,ikx);
         end
      else
         cup=data(iw,ikx);
         for iz=1:nz
            modl(iz,ikx)=modl(iz,ikx)+cup;
            cup=cup*conj(eiktau(dt,w,v(iz)*kx,qi));
         end
      end
   end
end
if(adj==0)	data=ft1axis(1, 1,data); data=ft2axis(1,-1,data); out=real(data);
else		modl=ft2axis(1,-1,modl); out=real(modl);
end


function eik = eiktau(dt, w, vkx, qi)

eik=exp( -dt*sqrt(complex(qi, -w)^2 + vkx*vkx/4) );
