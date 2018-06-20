function dAdx2dz2=fastfurdx2dz2(A,dx,dz)



%----------------------------------------------------------------

% This function calculate the mixed deriviative of a 2-D field U(x,z)

% by the fourier transform {i.e. dU/(dxdz)} 

% The free surface B.C are done by zero pading the displacement field. 

%Written by R. Bachrach
%-----------------------------------------------------------------

ia=0+i;

[Nx,Nz]=size(A);

dKx=2*pi/(Nx*dx);

dKz=2*pi/(Nz*dz);

f2U=fft2(A);


%for j=1:Nz
%if j <=Nz/2+1
%Kz(j)=dKz*(j-1);
%else
%Kz(j)=-dKz*(Nz-j+1);
%end,end
j=1:Nz;
Kz= dKz*(j-1).*(j<=(Nz/2+1)) -dKz*(Nz-j+1).*(j>(Nz/2+1));

%for i=1:Nx
%if i <=Nx/2+1
%Kx(i)=dKx*(i-1);
%else
%Kx(i)=-dKx*(Nx-i+1);
%end, end
i=1:Nx;
Kx= dKx*(i-1).*(i<=(Nx/2+1)) -dKx*(Nx-i+1).*(i>(Nx/2+1));


a1=Kz*0+1;
KKz=Kz'*a1;
KKx=a1'*Kx;

fdUdxdz=-(KKz.^2+KKx.^2).*f2U;




%_____________________________________________________________________

%   inverse fft2 to get the spatial deriviative in the space domain

%---------------------------------------------------------------------
dAdx2dz2=ifft2(fdUdxdz);
