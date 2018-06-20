function [mats,matf]=bornfx(mate,fmin,fmax);
% [MATS,MATF]=BORNFX(MATE,FMIN,FMAX);
%
%
% "Interactive" program to compute the filtered seismic image, given the true
% velocity image as the input, and the signal bandwidth. 
% The Born approximation is used to filter the
% input using the characterisitc transfer function appropriate for the
% measurement geometry. Assumes real data, and continuous lines of sources
% and receivers. Number of rows and columns of input matrix should be
% multiple (not necessarily power) of 2.
%
% input: MATE, input matrix of true velocities
%        FMIN, minimum signal frequency (Hz)
%        FMAX, maximim signal frequency (Hz)
% outputs: MATS, filtered seismic image
%          MATF, image response function

% written by Philippe Rio, 1994.



disp(' we assume real data ')
disp(' and continuous lines of receivers and sources')
disp(' and the angles are calculated only at the center of the image')
disp(' We assume size(input)=[2*n,2*m]')
disp(' In case of problem for the conversion  object function -> velocities,')
disp(' the error message "conversion problem" appears')

% to display the input matrix
simate=-1;
while (simate~=1)&(simate~=0)
   simate=input('display the input matrix ? YES = 1, NO = 0 ');
end
% to smooth the filter
simats=-1;
while (simats~=1)&(simats~=0)
   simats=input('smooth filter ? YES = 1, NO = 0 ');
end

%f1max = 100.;f2max = 100.;f3max = 1000.;f4max = 10000.;f5max = 10000.;
%f1min = 05.;f2min =   5.;f3min =  100.;f4min =  10000.;f5min =  1000.;
%f6maxGPR = 2000000.;f6maxSS = 50;
%f6minGPR =  500000.;f6minSS = 05;

f1max = fmax;f2max = fmax;f3max = fmax;f4max = fmax;f5max = fmax;
f1min = fmin;f2min =   fmin;f3min =  fmin;f4min =  fmin;f5min =  fmin;
f6maxGPR = fmax;f6maxSS = fmax;
f6minGPR =  fmin;f6minSS = fmin;


dx = -100;dz = -100;
while (dx<=0)
dx = input('horizontal length of the picture dx = ? ');
end
while (dz<=0)
dz = input('vertical length of the picture dz = ? ');
end
disp(' ')
disp('1 --- surface seismic (every source sends a ray to every receiver) ')
disp('2 --- vertical seismic profile')
disp('3 --- single-well reflection')
disp('4 --- well-to-well transmission')
disp('5 --- log')
disp('6 --- surface seismic (1 ray per couple of source and receiver)')
disp('7 --- GPR (1 ray per couple of source and receiver)')
filtre = -1;
while (filtre~=1)&(filtre~=2)&(filtre~=3)&(filtre~=4)&(filtre~=5)...
      &(filtre~=6)&(filtre~=7)
filtre = input('Type of filter ? ');
end
if filtre == 6
   filtre2 = 2;
end
if filtre == 7
   filtre2 = 1;
end

[nxmax,nzmax] = size(mate);
matf = zeros(nxmax,nzmax);
mats = zeros(nxmax,nzmax);
ave = mean(mean(mate));
% velocity -> object function
mate=1-ave*ave./(mate.*mate);

if filtre == 1
   fmax=f1max;
   m1 = 'OUTPUT DATA for Surface Seismic';
   m2 = 'FILTER for Surface Seismic';
   disp('            ')
   disp('              Y-axis            ')
   disp('                 |              ')
   disp('    dsmin <--  SOURCES  --> dsmax ')
   disp('    drmin <-- RECEIVERS --> drmax ')
   disp('                 |              ')
   disp('--------------------------------- surface')
   disp('  |              |              ')
   disp(' dsf  -----------------------')
   disp('  |   |xxxxxxxxxx|xxxxxxxxxx|   ')
   disp('  ----|----------0----------|--- X-axis')
   disp('      |xxxxxxxxxx|xxxxxxxxxx|')
   disp('      -----------------------   ')
   disp('                 |              ')
   dsf = -1;
   dsmin = 1;dsmax = -1;drmin = 1;drmax = -1;
   while (dsmax<=dsmin)
      dsmin = input('dsmin (<dsmax ) = ? ');
      dsmax = input('dsmax  = ? ');
   end
   while (drmax<=drmin)
      drmin = input('drmin (<drmax ) = ? ');
      drmax = input('drmax  = ? ');
   end
   while (dsf<(dz/2))
      dsf = input('dsf ( > dz/2 , positive) = ? ');
   end
   kmax = 2*pi*f1max/ave;
   kmin = 2*pi*f1min/ave;
   matf = f21(dx,dz,nxmax,nzmax,kmax,kmin,dsmax,dsmin,drmax,drmin,dsf);
elseif filtre == 2
   fmax=f2max;
   m1 = 'OUTPUT DATA for V.S.P.';
   m2 = 'FILTER for V.S.P.';
   disp('         ')
   disp('           Y-axis            ')
   disp('              |              ')
   disp('  dsmin<-- SOURCES --> dsmax    ')
   disp('              |                 ')
   disp('---------------------------------------R-- surface')
   disp('  |           |                        E    ')
   disp(' dsf  -----------------                C')
   disp('  |   |xxxxxxx|xxxxxxx|      |<- drmax E   ')
   disp('  ----|-------0-------|------|---------I-- X-axis')
   disp('      |xxxxxxx|xxxxxxx|      |         V   ')
   disp('      -----------------      |         E   ')
   disp('              |<---- dr ---->|         R   ')
   disp('              |              |<- drmin S     ')
   dr = -1;dsf = -1;
   while (dsf < (dz/2))
      dsf = input('dsf ( > dz/2 , positive) = ? ');
   end
   while (dr < 0)
      dr = input('dr ( > 0 ) = ? ');
   end
   dsmin = 1;dsmax = -1;drmin = 1;drmax = -1;
   while (dsmax<=dsmin)
      dsmin = input('dsmin (<dsmax) = ? ');
      dsmax = input('dsmax  = ? ');
   end
   while (drmax<=drmin)|(drmax>dsf)
      drmin = input('drmin (<drmax) = ? ');
      drmax = input('drmax (<dsf) = ? ');
   end
   kmax = 2*pi*f2max/ave;
   kmin = 2*pi*f2min/ave;
   matf = f22(dx,dz,nxmax,nzmax,kmax,kmin,dsmax,dsmin,drmax,drmin,dsf,dr);
elseif filtre == 3
   fmax=f3max;
   m1 = 'OUTPUT DATA for Single-Well Reflection';
   m2 = 'FILTER for Single-Well Reflection';
   disp('         ')
   disp('                          Y-axis            ')
   disp('                             |              ')
   disp('---------------------------------------- surface')
   disp('                             |              ')
   disp('R  dsmax,drmax ->|<--- dr -->|              ')
   disp('E S              |           |            ')
   disp('C O              |    ---------------   ')
   disp('E U              |    |xxxxxx|xxxxxx|   ')
   disp('I R              | ---|------0------|--- X-axis')
   disp('V C              |    |xxxxxx|xxxxxx|   ')
   disp('E E              |    ---------------   ')
   disp('R S              |           |        ')
   disp('S  dsmin,drmin ->|           |       ')
   dr = -1;
   while (dr < (dx/2))
      dr = input('dr ( > dx/2 , positive) = ? ');
   end
   dsmin = 1;dsmax = -1;drmin = 1;drmax = -1;
   while (dsmax<=dsmin)
      dsmin = input('dsmin (<dsmax) = ? ');
      dsmax = input('dsmax  = ? ');
   end
   while (drmax<=drmin)
      drmin = input('drmin (<drmax) = ? ');
      drmax = input('drmax  = ? ');
   end
   kmax = 2*pi*f3max/ave;
   kmin = 2*pi*f3min/ave;
   matf = f23(dx,dz,nxmax,nzmax,kmax,kmin,dsmax,dsmin,drmax,drmin,dr);
elseif filtre == 4
   fmax=f4max;
   m1 = 'OUTPUT DATA for Well-To-Well Transmission';
   m2 = 'FILTER for Well-To-Well Transmission';
   disp('      ')
   disp('             Y-axis            ')
   disp('                |              ')
   disp('-------------------------------------------------- surface')
   disp('                       |              ')
   disp(' dsmax ->|<--- drs --->|<--- drr --->|<- drmax R ')
   disp('S        |             |             |         E')
   disp('O        |      ---------------      |         C')
   disp('U        |      |xxxxxx|xxxxxx|      |         E')
   disp('R        |   ---|------0------|------|---------I-- X-axis')
   disp('C        |      |xxxxxx|xxxxxx|      |         V ')
   disp('E        |      ---------------      |         E')
   disp('S        |             |             |         R')
   disp(' dsmin ->|             |             |<- drmin S')
   drs = -1;drr = -1;
   while (drs < (dx/2))
      drs = input('drs ( > dx/2 , positive) = ? ');
   end
   while (drr < (dx/2))
      drr = input('drr ( > dx/2 , positive) = ? ');
   end
   dsmin = 1;dsmax = -1;drmin = 1;drmax = -1;
   while (dsmax<=dsmin)
      dsmin = input('dsmin (<dsmax) = ? ');
      dsmax = input('dsmax  = ? ');
   end
   while (drmax<=drmin)
      drmin = input('drmin (<drmax) = ? ');
      drmax = input('drmax  = ? ');
   end
   kmax = 2*pi*f4max/ave;
   kmin = 2*pi*f4min/ave;
   matf = f24(dx,dz,nxmax,nzmax,kmax,kmin,dsmax,dsmin,drmax,drmin,drs,drr);
elseif filtre == 5
   fmax=f5max;
   m1 = 'OUTPUT DATA for Log';
   m2 = 'FILTER for Log';
   disp('         ')
   disp('           Y-axis            ')
   disp('              |              ')
   disp('--------------------------- surface')
   disp('      dsmax ->|<-drmax  R      ')
   disp('   S          |         E    ')
   disp('   O  ----------------- C   ')
   disp('   U  |xxxxxxx|xxxxxxx| E   ')
   disp('---R--|-------0-------|-I--- X-axis')
   disp('   C  |xxxxxxx|xxxxxxx| V   ')
   disp('   E  ----------------- E   ')
   disp('   S          |         R   ')
   disp('      dsmin ->|<-drmin  S     ')
   dsmin = 1;dsmax = -1;drmin = 1;drmax = -1;
   while (dsmax<=dsmin)
      dsmin = input('dsmin (<dsmax) = ? ');
      dsmax = input('dsmax  = ? ');
   end
   while (drmax<=drmin)
      drmin = input('drmin (<drmax) = ? ');
      drmax = input('drmax  = ? ');
   end
   kmax = 2*pi*f5max/ave;
   kmin = 2*pi*f5min/ave;
   matf = f25(dx,dz,nxmax,nzmax,kmax,kmin,dsmax,dsmin,drmax,drmin);
elseif filtre == 6
   m1 = 'OUTPUT DATA for Surface Seismic / GPR';
   m2 = 'FILTER for Surface Seismic / GPR';
   disp('            ')
   disp('              Y-axis            ')
   disp('                 |              ')
   disp('    dsmin <--  SOURCES  --> dsmax ')
   disp('    drmin <-- RECEIVERS --> drmax ')
   disp('                 |              ')
   disp('--------------------------------- surface')
   disp('  |              |              ')
   disp(' dsf  -----------------------')
   disp('  |   |xxxxxxxxxx|xxxxxxxxxx|   ')
   disp('  ----|----------0----------|--- X-axis')
   disp('      |xxxxxxxxxx|xxxxxxxxxx|')
   disp('      -----------------------   ')
   disp('                 |              ')
   dsf = -1;
   dsmin = 1;dsmax = -1;drmin = 1;drmax = -1;
   while (dsmax<=dsmin)
      dsmin = input('dsmin (<dsmax) = ? ');
      dsmax = input('dsmax  = ? ');
   end
   while (drmax<=drmin)
      drmin = input('drmin (<drmax) = ? ');
      drmax = input('drmax  = ? ');
   end
   while (dsf<(dz/2))
      dsf = input('dsf ( > dz/2 , positive) = ? ');
   end
   if filtre2 == 1 
      fmax=f6maxGPR;
      kmax = 2*pi*f6maxGPR/ave;
      kmin = 2*pi*f6minGPR/ave;
   elseif filtre2 == 2
      fmax=f6maxSS;
      kmax = 2*pi*f6maxSS/ave;
      kmin = 2*pi*f6minSS/ave;
   end
   matf = f26(dx,dz,nxmax,nzmax,kmax,kmin,dsmax,dsmin,drmax,drmin,dsf);
end

% smoothing
if simats == 1;
% we assume [2*nxmax',2*nzmax']=size(matf)
disp ('smoothing')
% 125 => 200
rangei=min(nxmax,max(ceil(200*ave/(2*fmax*dx)),10));
rangej=min(nzmax,max(ceil(200*ave/(2*fmax*dz)),10));
% we want to preserve the mean value
mat_mil=zeros(2);
mat_mil=matf(nxmax/2:(nxmax/2+1),nzmax/2:(nzmax/2+1));
matf=smooth(matf,rangei,rangej,filtre);
matf(nxmax/2:(nxmax/2+1),nzmax/2:(nzmax/2+1))=mat_mil;
end

if (rem(nxmax,2)==0)&(rem(nzmax,2)==0)
   mats = real(ifft2(fftshift(fftshift(fft2(mate)).*matf)));
elseif (rem(nxmax,2)~=0)&(rem(nzmax,2)~=0)
   matint = fftshift(fft2(mate)).*matf;
   nx = round(nxmax/2);
   nz = round(nzmax/2);
   for i = 1:nx;
      for j = 1:nz;
      mats(i,j) = matint(i+nx-1,j+nz-1);
      end;
      for j=1:nz-1;
      mats(i,j+nz) = matint(i+nx-1,j);
      end;
   end
   for i = 1:nx-1;
      for j = 1:nz;
      mats(i+nx,j) = matint(i,j+nz-1);
      end;
      for j=1:nz-1;
      mats(i+nx,j+nz) = matint(i,j);
      end;
   end
   matint = mats;
   mats = real(ifft2(matint));
   clear matint;
elseif (rem(nxmax,2)~=0)&(rem(nzmax,2)==0)
   matint = fftshift(fft2(mate)).*matf;
   nx = round(nxmax/2);
   nz = nzmax/2;
   for j = 1:nz;
      for i = 1:nx;
      mats(i,j) = matint(i+nx-1,j+nz);
      mats(i,j+nz) = matint(i+nx-1,j);
      end;
      for i=1:nx-1;
      mats(i+nx-1,j) = matint(i,j+nz);
      mats(i+nx-1,j+nz) = matint(i,j);
      end;
   end
   matint = mats;
   mats = real(ifft2(matint));
   clear matint;
else 
   matint = fftshift(fft2(mate)).*matf;
   nz = round(nzmax/2);
   nx = nxmax/2;
   for i = 1:nx;
      for j = 1:nz;
      mats(i,j) = matint(i+nx,j+nz-1);
      mats(i+nx,j) = matint(i,j+nz-1);
      end;
      for j=1:nz-1;
      mats(i,j+nz-1) = matint(i+nx,j);
      mats(i+nx,j+nz-1) = matint(i,j);
      end;
   end
   matint = mats;
   mats = real(ifft2(matint));
   clear matint;
end


% object function -> velocity
mats=1-mats;
ij=[];
while min(min(mats)) <= 0
   disp('conversion problem')
   [y,j]=min(min(mats));
   [y,i]=min(mats(:,j));
   ij=[ij;i j];
   mats(i,j)=9;
end
mats=ave./sqrt(mats);
mate=ave./sqrt(1-mate);
ij
minmate=0.5*min(min(mate));
maxmate=2*max(max(mate));
[nx1,nz1]=size(mate);
for i=1:nx1
for j=1:nz1
mats(i,j)=max(mats(i,j),minmate);
mats(i,j)=min(mats(i,j),maxmate);
end
end

figure(1);
newplot;
colormap(1-gray);
set(gca,'ydir','reverse')
imagesc(matf);
%contour(matf);
title(m2);
xlabel('horizontal axis (j-coordinate for a matrix (i,j))')
ylabel('vertical axis (i-coordinate for a matrix (i,j))')

figure(2);
newplot;
colormap(gray);
imagesc(mats);
title(m1);
xlabel('horizontal axis (j-coordinate for a matrix (i,j))')
ylabel('vertical axis (i-coordinate for a matrix (i,j))')

if simate==1
   figure(3);
   colormap(gray);
   newplot;
   imagesc(mate);
   title('INPUT');
   xlabel('horizontal axis (j-coordinate for a matrix (i,j))')
   ylabel('vertical axis (i-coordinate for a matrix (i,j))')
end
