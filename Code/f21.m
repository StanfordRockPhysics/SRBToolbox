function [matf] = f21(dx1,dz1,nxmax,nmax,kmax,kmin,dsmax,dsmin,drmax,drmin,dsf)
% function used by bornfilt. Usually need not be called seperately by user.
% see BORNFILT

% finite lines of sources and receivers

% surface seismic

dx = 2*pi*(nmax-1)/(nmax*dx1);
dz = 2*pi*(nxmax-1)/(nxmax*dz1);

% angles
d = sqrt(drmax*drmax+dsf*dsf);
crmax = drmax/d;
srmax = sqrt(1-crmax*crmax);
rmax = acos(drmax/d);
d = sqrt(drmin*drmin+dsf*dsf);
rmin = acos(drmin/d);
crmin = drmin/d;
srmin = sqrt(1-crmin*crmin);
d = sqrt(dsmax*dsmax+dsf*dsf);
csmin = -dsmax/d;
smin = -pi+acos(dsmax/d);
ssmin = -sqrt(1-csmin*csmin);
d = sqrt(dsmin*dsmin+dsf*dsf);
csmax = -dsmin/d;
smax = -pi+acos(dsmin/d);
ssmax = -sqrt(1-csmax*csmax);
xd = crmax-csmin;
yd = srmax-ssmin;
xg = crmin-csmax;
yg = srmin-ssmax;

matf = zeros(nxmax,nmax);
% 4 cases

% first case
if ((rmax <= (pi+smin))&(rmin >= (pi+smax)))
% big loop
for i=1:round(nxmax/2);
   for j=1:nmax;
% matrix (i,j) => spatial coordinates (ii,jj)
ii=(j-(nmax+1)/2)*dx;
jj=(-i+(nxmax+1)/2)*dz;

% first test
if (ii*ii+jj*jj) <= (4*kmax*kmax)
% second and third tests
if (ii*yd-jj*xd) <= 0
if (ii*yg-jj*xg) >= 0
   d1 = (ii-kmin*crmax)*(ii-kmin*crmax) ... 
        +(jj-kmin*srmax)*(jj-kmin*srmax)-kmin*kmin;
   d2 = (ii+kmax*csmin)*(ii+kmax*csmin) ... 
        +(jj+kmax*ssmin)*(jj+kmax*ssmin)-kmax*kmax;
   d5 = (ii+kmin*csmin)*(ii+kmin*csmin) ... 
        +(jj+kmin*ssmin)*(jj+kmin*ssmin)-kmin*kmin;
   d3 = (ii-kmin*crmin)*(ii-kmin*crmin) ... 
        +(jj-kmin*srmin)*(jj-kmin*srmin)-kmin*kmin;
   d4 = (ii+kmax*csmax)*(ii+kmax*csmax) ...
        +(jj+kmax*ssmax)*(jj+kmax*ssmax)-kmax*kmax;
   d6 = (ii+kmin*csmax)*(ii+kmin*csmax) ...
        +(jj+kmin*ssmax)*(jj+kmin*ssmax)-kmin*kmin;
   if (((d6>=0)&(d1>=0))|((d3>=0)&(d5>=0)))&((d2<=0)|(d4<=0))
      matf(i,j) = 1;
   else
   if ((rmin>=(pi+smin))&(rmax<=(smax+pi)))
   theta = acos(ii/sqrt(ii*ii+jj*jj));
   if (d4>=0)&(d2>=0)&(theta>=(pi+smin))&(theta<=(pi+smax))
      matf(i,j) = 1;
   end
   end
   end
% end of the three tests
end
end
end

% end of the big loop
   end
end

% second study
elseif ((rmax >= (pi+smin))&(rmin <= (pi+smax)))
% big loop
for i=1:round(nxmax/2);
   for j=1:nmax;

% matrix (i,j) => spatial coordinates (ii,jj)
ii=(j-(nmax+1)/2)*dx;
jj=(-i+(nxmax+1)/2)*dz;

% first test
if (ii*ii+jj*jj) <= (4*kmax*kmax)
% second and third test
if (ii*yd-jj*xd) <= 0
if (ii*yg-jj*xg) >= 0
   d1 = (ii-kmin*crmax)*(ii-kmin*crmax) ... 
        +(jj-kmin*srmax)*(jj-kmin*srmax)-kmin*kmin;
   d2 = (ii-kmax*crmax)*(ii-kmax*crmax) ... 
        +(jj-kmax*srmax)*(jj-kmax*srmax)-kmax*kmax;
   d5 = (ii+kmin*csmin)*(ii+kmin*csmin) ... 
        +(jj+kmin*ssmin)*(jj+kmin*ssmin)-kmin*kmin;
   d3 = (ii-kmin*crmin)*(ii-kmin*crmin) ... 
        +(jj-kmin*srmin)*(jj-kmin*srmin)-kmin*kmin;
   d4 = (ii-kmax*crmin)*(ii-kmax*crmin) ... 
        +(jj-kmax*srmin)*(jj-kmax*srmin)-kmax*kmax;
   d6 = (ii+kmin*csmax)*(ii+kmin*csmax) ...
        +(jj+kmin*ssmax)*(jj+kmin*ssmax)-kmin*kmin;
   if (((d6>=0)&(d1>=0))|((d3>=0)&(d5>=0)))&((d2<=0)|(d4<=0))
      matf(i,j) = 1;
   else
   if ((rmin<=(smax+pi))&(rmax>=(pi+smin)))
   theta = acos(ii/sqrt(ii*ii+jj*jj));
   if (d4>=0)&(d2>=0)&(theta>=rmax)&(theta<=rmin)
      matf(i,j) = 1;
   end
   end
   end
% end of the three tests
end
end
end

% end of principal loop
   end
end

% third study
elseif ((rmax >= (pi+smin))&(rmin >= (pi+smax)))
% big loop
for i=1:round(nxmax/2);
   for j=1:nmax;

% matrix (i,j) => spatial coordinates (ii,jj)
ii=(j-(nmax+1)/2)*dx;
jj=(-i+(nxmax+1)/2)*dz;

% first test
if (ii*ii+jj*jj) <= (4*kmax*kmax)
% second and third test
if (ii*yd-jj*xd) <= 0
if (ii*yg-jj*xg) >= 0
   d1 = (ii-kmin*crmin)*(ii-kmin*crmin) ... 
        +(jj-kmin*srmin)*(jj-kmin*srmin)-kmin*kmin;
   d2 = (ii-kmax*crmax)*(ii-kmax*crmax) ... 
        +(jj-kmax*srmax)*(jj-kmax*srmax)-kmax*kmax;
   d3 = (ii+kmin*csmin)*(ii+kmin*csmin) ... 
        +(jj+kmin*ssmin)*(jj+kmin*ssmin)-kmin*kmin;
   d4 = (ii+kmax*csmax)*(ii+kmax*csmax) ... 
        +(jj+kmax*ssmax)*(jj+kmax*ssmax)-kmax*kmax;
   if (d1>=0)&(d3>=0)&(d2<=0)&(d4<=0)
      matf(i,j) = 1;
   else
   if (rmax<=(pi+smax))
   if ((d1>=0)&(d3>=0))&((d2<=0)|(d4<=0))
      matf(i,j) = 1;
   end
   theta = acos(ii/sqrt(ii*ii+jj*jj));
   if (d4>=0)&(d2>=0)&(theta>=rmax)&(theta<=(pi+smax))
      matf(i,j) = 1;
   end
   end
   end
% end of the three tests
end
end
end

% end of principal loop
   end
end
% fourth study
elseif ((rmax <= (pi+smin))&(rmin <= (pi+smax)))
% big loop
for i=1:round(nxmax/2);
   for j=1:nmax;

% matrix (i,j) => spatial coordinates (ii,jj)
ii=(j-(nmax+1)/2)*dx;
jj=(-i+(nxmax+1)/2)*dz;

% first test
if (ii*ii+jj*jj) <= (4*kmax*kmax)
% second and third test
if (ii*yd-jj*xd) <= 0
if (ii*yg-jj*xg) >= 0
   d1 = (ii-kmin*crmax)*(ii-kmin*crmax) ... 
        +(jj-kmin*srmax)*(jj-kmin*srmax)-kmin*kmin;
   d2 = (ii-kmax*crmin)*(ii-kmax*crmin) ... 
        +(jj-kmax*srmin)*(jj-kmax*srmin)-kmax*kmax;
   d3 = (ii+kmin*csmax)*(ii+kmin*csmax) ... 
        +(jj+kmin*smax)*(jj+kmin*ssmax)-kmin*kmin;
   d4 = (ii+kmax*csmin)*(ii+kmax*csmin) ... 
        +(jj+kmax*ssmin)*(jj+kmax*ssmin)-kmax*kmax;
   if (d1>=0)&(d3>=0)&(d2<=0)&(d4<=0)
      matf(i,j) = 1;
   else
   if (rmin>=(smin+pi))
   if ((d1>=0)&(d3>=0))&((d2<=0)|(d4<=0))
      matf(i,j) = 1;
   end
   theta = acos(ii/sqrt(ii*ii+jj*jj));
   if (d4>=0)&(d2>=0)&(theta>=(smin+pi))&(theta<=rmin)
      matf(i,j) = 1;
   end
   end
   end
% end of the three tests
end
end
end

% end of principal loop
   end
end
else
  disp(' Erreur ')
  rmax
  rmin
  smax
  smin
end

% symmetry
for i=round(nxmax/2)+1:nxmax;
    for j=1:nmax;
        matf(i,j) = matf(nxmax-i+1,nmax-j+1);
    end
end
%end

% we arbitrarely put 1 at the origin to have a good average value
if rem(nmax,2) == 0
   if rem(nxmax,2) == 0
      matf((nxmax/2),nmax/2) = 1;
      matf((nxmax/2)+1,nmax/2) = 1;
      matf((nxmax/2)+1,nmax/2+1) = 1;
      matf((nxmax/2),nmax/2+1) = 1;
    else
      matf((nxmax+1)/2,nmax/2) = 1;
      matf((nxmax+1)/2,nmax/2+1) = 1;
    end
else
   if rem(nxmax,2) == 0
      matf((nxmax/2),(nmax+1)/2) = 1;
      matf((nxmax/2+1),(nmax+1)/2) = 1;
   else
      matf((nxmax+1)/2,(nmax+1)/2) = 1;
   end
end 
%end
