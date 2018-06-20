function[matf]=f24(dx1,dz1,nxmax,nmax,kmax,kmin,dsmax,dsmin,drmax,drmin,drs,drr)
% function used by bornfilt. Usually need not be called seperately by user.
% see BORNFILT


% finite lines of sources and receivers

% WTW transmission

im = sqrt(-1);
dx = 2*pi*(nmax-1)/(nmax*dx1);
dz = 2*pi*(nxmax-1)/(nxmax*dz1);

% angles
d = sqrt(drmax*drmax+drr*drr);
srmax = drmax/d;
crmax = sqrt(1-srmax*srmax);
d = sqrt(drmin*drmin+drr*drr);
srmin = drmin/d;
crmin = sqrt(1-srmin*srmin);
d = sqrt(dsmin*dsmin+drs*drs);
ssmax = -dsmin/d;
csmax = sqrt(1-ssmax*ssmax);
d = sqrt(dsmax*dsmax+drs*drs);
ssmin = -dsmax/d;
csmin = sqrt(1-ssmin*ssmin);
xa = crmax - csmax;
ya = srmax - ssmax;
xb = crmin - csmin;
yb = srmin - ssmin;

matf = zeros(nxmax,nmax);

% big loop
for i=1:round(nxmax/2);
    for j=1:nmax;
% matrix (i,j) => spatial coordinates (ii,jj)
ii=(j-(nmax+1)/2)*dx;
jj=(-i+(nxmax+1)/2)*dz;

% first test
if (ii*ii+jj*jj) <= (4*kmax*kmax)
% second and third tests
if (ii*ya-jj*xa)*ya >= 0
if (ii*yb-jj*xb)*yb <= 0
   d1 = (ii+kmin*csmax)*(ii+kmin*csmax) ...
        +(jj+kmin*ssmax)*(jj+kmin*ssmax)-kmin*kmin;
   d2 = (ii+kmax*csmin)*(ii+kmax*csmin) ...
        +(jj+kmax*ssmin)*(jj+kmax*ssmin)-kmax*kmax;
   d3 = (ii-kmin*crmin)*(ii-kmin*crmin) ...
        +(jj-kmin*srmin)*(jj-kmin*srmin)-kmin*kmin;
   d4 = (ii-kmax*crmax)*(ii-kmax*crmax) ...
        +(jj-kmax*srmax)*(jj-kmax*srmax)-kmax*kmax;
   d5 = (ii+kmin*crmax)*(ii+kmin*crmax) ...
        +(jj+kmin*srmax)*(jj+kmin*srmax)-kmin*kmin;
   d6 = (ii+kmax*crmin)*(ii+kmax*crmin) ...
        +(jj+kmax*srmin)*(jj+kmax*srmin)-kmax*kmax;
   d7 = (ii-kmin*csmin)*(ii-kmin*csmin) ...
        +(jj-kmin*ssmin)*(jj-kmin*ssmin)-kmin*kmin;
   d8 = (ii-kmax*csmax)*(ii-kmax*csmax) ...
        +(jj-kmax*ssmax)*(jj-kmax*ssmax)-kmax*kmax;
   if ((d2<=0)&(d4<=0)&(d1>=0)&(d3>=0))| ...
      ((d6<=0)&(d8<=0)&(d5>=0)&(d7>=0))
      matf(i,j) = 1;
   end
% end of the three tests
end
end
end
% end of the big loop
end
end

% symmetry
for i=round(nxmax/2)+1:nxmax;
    for j=1:nmax;
        matf(i,j) = matf(nxmax-i+1,nmax-j+1);
    end
end
end


% we arbitrarely put 1 at the origin
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
end
