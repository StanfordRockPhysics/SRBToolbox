function[matf]=f22(dx1,dz1,nxmax,nmax,kmax,kmin,dsmax,dsmin,drmax,drmin,dsf,dr)
% function used by bornfilt. Usually need not be called seperately by user.
% see BORNFILT


% finite lines of sources and receivers

% VSP

im = sqrt(-1);
dx = 2*pi*(nmax-1)/(nmax*dx1);
dz = 2*pi*(nxmax-1)/(nxmax*dz1);

% angles
d = sqrt(drmax*drmax+dr*dr);
srmax = drmax/d;
crmax = sqrt(1-srmax*srmax);
rmax = asin(drmax/d);
d = sqrt(drmin*drmin+dr*dr);
rmin = asin(drmin/d);
srmin = drmin/d;
crmin = sqrt(1-srmin*srmin);
d = sqrt(dsmax*dsmax+dsf*dsf);
csmin = -dsmax/d;
smin = -pi+acos(dsmax/d);
ssmin = -sqrt(1-csmin*csmin);
d = sqrt(dsmin*dsmin+dsf*dsf);
csmax = -dsmin/d;
smax = -pi+acos(dsmin/d);
ssmax = -sqrt(1-csmax*csmax);
xd = crmin-csmin;
yd = srmin-ssmin;
xg = crmax-csmax;
yg = srmax-ssmax;

matf = zeros(nxmax,nmax);

% first case
if ((smax+pi)>=rmax)&(rmax>=(smin+pi))&((smin+pi)>=rmin)
for i=1:nxmax;
c = round((nmax+1)/2-dz*(nmax-1)*(1+nxmax-2*i)/(2*dx*(nxmax-1)));
if c<1
   c = 1;
elseif c>nmax
   c = nmax;
end
   for j=c:nmax;
% matrix (i,j) => spatial coordinates (ii,jj)
ii=(j-(nmax+1)/2)*dx;
jj=(-i+(nxmax+1)/2)*dz;
% first test
if (ii*ii+jj*jj) <= (4*kmax*kmax)
% first case
if (smax-rmin) >= 0
   d2 = (ii+kmax*crmin)*(ii+kmax*crmin) ...
        +(jj+kmax*srmin)*(jj+kmax*srmin)-kmax*kmax;
   d1 = (ii-kmax*csmax)*(ii-kmax*csmax) ...
        +(jj-kmax*ssmax)*(jj-kmax*ssmax)-kmax*kmax;
   if (d2<=0)&(d1<=0)
      matf(i,j) = 1;
   end
end
% second and third test
if (ii*yd-jj*xd) <= 0
if (ii*yg-jj*xg) >= 0
   d2 = (ii+kmax*csmin)*(ii+kmax*csmin) ...
        +(jj+kmax*ssmin)*(jj+kmax*ssmin)-kmax*kmax;
   d3 = (ii-kmin*crmin)*(ii-kmin*crmin) ...
        +(jj-kmin*srmin)*(jj-kmin*srmin)-kmin*kmin;
   d4 = (ii-kmax*crmax)*(ii-kmax*crmax) ...
        +(jj-kmax*srmax)*(jj-kmax*srmax)-kmax*kmax;
   d1 = (ii+kmin*csmax)*(ii+kmin*csmax) ...
        +(jj+kmin*ssmax)*(jj+kmin*ssmax)-kmin*kmin;
   if ((d1>=0)&(d3>=0)&((d2<=0)|(d4<=0)))
      matf(i,j) = 1;
   else
   theta = phase(ii+im*jj);
   if (d4>=0)&(d2>=0)&(theta>=(pi+smin))&(theta<=(rmax))
      matf(i,j) = 1;
   end
   end
% end of the 3 tests
end
end
end
end
end

% second case
elseif ((smax+pi)<=rmax)&((smin+pi)>=rmin)
for i=1:nxmax;
c = round((nmax+1)/2-dz*(nmax-1)*(1+nxmax-2*i)/(2*dx*(nxmax-1)));
if c<1
   c = 1;
end
   for j=c:nmax;
% matrix (i,j) => spatial coordinates (ii,jj)
ii=(j-(nmax+1)/2)*dx;
jj=(-i+(nxmax+1)/2)*dz;
% first test
if (ii*ii+jj*jj) <= (4*kmax*kmax)
% second et third test
if (ii*yd-jj*xd) <= 0
if (ii*yg-jj*xg) >= 0
   d2 = (ii+kmax*csmin)*(ii+kmax*csmin) ...
        +(jj+kmax*ssmin)*(jj+kmax*ssmin)-kmax*kmax;
   d3 = (ii-kmin*crmin)*(ii-kmin*crmin) ...
        +(jj-kmin*srmin)*(jj-kmin*srmin)-kmin*kmin;
   d6 = (ii+kmin*csmin)*(ii+kmin*csmin) ...
        +(jj+kmin*ssmin)*(jj+kmin*ssmin)-kmin*kmin;
   d4 = (ii+kmax*csmax)*(ii+kmax*csmax) ...
        +(jj+kmax*ssmax)*(jj+kmax*ssmax)-kmax*kmax;
   d1 = (ii+kmin*csmax)*(ii+kmin*csmax) ...
        +(jj+kmin*ssmax)*(jj+kmin*ssmax)-kmin*kmin;
   d5 = (ii-kmin*crmax)*(ii-kmin*crmax) ...
        +(jj-kmin*srmax)*(jj-kmin*srmax)-kmin*kmin;
   if ((d5>=0)&(d6>=0)&(d4<=0))|((d3>=0)&(d1>=0)&(d2<=0))
      matf(i,j) = 1;
   else
   theta = phase(ii+im*jj);
   if (d4>=0)&(d2>=0)&(theta>=(pi+smin))&(theta<=(pi+smax))
      matf(i,j) = 1;
   end
   end
% end of the 3 tests
end
end
end
end
end

% third case
elseif ((smax+pi)>=rmax)&((smin+pi)<=rmin)
for i=1:nxmax;
c = round((nmax+1)/2-dz*(nmax-1)*(1+nxmax-2*i)/(2*dx*(nxmax-1)));
if c<1
   c = 1;
end
   for j=c:nmax;
% matrix (i,j) => spatial coordinates (ii,jj)
ii=(j-(nmax+1)/2)*dx;
jj=(-i+(nxmax+1)/2)*dz;
% first test
if (ii*ii+jj*jj) <= (4*kmax*kmax)
% second and third test
if (ii*yd-jj*xd) <= 0
if (ii*yg-jj*xg) >= 0
   d2 = (ii-kmax*crmin)*(ii-kmax*crmin) ...
        +(jj-kmax*srmin)*(jj-kmax*srmin)-kmax*kmax;
   d3 = (ii-kmin*crmin)*(ii-kmin*crmin) ...
        +(jj-kmin*srmin)*(jj-kmin*srmin)-kmin*kmin;
   d6 = (ii+kmin*csmin)*(ii+kmin*csmin) ...
        +(jj+kmin*ssmin)*(jj+kmin*ssmin)-kmin*kmin;
   d4 = (ii-kmax*crmax)*(ii-kmax*crmax) ...
        +(jj-kmax*srmax)*(jj-kmax*srmax)-kmax*kmax;
   d1 = (ii+kmin*csmax)*(ii+kmin*csmax) ...
        +(jj+kmin*ssmax)*(jj+kmin*ssmax)-kmin*kmin;
   d5 = (ii-kmin*crmax)*(ii-kmin*crmax) ...
        +(jj-kmin*srmax)*(jj-kmin*srmax)-kmin*kmin;
   if ((d5>=0)&(d6>=0)&(d2<=0))|((d3>=0)&(d1>=0)&(d4<=0))
      matf(i,j) = 1;
   else
   theta = phase(ii+im*jj);
   if (d4>=0)&(d2>=0)&(theta>=rmin)&(theta<=rmax)
      matf(i,j) = 1;
   end
   end
% end of the 3 tests
end
end
end
end
end

% fourth case
elseif ((smax+pi)<=rmax)&(rmin<=(smax+pi))&((smin+pi)<=rmin)
for i=1:nxmax;
c = round((nmax+1)/2-dz*(nmax-1)*(1+nxmax-2*i)/(2*dx*(nxmax-1)));
if c<1
   c = 1;
end
   for j=c:nmax;
% matrix (i,j) => spatial coordinates (ii,jj)
ii=(j-(nmax+1)/2)*dx;
jj=(-i+(nxmax+1)/2)*dz;
% first test
if (ii*ii+jj*jj) <= (4*kmax*kmax)
% second and third test
if (ii*yd-jj*xd) <= 0
if (ii*yg-jj*xg) >= 0
   d2 = (ii-kmax*crmin)*(ii-kmax*crmin) ...
        +(jj-kmax*srmin)*(jj-kmax*srmin)-kmax*kmax;
   d6 = (ii+kmin*csmin)*(ii+kmin*csmin) ...
        +(jj+kmin*ssmin)*(jj+kmin*ssmin)-kmin*kmin;
   d4 = (ii+kmax*csmax)*(ii+kmax*csmax) ...
        +(jj+kmax*ssmax)*(jj+kmax*ssmax)-kmax*kmax;
   d5 = (ii-kmin*crmax)*(ii-kmin*crmax) ...
        +(jj-kmin*srmax)*(jj-kmin*srmax)-kmin*kmin;
   if (d5>=0)&(d6>=0)&((d2<=0)|(d4<=0))
      matf(i,j) = 1;
   else
   theta = phase(ii+im*jj);
   if (d4>=0)&(d2>=0)&(theta>=rmin)&(theta<=(pi+smax))
      matf(i,j) = 1;
   end
   end
% end of the 3 tests
end
end
end
end
end

% fifth case
elseif ((smax+pi)<=rmax)&(rmin>=(smax+pi))&((smin+pi)<=rmin)
for i=1:nxmax;
c = round((nmax+1)/2-dz*(nmax-1)*(1+nxmax-2*i)/(2*dx*(nxmax-1)));
if c<1
   c = 1;
end
   for j=c:nmax;
% matrix (i,j) => spatial coordinates (ii,jj)
ii=(j-(nmax+1)/2)*dx;
jj=(-i+(nxmax+1)/2)*dz;
% first test
if (ii*ii+jj*jj) <= (4*kmax*kmax)
% second and third test
if (ii*yd-jj*xd) <= 0
if (ii*yg-jj*xg) >= 0
   d2 = (ii-kmax*crmin)*(ii-kmax*crmin) ...
        +(jj-kmax*srmin)*(jj-kmax*srmin)-kmax*kmax;
   d6 = (ii+kmin*csmin)*(ii+kmin*csmin) ...
        +(jj+kmin*ssmin)*(jj+kmin*ssmin)-kmin*kmin;
   d4 = (ii+kmax*csmax)*(ii+kmax*csmax) ...
        +(jj+kmax*ssmax)*(jj+kmax*ssmax)-kmax*kmax;
   d5 = (ii-kmin*crmax)*(ii-kmin*crmax) ...
        +(jj-kmin*srmax)*(jj-kmin*srmax)-kmin*kmin;
   if (d5>=0)&(d6>=0)&(d2<=0)&(d4<=0)
      matf(i,j) = 1;
   end
% end of the 3 tests
end
end
end
end
end

% sixth case
elseif ((smax+pi)>=rmax)&(rmax<=(smin+pi))
for i=1:nxmax;
c = round((nmax+1)/2-dz*(nmax-1)*(1+nxmax-2*i)/(2*dx*(nxmax-1)));
if c<1
   c = 1;
end
for j=c:nmax;
% matrix (i,j) => spatial coordinates (ii,jj)
ii=(j-(nmax+1)/2)*dx;
jj=(-i+(nxmax+1)/2)*dz;
% first test
if (ii*ii+jj*jj) <= (4*kmax*kmax)
% 4 cases
if (xd>=0)&(yg>=0)
   if (ii*yd-jj*xd) <= 0
   if (ii*yg-jj*xg) >= 0
   d2 = (ii-kmax*crmax)*(ii-kmax*crmax) ...
        +(jj-kmax*srmax)*(jj-kmax*srmax)-kmax*kmax;
   d6 = (ii+kmin*csmax)*(ii+kmin*csmax) ...
        +(jj+kmin*ssmax)*(jj+kmin*ssmax)-kmin*kmin;
   d4 = (ii+kmax*csmin)*(ii+kmax*csmin) ...
        +(jj+kmax*ssmin)*(jj+kmax*ssmin)-kmax*kmax;
   d5 = (ii-kmin*crmin)*(ii-kmin*crmin) ...
        +(jj-kmin*srmin)*(jj-kmin*srmin)-kmin*kmin;
   if (d2<=0)&(d4<=0)&(d6>=0)&(d5>=0)
      matf(i,j) = 1;
   end
   end
   end 
elseif (xd<=0)&(yg>=0)
   if (ii*yd-jj*xd) >= 0
   if (ii*yg-jj*xg) >= 0
   d1 = (ii-kmax*crmax)*(ii-kmax*crmax) ...
        +(jj-kmax*srmax)*(jj-kmax*srmax)-kmax*kmax;
   d2 = (ii+kmin*csmax)*(ii+kmin*csmax) ...
        +(jj+kmin*ssmax)*(jj+kmin*ssmax)-kmin*kmin;
   d3 = (ii+kmax*csmax)*(ii+kmax*csmin) ...
        +(jj+kmax*ssmin)*(jj+kmax*ssmin)-kmax*kmax;
   d4 = (ii+kmax*crmin)*(ii+kmax*crmin) ...
        +(jj+kmax*srmin)*(jj+kmax*srmin)-kmax*kmax;
   d5 = (ii-kmin*csmin)*(ii-kmin*csmin) ...
        +(jj-kmin*ssmin)*(jj-kmin*ssmin)-kmin*kmin;
   d6 = (ii-kmax*csmax)*(ii-kmax*csmax) ...
        +(jj-kmax*ssmax)*(jj-kmax*ssmax)-kmax*kmax;
   if ((d1<=0)&(d3<=0)&(d2>=0))|((d6<=0)&(d4<=0)&(d5>=0)) 
       matf(i,j) = 1;
   end
   end
   end
elseif (xd>=0)&(yg<=0)
   if (ii*yd-jj*xd) <= 0
   if (ii*yg-jj*xg) <= 0
   d1 = (ii-kmax*crmax)*(ii-kmax*crmax) ...
        +(jj-kmax*srmax)*(jj-kmax*srmax)-kmax*kmax;
   d2 = (ii-kmin*crmin)*(ii-kmin*crmin) ...
        +(jj-kmin*srmin)*(jj-kmin*srmin)-kmin*kmin;
   d3 = (ii+kmax*csmax)*(ii+kmax*csmin) ...
        +(jj+kmax*ssmin)*(jj+kmax*ssmin)-kmax*kmax;
   d4 = (ii+kmax*crmin)*(ii+kmax*crmin) ...
        +(jj+kmax*srmin)*(jj+kmax*srmin)-kmax*kmax;
   d5 = (ii+kmin*crmax)*(ii+kmin*crmax) ...
        +(jj+kmin*srmax)*(jj+kmin*srmax)-kmin*kmin;
   d6 = (ii-kmax*csmax)*(ii-kmax*csmax) ...
        +(jj-kmax*ssmax)*(jj-kmax*ssmax)-kmax*kmax;
   if ((d1<=0)&(d3<=0)&(d2>=0))|((d6<=0)&(d4<=0)&(d5>=0)) 
       matf(i,j) = 1;
   end
   end
   end
elseif (xd<=0)&(yg<=0)
   if (ii*yd-jj*xd) >= 0
   if (ii*yg-jj*xg) <= 0
   d4 = (ii+kmax*crmin)*(ii+kmax*crmin) ...
        +(jj+kmax*srmin)*(jj+kmax*srmin)-kmax*kmax;
   d5 = (ii+kmin*crmax)*(ii+kmin*crmax) ...
        +(jj+kmin*srmax)*(jj+kmin*srmax)-kmin*kmin;
   d3 = (ii-kmin*csmin)*(ii-kmin*csmin) ...
        +(jj-kmin*ssmin)*(jj-kmin*ssmin)-kmin*kmin;
   d6 = (ii-kmax*csmax)*(ii-kmax*csmax) ...
        +(jj-kmax*ssmax)*(jj-kmax*ssmax)-kmax*kmax;
   if (d6<=0)&(d4<=0)&(d5>=0)&(d3>=0) 
       matf(i,j) = 1;
   end
   end
   end
else
   disp('erreur')
end
% end of test
end
end
end
end

% symmetry
for i=1:nxmax;
c = round((nmax+1)/2-dz*(nmax-1)*(1+nxmax-2*i)/(2*dx*(nxmax-1)));
if c>nmax
   c = nmax;
end
    for j=1:c-1;
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
