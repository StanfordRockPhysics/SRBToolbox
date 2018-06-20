function [matf] = f25(dx1,dz1,nxmax,nzmax,kmax,kmin,dsmax,dsmin,drmax,drmin)
% function used by bornfilt. Usually need not be called seperately by user.
% see BORNFILT


% dx1,dz1 : size of the picture
% nxmax,nzmax : dimensions of the matrix

% log

dz = 2*pi*(nxmax-1)/(nxmax*dz1);
matf = zeros(nxmax,nzmax);
if ((dsmin<0)&(drmin<0))|((drmax>0)&(dsmax>0))
if rem(nzmax,2) == 0
   for i=1:nxmax
       ii=abs((-i+(nxmax+1)/2)*dz);
       if ii <= 2*kmax
          if ii >= 2*kmin
             matf(i,nzmax/2) = 1;
             matf(i,nzmax/2+1) = 1;
          end
       end
   end
% we put 1 at the center of the filter
   if rem(nxmax,2) == 0
      matf(nxmax/2,nzmax/2+1) = 1;
      matf(nxmax/2,nzmax/2) = 1;
      matf(nxmax/2+1,nzmax/2+1) = 1;
      matf(nxmax/2+1,nzmax/2) = 1;
   else
      matf((1+nxmax)/2,(nzmax)/2+1) = 1;
      matf((1+nxmax)/2,(nzmax)/2) = 1;
   end
else
   for i=1:nxmax
       ii=abs((-i+(nxmax+1)/2)*dz);
       if ii <= 2*kmax
          if ii >= 2*kmin
             matf(i,(nzmax+1)/2) = 1;
          end
       end
   end
% we put 1 at the center of the filter
   if rem(nxmax,2) == 0
      matf(nxmax/2+1,(nzmax+1)/2) = 1;
      matf(nxmax/2,(nzmax+1)/2) = 1;
   else
      matf((nxmax+1)/2,(1+nzmax)/2) = 1;
   end
end
else
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
end
