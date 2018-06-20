function [matt] = f23(dx1,dz1,nxmax,nmax,kmax,kmin,dsmax,dsmin,drmax,drmin,dr)
% function used by bornfilt. Usually need not be called seperately by user.
% see BORNFILT


% finite lines of sources and receivers

% WTW reflection

matt = zeros(nxmax,nmax);
matd = f21(dz1,dx1,nmax,nxmax,kmax,kmin,-dsmin,-dsmax,-drmin,-drmax,dr);
for i = 1:nxmax;
   for j= 1:nmax;
       matt(i,j) = matd(nmax+1-j,i);
   end
end
clear matd;
end
