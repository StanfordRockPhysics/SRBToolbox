function zz=nmo0(Tr,V,t0,dt,x,flag);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function zz=nmo0(Tr,V,t0,dt,x,flag);
%
%This takes a trace and NMO correct it or spray a zero offset trace 
%into an hyperbola
%
% Tr:  Trace at distance x
% t0- two way travel time at x=0
% V- velocity (size(Tr))
% dt- Sampling rate
% zz-  nmo'ed trace
% flag- 0 for nmo correction (make hyperbola flat), 1 for spraying trace 
% Note:  Due to NMO stretch the amplitude are not consistent from
%input/output
% Best to use agc after nmo correction.  See Jon's BEI for details

%written by Ran Bachrach, 1999
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=length(Tr);
   zz=Tr*0;
for iz=1:n
   z=t0+dt*(iz-1);
   t=sqrt(z^2+(x/V(iz))^2);
   it=1+round((t-t0)/dt);% Round to nearest neighbor
   if it<=n
      if flag==1
         zz(it)=zz(it)+Tr(iz);
         else
        zz(iz)=zz(iz)+Tr(it);
        end
   end
end


