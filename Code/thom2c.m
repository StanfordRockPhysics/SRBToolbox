function c = thom2c(ex,gx,dx,ey,gy,dy,gxy,d3,vp,vs,rho);
% c = thom2c(ex,gx,dx,ey,gy,dy,gxy,d3,vp,vs,rho);
%
% converts Tsvankin's extended Thomsen parameters to 6x6 elastic stiffness tensor.
% Assumes an orthorhombic symmetry or higher.
%
% inputs
%   ex,gx,dx,   Tsvankin's extended Thomsen parameters in y-z symmetry planes
%   ey,gy,dy,   Tsvankin's extended Thomsen parameters in x-z symmetry planes
%   gxy         Tsvankin's extended Thomsen parameters parameters in xy symmetry planes
%   d3          extra parameter defined in the xy plane
%   vp,vs       vertical velocities sqrt(c33/rho);sqrt(c44/rho);
%   rho     	density, used only to compute velocities
% outputs
%   c           elastic 6x6 stiffness tensor (Voigt notation)

% Gary Mavko, August 2003

c33 = rho.*vp.^2;
c11 = c33 + 2*c33.*ey;
c22 = c33 + 2*c33.*ex;
c44 = rho.*vs.^2;
c55 = c44./(2*gxy+1);
c66 = 2*c44.*gy + c44;
c13 = -c55 + sqrt(dy.*(2*c33.*(c33-c55)) + (c33-c55).^2);
c23 = -c44 + sqrt(dx.*(2*c33.*(c33-c44)) + (c33-c44).^2);
c12 = -c66 + sqrt(d3.*(2*c11.*(c11-c66)) + (c11-c66).^2);

c = zeros(6,6);
c(1,1,:) = c11; 
c(2,2,:) = c22;
c(3,3,:) = c33;
c(4,4,:) = c44;
c(5,5,:) = c55;
c(6,6,:) = c66;
c(1,2,:) = c12;
c(2,1,:) = c12;
c(1,3,:) = c13;
c(3,1,:) = c13;
c(2,3,:) = c23;
c(3,2,:) = c23;

c = squeeze(c);