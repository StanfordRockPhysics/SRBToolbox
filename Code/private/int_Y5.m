function Y5_integration=int_Y5(theta,phi,E)
I(1,:)=sin(theta).*cos(phi);
I(2,:)=sin(theta).*sin(phi);
I(3,:)=cos(theta);
Y5=sqrt(E(1,1).*I(1,:).*I(1,:) + E(1,2).*I(1,:).*I(2,:) + E(1,3).*I(1,:).*I(3,:) ...
    + E(2,1).*I(2,:).*I(1,:) + E(2,2).*I(2,:).*I(2,:) + E(2,3).*I(2,:).*I(3,:) ...
    + E(3,1).*I(3,:).*I(1,:) + E(3,2).*I(3,:).*I(2,:) + E(3,3).*I(3,:).*I(3,:));
Y5_integration = (1/(4*pi)).*sin(theta).*Y5.*2.*(-E(3,1).*I(1,:).*I(3,:) -E(3,2).*I(2,:).*I(3,:) -E(3,3).*I(3,:).*I(3,:));

% test=dblquad(@(theta,phi)int_X_I1sq(theta,phi,E),0,pi,0,2*pi,1,@quadl);
% test=dblquad(@(theta,phi)int_X_I1sq(theta,phi,E),0,pi,0,2*pi);