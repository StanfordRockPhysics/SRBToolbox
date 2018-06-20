function [A,B,C,D]=TransConfig1(f,Gam_c,Gam_ml,Cprime,Co,phi,...
    d_c,d_oc,t_ml,Z_c,Z_bl,Z_ml)

j = sqrt(-1);
arg_c = Gam_c*d_c/2;
arg_oc = Gam_c*d_oc;
arg_ml = Gam_ml*t_ml;

%Pre-allocate variable sizes
M_Cprime=zeros(2,2,length(f)); M_Co=zeros(2,2,length(f));
M_phi=zeros(2,2,length(f)); M_f=zeros(2,2,length(f));
Z_ba=length(f); M_ba=zeros(2,2,length(f));
M_oc=zeros(2,2,length(f)); M_ml=zeros(2,2,length(f));
M=zeros(2,2,length(f));

for i=1:length(f)
    M_Cprime(:,:,i) = [1 -j./(2*pi*f(i).*Cprime(i))
        0 1];
    M_Co(:,:,i) = [1 -j./(2*pi*f(i).*Co)
        0 1];
    M_phi(:,:,i) = [phi(i) 0
        0 1./phi(i)];

    M_f(:,:,i) = [cosh(arg_c(i)) Z_c*sinh(arg_c(i))
        1./Z_c.*sinh(arg_c(i)) cosh(arg_c(i))];

    Z_ba(i) = (M_f(1,1,i)*Z_bl+M_f(1,2,i))/...
        (M_f(2,1,i)*Z_bl+M_f(2,2,i));
    M_ba(:,:,i) = [1 0
        1./Z_ba(i) 1];
    
    M_oc(:,:,i) = [cosh(arg_oc(i)) Z_c*sinh(arg_oc(i))
        1/Z_c*sinh(arg_oc(i)) cosh(arg_oc(i))];
    
    M_ml(:,:,i) = [cosh(arg_ml(i)) Z_ml*sinh(arg_ml(i))
        1./Z_ml*sinh(arg_ml(i)) cosh(arg_ml(i))];
    
    M(:,:,i) = squeeze(M_Cprime(:,:,i))*squeeze(M_Co(:,:,i))...
        *squeeze(M_phi(:,:,i))*squeeze(M_ba(:,:,i))*...
        squeeze(M_f(:,:,i))*squeeze(M_ml(:,:,i));
    
    A = squeeze(M(1,1,:));
    B = squeeze(M(1,2,:));
    C = squeeze(M(2,1,:));
    D = squeeze(M(2,2,:));
end
