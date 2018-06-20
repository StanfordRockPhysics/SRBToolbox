function [A,B,C,D]=TransConfig2(f,Gam_c,Gam_ml,Cprime,Co,phi,...
    d_c,d_oc,t_ml,Z_c,Z_bl,Z_ml)

j = sqrt(-1);
arg_c = Gam_c*d_c/2;
arg_oc = Gam_c*d_oc;
arg_ml = Gam_ml*t_ml;

%Pre-allocate variable sizes
M_Cprime=zeros(2,2,length(f)); M_Co=zeros(2,2,length(f));
M_phi=zeros(2,2,length(f)); M_f=zeros(2,2,length(f));
Z_bocf=length(f); M_ocf=zeros(2,2,length(f));
M_oc=zeros(2,2,length(f)); M_ml=zeros(2,2,length(f));
M_bocf=zeros(2,2,length(f)); M=zeros(2,2,length(f));

for i=1:length(f)
    M_Cprime(:,:,i) = [1 -j./(2*pi*f(i).*Cprime(i))
        0 1];
    M_Co(:,:,i) = [1 -j./(2*pi*f(i).*Co)
        0 1];
    M_phi(:,:,i) = [phi(i) 0
        0 1./phi(i)];

    M_f(:,:,i) = [cosh(arg_c(i)) Z_c*sinh(arg_c(i))
        1./Z_c.*sinh(arg_c(i)) cosh(arg_c(i))];

    M_oc(:,:,i) = [cosh(arg_oc(i)) Z_c*sinh(arg_oc(i))
        1/Z_c*sinh(arg_oc(i)) cosh(arg_oc(i))];
    M_ocf(:,:,i) = squeeze(M_f(:,:,i))*squeeze(M_oc(:,:,i));

    Z_bocf(i) = (M_ocf(1,1,i)*Z_bl+M_ocf(1,2,i))/...
        (M_ocf(2,1,i)*Z_bl+M_ocf(2,2,i));
    M_bocf(:,:,i) = [1 0
        1./Z_bocf(i) 1];
    
    M_ml(:,:,i) = [cosh(arg_ml(i)) Z_ml*sinh(arg_ml(i))
        1./Z_ml*sinh(arg_ml(i)) cosh(arg_ml(i))];
    
    M(:,:,i) = squeeze(M_Cprime(:,:,i))*squeeze(M_Co(:,:,i))...
        *squeeze(M_phi(:,:,i))*squeeze(M_bocf(:,:,i))*...
        squeeze(M_f(:,:,i))*squeeze(M_ml(:,:,i));
    
    A = squeeze(M(1,1,:));
    B = squeeze(M(1,2,:));
    C = squeeze(M(2,1,:));
    D = squeeze(M(2,2,:));
end
