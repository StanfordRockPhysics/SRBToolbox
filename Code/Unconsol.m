function result = Unconsol(varargin)
%function result = Unconsol(PhiC,C,Gs,Nus,Kf,P)
%Modified Hashin-Shtrikman lower bound with Hertz-Mindlin end point,
% for Marine Sediment Skeleton. Used for modeling unconsolidated sediments.   
%Entire porosity range.  Sediment frame fully water saturated
%Inputs by dialog box, if called without any input arguments.
%PhiC		Porosity of a sphere pack  
%C		Coordination number 
%Gs		Solid shear modulus (GPa)  
%Nus		Solid Poisson's ratio 
%Kf		Fluid bulk modulus (GPa)  
%P		Effective pressure (MPa)
%output result matrix: result=[porosity M-modulus shear-modulus];
%Plots porosity vs. modluli if called with no output arguments;

%Written by Jack Dvorkin
%I/O modifications T. Mukerji

prompt={'PhiC','Coord.#','Gs (GPa)','Nus','Kf (GPa)','P (MPa)'};
defans={'.38','8.5','6.58','.352','2.25','20'};

if nargin==0
getpar=inputdlg(prompt,'Unconsolidated Model',1,defans);
for k=1:length(getpar), param(k)=str2num(getpar{k}); end;
PhiC=param(1); C=param(2); G=param(3); nu=param(4); Kf=param(5); P=param(6);
else
PhiC=varargin{1};C=varargin{2};G=varargin{3};nu=varargin{4};
Kf=varargin{5}; P=varargin{6};
end;
format short
P = P./1000;
Phi0=PhiC;

%Porosity loop 
M = 2.*G.*(1-nu)./(1-2.*nu);
K = M-4.*G./3;
i = (1:100)';
Phi = 0.01.*i;
%Effective K and G at Phi0     
b = (3.*3.14.*(1.-nu).*P./(2.*C.*(1.-Phi0).*G)).^(1./3);
Keff = C.*(1.-Phi0).*G.*b./(3.*3.14.*(1.-nu));
Geff = C.*(1.-Phi0).*G.*b.*(5.-4.*nu)./(5.*3.14.*(1.-nu).*(2.-nu));
Khat = Keff; Ghat = Geff; 
%Effective bulk and shear moduli at porosity Phi(i)   

Keff = 1./((Phi./Phi0)./(Khat+4.*Ghat./3)+((Phi0-Phi)./Phi0)./(K+4.*Ghat./3))-4.*Ghat./3;
ZZ = (Ghat./6).*(9.*Khat+8.*Ghat)./(Khat+2.*Ghat);
Geff = 1./((Phi./Phi0)./(Ghat+ZZ)+((Phi0-Phi)./Phi0)./(G+ZZ))-ZZ; 

%Gassmann water 
KgassmW = K.*(Phi.*Keff-(1+Phi).*Kf.*Keff./K+Kf)./((1-Phi).*Kf+Phi.*K-Kf.*Keff./K);    
Mdry = Keff+(4./3).*Geff; 
Meff = KgassmW+(4./3).*Geff;    

if nargout==0
subplot(1,2,1)
plot(Phi,Meff,'blue-')
axis([0.1 0.4 0 45])
set(gca,'fontname','bookman','fontsize',9)
xlabel('Porosity','fontname','bookman','fontsize',11)
ylabel('M-Modulus (GPa)','fontname','bookman','fontsize',11)
hold on
subplot(1,2,2)
plot(Phi,Geff,'blue-')
axis([0.1 0.4 0 20])
set(gca,'fontname','bookman','fontsize',9)
xlabel('Porosity','fontname','bookman','fontsize',11)
ylabel('G-Modulus (GPa)','fontname','bookman','fontsize',11)
hold on
end;

result=[Phi Meff Geff];

