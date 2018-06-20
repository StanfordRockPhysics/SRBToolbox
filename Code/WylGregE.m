function e = WylGregE(Phi,d) 
%function e = WylGregE(Phi,d) 
%
%This equation relates the permeability of a porous
%medium to its specific surface area (S)
%and porosity (Phi).
% 
%   K= (Phi.^3)./(C.*S.^2);
%
% substitute S with its expression 6.*(1-Phi)./d:
%
%  K = (1/(36*C)).*(d.^2).*(Phi.^3)/(1-Phi).^2;
%
%where:Phi is the porosity, scalar or vector,
%      C is the Kozeny-Carman constant,
%      d is pore diameter in micron-meter, and
%      S is the surface area per unit volume.
%
%The output is two column matrix, porosity in fraction
%in first column, and permeability (in md) in second
%column.
%
%VALIDITY:The equation is applicable only under 
%         conditions of viscous flow. 
%REFERENCE:"Fluid Flow through Unconsolidated Porous 
%          Aggregates," M.R.J. Wyllie and A. R. Gregory,
%          Industrial and Engineering Chemistry, 1997.
%
% See also BERNABE, BLOCH, COATES, COATDUM, FREDRICH, 
%          KOZCARM, MODKOZCARM, PANDALAKE,PANDALAKEKC,
%          TIMUR, TXIER, OWOLABI, REVIL 

% Written by R. E. Rasolovoahangy, June 2000

if nargin==0   
   prompt={'Porosity (Phi)','Pore Diameter (d), micronmeter'};
   defans={'[0.01:0.01:0.35]','250'};
   getpar=inputdlg(prompt,'Wyllie-Gregory Model Default Input Parameters',1,defans);
   if isempty(getpar)
	return;
   end
   for k=1:length(getpar),
       if isempty(str2num(getpar{k}))==1
          param{k} = evalin('base',getpar{k});
       else
          param{k} = str2num(getpar{k});
       end
   end
   e = WG(param{1},param{2});  
else
e = WG(Phi,d);
end
   
function e = WG(Phi,d)

K = 1.852.*(d.^2).*(Phi.^3)./(1-Phi).^2;

semilogy(Phi,K,'b-');
set(gca,'fontsize',9);
xlabel('Porosity','fontsize',9);
ylabel('Permeability (md)','fontsize',9);
title('Wyllie-Gregory Permeability versus Porosity','fontsize',10);
grid on;
hold on;
e = [Phi K];

