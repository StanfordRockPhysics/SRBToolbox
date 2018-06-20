function e = Tixier(Phi,Swr) 
% function e = Tixier(Phi,Swr)
%
% Timur relation provides a way to estimate the 
% permeability of sandstone from in situ measurements 
% of porosity(Phi)and residual fluid saturation (Swr). 
%
%  K=(A*Phi^B)/Swr^C
% where: K is the Permeability,
%        Phi is the porosity, scalar or vector, and 
%        Swr is the residual water saturation, scalar or vector.
%
% Tixier, Coates-Dumanoir, Coates, and Timur have the same form
% as above, but the coefficient A, B, C are different in each case.
%
% The output is two column matrix, porosity (in fraction) 
% in first column, and permeability (in md) in second column.
%
% ASSUMPTIONS: A,B,C are determined empirically 
% VALIDITY: Unconsolidated Sands
% REFERENCE:"Evaluation of Permeability from Eletrical-Log,"
%            Tixier et al., 1994, OGJ Vol. 48, n.6, p.113-123.
%           "An Empirical Expression for Permeability
%            in Unconsolidated Sands of Eastern Niger
%            Delta,"Journal Of petroleum Geology, 
%            Vol. 17 (1),January 1994.
%
% See also BERNABE, BLOCH, COATES, COATDUM, FREDRICH, 
%          KOZCARM, MODKOZCARM, PANDALAKE,PANDALAKEKC,
%          TIMUR, OWOLABI, REVIL, WYLGREG 

% Written by R. E. Rasolovoahangy, June 2000

if nargin==0   
   prompt={'Porosity (Phi)','Irreducible Water Saturation (Swr)'};
   defans={'[0.01:0.01:0.35]','0.15'};
   getpar=inputdlg(prompt,'Tixier Model Default Input Parameters',1,defans);
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
   e = Tix(param{1},param{2});  
else
e = Tix(Phi,Swr);
end
   
function e = Tix(Phi,Swr)    

K = 62500.*((Phi.^6)./(Swr.^2));

semilogy(Phi,K,'g-');
set(gca,'fontsize',9);
xlabel('Porosity','fontsize',9);
ylabel('Permeability (md)','fontsize',9);
title('Tixier Permeability versus Porosity','fontsize',10);
grid on;
hold on;
e = [Phi K];

