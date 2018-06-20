function e = Timur(Phi,Swr) 
% function e=Timur(Phi,Swr) 
%
% Timur relation provides a way to estimate the 
% permeability of sandstone from in situ measurements 
% of porosity(Phi)and residual fluid saturation (Swr). 
%
%  K=(A*Phi^B)/Swr^C
% where: K is the Permeability,
%        Phi is the porosity, scalar or vector, and 
%        Swr is the irreducible water saturation, scalar or vector.
%
% Timur, Tixier, Coates-Dumanoir, and Coates have the same form
% as above, but the coefficient A, B, C are different in each case.
%
% The output is two column matrix, porosity in fraction in first 
% column, and permeability (in md) in second column.
%
% ASSUMPTIONS: A,B,C are determined empirically 
% LIMITATIONS:Homogeneous and uniform reservoir; Poorly  
%             consolidated; well sorted with relative high 
%             Phi sandstones of the Gulf Coast; well-consolidated,
%             tight sandstones from Colorado field; moderately
%             well-consolidated but poorly sorted micaceous
%             sandstones from California Field.
% REFERENCE: "An investigation of Permeability, Porosity, and
%             Residual Water Saturation Relationship for 
%             Sandstone Reservoirs," A. Timur, 9th Annual 
%             SPWLA Logging Symposium.
%
% See also BERNABE, BLOCH,COATDUM, COATES, FREDRICH, KOZCARM, 
%          MODKOZCARM, PANDALAKE,PANDALAKEKC,TIXIER, 
%          OWOLABI, REVIL, WYLGREG 

% Written by R. E. Rasolovoahangy, June 2000

if nargin==0   
   prompt={'Porosity (Phi)','Irreducible Water Saturation (Swr)'};
   defans={'[0.01:0.01:0.35]','0.15'};
   getpar=inputdlg(prompt,'Timur Model Default Input Parameters',1,defans);
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
   e = Tim(param{1},param{2});  
else
e = Tim(Phi,Swr);
end
   
function e = Tim(Phi,Swr)    
   
K = 0.136.*(((100.*Phi).^4.4)./((100.*Swr).^2));

semilogy(Phi,K,'k-');
set(gca,'fontsize',9);
xlabel('Porosity','fontsize',9);
ylabel('Permeability (md)','fontsize',9);
title('Timur Permeability versus Porosity','fontsize',10);
grid on;
hold on;
e = [Phi K];
