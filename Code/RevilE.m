function e = RevilE(Phi,d) 
%function e = RevilE(Phi,d) 
%
%Revil and al model permit to predict permeability
%in very shaly rock
%
%K= (d.^2).*(Phi.^(3.*m))./24
%
%where: d = average grain diameter (micron-m),
%       Phi = porosity, scalar or vector, and
%       m = cementation exponent, 1.5.
%
%The output is two column matrix, porosity in fraction
%in first column, and permeability (in md) in second
%column.
%
%VALIDITY:Agree with very shaly rocks, Waxaman and 
%         Smith 1968.
%REFERENCE:"Electrical Conductivity, Permeability, 
%          streaming potential and electro-osmosis in
%          granual porous media: a unified approach'"
%          Revil, Glover, Pezard, Zamora, 1997
%
%See also BERNABE, BLOCH,COATDUM, COATES,COATDUM, 
%         KOZCARM, MODKOZCARM, PANDALAKE,PANDALAKEKC,
%         TIMUR, TIXIER, OWOLABI, WYLGREG 

% Written by R. E. Rasolovoahangy, June 2000

if nargin==0   
   prompt={'Porosity (Phi)','Average Grain Diameter (d), micronmeter'};
   defans={'[0.01:0.01:0.35]','250'};
   getpar=inputdlg(prompt,' Revil Model Default Input Parameters',1,defans);
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
   e = Rev(param{1},param{2});  
else
e = Rev(Phi,d);
end
   
function e = Rev(Phi,d) 

K= 1000.*(d.^2).*(Phi.^4.5)./24;

semilogy(Phi,K,'c-');
set(gca,'fontsize',9);
xlabel('Porosity','fontsize',9);
ylabel('Permeability (md)','fontsize',9);
title('Revil Permeability versus Porosity','fontsize',10);
grid on;
hold on;
e = [Phi K];

