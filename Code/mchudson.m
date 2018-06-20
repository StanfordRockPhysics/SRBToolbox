function [C,den,defans,armc,cdmc]=mchudson(vp,vs,ro,defans,cd,ar,dcd,dar,kfl,rofl,ax,N)
% function [C,den,defans,armc,cdmc]=mchudson(vp,vs,ro,defans,cd,ar,dcd,dar,kfl,rofl,ax,N)
%Monte-Carlo simulation of aligned cracked rocks over distribution of background
%rock matrix properties, and distributions of crack density and aspect ratio.
%Generates Monte Carlo realizations for vp,vs,ro,cd,ar and then does
%Hudson fracture modeling.
%Required input arguments:
% vp, vs   - Rock matrix P and S wave velocities (vectors)
%            to estimate cdf for MC draw.
% ro       - unfractured rock matrix density (vector)
%optional inputs (if not given, brings up a dialog box):
% cd       - crack density of the single set of aligned cracks
% ar       - aspect ratio
% dcd,dar  - type of distribution for cd and ar (0 log-uniform, 1 normal,
%          2 nonparametric, 3 scalar);
%          NOTE:- for uniform distribution cd (or ar) is a 2 element vector
%          (min and max value for cd or ar)
%               - for normal distribution cd (or ar) is a 2 element
%                vector (average and standard deviation)
% kfl     - fluid bulk modulus
% rofl    - fluid density
% ax      - axis of crack normals (1 or 3 usually)
% N       - number of Monte Carlo realizations
% defans  - default answers to use in the dialog box. Pass in [] to
%           use built-in defaults.
% output arguments
% C          - Stiffness matrix, 6x6xN array
% den        - density of fractured fluid filled rock
% defans     - current crack and fluid parameters used in the MC simulations
% armc       - Monte Carlo simulated draws of aspect ratios
% cdmc       - Monte Carlo simulated draws of crack densities
%With no output arguments, plots simulated anisotropy parameters.
%
% See also: HUDSON, MONTE, CTI2V

%Written by Diana Sava, July 2000.


if nargin <12
  prompt ={'cd distrib: 0-unif [min max] 1-norm [avg std] 2-nonparam 3-scalar','ar distrib: 0-unif [min max] 1-norm [avg std] 2-nonparam 3-scalar','crack density','aspect ratio','Kfl fluid bulk modulus','rofl-fluid density','ax -axis of symmetry','N -number of Monte Carlo realizations'};
 if nargin ==3
  defans={'0','0','[0.001 0.1]','[0.001 0.1]','1.3e6','6.5','3','100'};
end
  if isempty(defans);
  defans={'0','0','[0.001 0.1]','[0.001 0.1]','1.3e6','6.5','3','100'};
end


  answer   = inputdlg(prompt,'MChudson',1,defans);
  dcd      = str2num(answer{1});
  dar      = str2num(answer{2});
  cd       = str2num(answer{3});
  ar       = str2num(answer{4});
  kfl      = str2num(answer{5});
  rofl     = str2num(answer{6});
  ax       = str2num(answer{7});
  N        = str2num(answer{8});
end
defans=answer;

if length(cd)==1
  cdmc = cd*ones(N,1);
end

if length(ar)==1
  armc = ar*ones(N,1);  
end

if length(cd)==2
  if dcd==0;
    ymin=-log10(cd(1));
    ymax=-log10(cd(2));
    y=ymax+(ymin-1)*rand(N,1);
    cdmc=10.^-y;
  elseif dcd==1;
    cdmc=cd(1)+cd(2)*randn(N,1)
  end
 end 

 if length(ar)==2
if dar==0;
    ymin=-log10(ar(1));
    ymax=-log10(ar(2));
    y=ymax+(ymin-1)*rand(N,1);
    armc=10.^-y;
  elseif dar==1;
    armc=ar(1)+ar(2)*randn(N,1)
  end

end 

if length(cd)>2
  cdmc=monte(cd(:),N);
end

if length(ar)>2
  armc=monte(ar(:),N);
end

vpvsrosim=monte([vp(:),vs(:),ro(:)],N);
vpmc=vpvsrosim(:,1); vsmc=vpvsrosim(:,2); romc=vpvsrosim(:,3);
[K,G]=v2ku(vpmc,vsmc,romc);
[C,den]=hudson(cdmc,armc,kfl,rofl,K,G,romc,ax); 

if nargout==0
subplot(321)
plot(e,cdmc,'.'); xlabel('crack density'); ylabel('Thomsens epsilon');
subplot(322)
plot(e,armc,'.');xlabel('crack aspect ratio'); ylabel('Thomsens epsilon');
subplot(323)
plot(g,cdmc,'.'); xlabel('crack density'); ylabel('Thomsens gamma');
subplot(324)
plot(g,armc,'.');xlabel('crack aspect ratio'); ylabel('Thomsens gamma');
subplot(325)
plot(d,cdmc,'.'); xlabel('crack density'); ylabel('Thomsens delta');
subplot(326)
plot(d,armc,'.');xlabel('crack aspect ratio'); ylabel('Thomsens delta');
end;

  
  
 
  
  

