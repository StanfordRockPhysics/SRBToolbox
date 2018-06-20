function agx=agc(xx,t,dt,opt);
% function agx=AGC(xx,t,dt,opt);
%
% This function calculates the smooth AGC display of the data.
% xx - data matrix size(nt,nx) (Can be shot gather or cdp or just one trace)
% t  - time length of the agc operator
% dt - sampling rate, or t=#of samples to be averaged dt=1.
% opt- normalization option.  opt=2 gives RMS; opt=1 normalize by sum of 
% abs values. Defualt is RMS
% 
% see Yilmaz for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Written by R. Bachrach
%Modified by T. Mukerji

if nargin==3 opt=2; end
ng=fix(t/dt);

[nt,ny]=size(xx);
 
for i=1:fix((nt/ng))

   v=[ng*(i-1)+1:ng*i];

  if opt==2
   M(i,:)=sqrt((sum(xx(v,:).^2))/ng);
   elseif opt==1
   M(i,:)=(sum(abs(xx(v,:))))/ng;
  end
end

x1=[1; [ng/2:ng:ng*fix(nt/ng)]'; nt]; M1=[M(1,:); M; M(end,:)]; xI=[1:nt]';
MI=interp1(x1,M1,xI);   
agx=xx./MI;



