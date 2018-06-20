function agx=AGC(xx,t,dt,opt);
% function agx=agc(xx,t,dt,opt);
%
% This function calculates the smooth AGC display of the data
% see Yilmaz for details.
% xx - data matrix size(nt,nx) (Can be shot gather or cdp or just one trace)
% t  - time length of the agc operator
% dt - sampling rate, or t=#of samples to be averaged dt=1.
% opt- normalization option.  opt=2 gives RMS; opt=1 normalize by sum of 
% abs values. Defualt is RMS
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==3
opt=2;
end
ng=fix(t/dt);

[nt,ny]=size(xx);
 
for j=1:ny

for i=1:fix((nt/ng))

   v=[ng*(i-1)+1:ng*i];

if opt==2
 M(i)=sqrt((sum(xx(v,j).^2))/ng);
 elseif opt==1
 M(i)=(sum(abs(xx(v,j))))/ng;
end
end

nM=length(M);
x=ng/2:ng:nt;
x1=[1,x,nt];M1=[M(1),M,M(nM)];
xI=1:nt;
MI=interp1(x1,M1,xI);   
agx(:,j)=xx(:,j)./MI';
end



