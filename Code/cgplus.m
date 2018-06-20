function [x,g,s,r,gg,ss]=cgplus(itr,x0,g,s,r,gg,ss);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CGPLUS: Conjugate Gradient Solver
%  function [x,g,s,r,gg,ss]=cgplus(itr,x0,g,s,r,gg,ss);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates 1 step in Conjugate Gradient descent.
%  Modified from ratfor program by Jon Claerbout, from his book
%  3 dimensional filtering.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  x=solution,r=residual
%  g=gradient,gg=conjugate gradient
%  s=step,ss=conj. step
%
% See also MINRESCG

%  last modified dec 31,1996, Ran Bachrach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
if itr<1
if gg'*gg==0 error('first grad=0!!!')
end
alfa=-r'*gg/(gg'*gg);
beta=0;
else
  gdg=gg'*gg;
  sds=ss'*ss;
  gds=gg'*ss;
   if gdg==0 error('Grad =0!!!')
    elseif sds==0 
     alfa=-gg'*r/(gg'*gg);,beta=0;
    
   else
determ=gdg*sds*max(1.-gds/gdg*gds/sds,1e-12);
gdr=-gg'*r;sdr=-ss'*r;
alfa=(sds*gdr-gds*sdr)/determ;
beta=-(gds*gdr+gdg*sdr)/determ;
end
end
s=alfa*g+beta*s;
ss=alfa*gg+beta*ss;
x=x0+s;
r=r+ss;

