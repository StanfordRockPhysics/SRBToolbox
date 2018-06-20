function [x,r]=minrescg(F,x0,d,itr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  minrescg.m -minimum residual CG linear template.
%  This function is a template that minimizes the 
% residual r=Fx-d, using the conjugate gradient methods.
% {in other word, if you have an equation Fx=d, then you will get x=inv(F)d.}
%  
%  function [x,r]=minrescg(F,x0,d,itr)
%  x0: initial guess -size [m,1]
%  F: matrix.(operator in matrix form) -size [n,m]
%  d: data -size [n,1] 
%  itr: # of iterations
%  minrescg calls cgplus.m function for 1 conjugate gradient step.

%  last modified Jan, 3 1997. Ran Bachrach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sx=0;sr=0;
r=F*x0-d;
for i=0:itr
dx=F'*r;
dr=F*dx;
[x,dx,sx,r,dr,sr]=cgplus(itr,x0,dx,sx,r,dr,sr);
x0=x;
end

