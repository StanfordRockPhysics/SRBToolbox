function [ent]=entropdf(p,units);
%ENTROPDF  Shannon's Information Entropy
%  ENTROPDF(p) calculates Shannon's information entropy from 
%  probability distribution function. 
%
%  p:	   Input probability distribution function.
%      	   (Can have any number of dimensions)
%          Scaling will be applied if sum(p) is not equal to 1.
%  units:  Defines the unit of entropy. Can be either a number or
%	   one of the characters, 'nat','e', 'bit', or 'dit.'
%
%  Entropy is defined as
%  	                      N
%  	ENTROPY(p,units) = - sum p(n) log p(n) / log (units);
%  	                     n=1 
%  
%  ENTROPDF(p), ENTROPDF(p,'nat') or ENTROPDF(p,'e') compute the entropy with 
%  base e, which is equal to ENTROPY(p,exp(1)). Its unit is called 'nat.'
%  
%  ENTROPDF(p,'bit') computes the entropy with base 2, which is equal to
%  ENTROPDF(p,2). Its unit is called 'bit.'
%  
%  ENTROPDF(p,'dit') computes the entropy with base 10, which is equal to
%  ENTROPDF(p,10). Its unit is called 'dit.'
%

%  Written by Isao Takahashi, 1/20/1999.
%  Name changed to ENTROPDF from older version ENTROPY so as not to clash
%  with MATLAB 7 Image Processing toolbox function ENTROPY. Tapan Mukerji,
%  2/3/2005

if nargin==1, units='e'; end;
if any(lower(units(1))=='en'),	base=exp(1);
elseif	lower(units(1))=='b',	base=2;
elseif	lower(units(1))=='d',	base=10;
elseif	~ischar(units(1)),	base=units;
else	error('The second argument should be either of nat, bit, dit or a number');
end
p=p(:);
p=p/sum(p);
flag=p~=0;
[ent]=-sum(log(p(flag)).*p(flag));
ent=ent/log(base);
