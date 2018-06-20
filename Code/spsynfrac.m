function y=spsynfrac(beta,n)
% Y = SPSYNFRAC(BETA,N) generates a 1-D fractal series Y of length N and 
% a power spectrum that falls off as frequency^BETA. 
% mean(Y) = 0; var(Y)= 1.
% Uses spectral synthesis method.
%
% See also SPSYNEXP, SPSYNGS, SPSYNFRAC2, SPSYNGS2

%Written by T. Mukerji

	randn('seed', sum(100*clock));
%	n=2*n;
	f=(1/n)*[1:1:(n/2-1)]; f=f';
	x=randn(n,1); ph=rand(n,1);
	rad= (f .^(beta/2)).*x(1:(n/2-1));
	phase=2*3.1415927.*ph;
	y=rad .*cos(phase(1:n/2-1))+sqrt(-1)*rad .*sin(phase(1:n/2-1));
	y=[0.0;y;0.0;conj(flipud(y))];
	y=ifft(y); 
%	y=y(n/4+1:3*n/4);
	 y=real(y);
	 y=y-mean(y); y=(1/std(y))*y;
