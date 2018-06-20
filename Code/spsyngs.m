function y=spsyngs(beta,n)
% Y = SPSYNGS(BETA,N) generates a 1-D series Y of length N and 
% a gaussian autocorrelation function with correlation length BETA
% mean(Y) = 0; var(Y)= 1.
% Uses spectral synthesis method.
%
% See also SPSYNEXP, SPSYNFRAC, SPSYNFRAC2, SPSYNGS2

%Written by T. Mukerji

	randn('seed', sum(100*clock));
%	n=2*n;
	f=(1/n)*[1:1:(n/2-1)]; f=2*pi*f';
	ph=rand(n,1);
        x=1+0.2*randn(size(f));
        powsp=beta*sqrt(pi)*exp(-(beta^2/4)*f.^2);
        rad= sqrt(powsp).*x;
	phase=2*3.1415927.*ph;
	y=rad .*cos(phase(1:n/2-1))+sqrt(-1)*rad .*sin(phase(1:n/2-1));
	y=[0.0;y;0.0;conj(flipud(y))];
	y=ifft(y); 
%	y=y(n/4+1:3*n/4);
	 y=real(y);
	 y=y-mean(y); y=(1/std(y))*y;
