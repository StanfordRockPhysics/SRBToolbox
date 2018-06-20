function y=spsynfrac2(beta,n)
% Y = SPSYNFRAC(BETA,N) generates a 2-D fractal random field  Y of size
% 2n+1 by 2n+1 and a power spectrum that falls off as frequency^BETA. 
% mean(Y) = 0; var(Y)= 1.
% Uses spectral synthesis method.
%
% See also SPSYNEXP, SPSYNFRAC, SPSYNGS, SPSYNGS2


% Written by Philippe Rio

        val1=4000;
        val2=750;
	randn('seed', sum(100*clock));
%	n=2*n;
	b=(2*pi)*(1/(2*n+1))*[-n:1:n];
	a=(2*pi)*(1/(2*n+1))*[-n:1:n];
        f=zeros(n+1,(2*n+1));
        for j=1:(2*n+1);
            for i=1:(n+1);
                f(i,j)=sqrt(a(i)*a(i)+b(j)*b(j));
            end
        end
	x=(1+0.2*randn((n+1),(2*n+1)));
        ph=rand((n+1),(2*n+1));
	rad= (f .^(beta/2)).*x;
        rad(n+1,n+1)=0;
	phase=2*3.1415927.*ph;
	y=rad .*cos(phase)...
           +sqrt(-1)*rad .*sin(phase);
	y=[y;conj(fliplr(flipud(y(1:n,:))))];
        for i=1:n;
            y(n+1,-i+2*n+2)=conj(y(n+1,i));
        end
        m=zeros(2*n+1);
        for i=1:n+1;
            for j=1:n+1;
               m(i,j)=y(i+n,j+n);
            end
        end
        for i=n+2:2*n+1;
            for j=n+2:2*n+1;
               m(i,j)=y(i-n-1,j-n-1);
            end
        end
        for i=1:n+1;
            for j=n+2:2*n+1;
               m(i,j)=y(i+n,j-n-1);
            end
        end
        for i=n+2:2*n+1;
            for j=1:n+1;
               m(i,j)=y(i-1-n,j+n);
            end
        end
	y=ifft2(m); 
        clear m;
        y=real(y);
	y=y-mean(mean(y)); y=(1/std2(y))*y;
%        y=val1+val2*y;
