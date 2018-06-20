function y=spsyngs3(beta,n)
% Y = SPSYNGS3(BETA,N) generates a 3-D random field Y of size 2n+1 by 2n+1
% by 2n+1 and
% a Gaussian autocorrelation function with correlation length BETA.
% mean(Y) = 0; var(Y)= 1.
% Uses spectral synthesis method.
%
% See also SPSYNEXP, SPSYNFRAC, SPSYNGS, SPSYNFRAC2

% Written by Ran Bachrach.  Based on spsyngs.
        randn('seed', sum(100*clock));
%       n=2*n;

        fx=(2*pi)*(1/(2*n+1))*[-n:1:n];nx=length(fx);
        fy=(2*pi)*(1/(2*n+1))*[-n:1:n];ny=length(fy);
    fz=(2*pi)*(1/(2*n+1))*[-n:1:n];nz=length(fz);
        f=zeros(nx,ny,nz);
[fxx,fyy,fzz]=ndgrid(fx,fy,fz);
f=sqrt(fxx.^2+fyy.^2+fzz.^2);
x=(1+0.2*randn([nx,ny,nz]));
        ph=rand(nx,ny,nz);
        powsp=beta*sqrt(pi)*exp(-(beta^2/4)*f.^2);
        rad= sqrt(powsp).*x;
%       rad= (f .^(beta/2)).*x;
        %rad(n+1,n+1)=0;
        phase=2*3.1415927.*ph;
        y=rad .*cos(phase)...
           +sqrt(-1)*rad .*sin(phase);
        %y=[y;conj(fliplr(flipud(y(1:n,:))))];
    m=zeros(size(y));
    n2=(nx-1)/2;
    m(n2+1:end,n2+1:end,n2+1:end)=y(n2+1:end,n2+1:end,n2+1:end);
    %FFT symetries
    %m(n2+1,n2+1,n2+1)=0;%zero DC



    m(1:n2,:,:)=flipdim(conj(m(n2+2:end,:,:)),1);
    m(:,1:n2,:)=flipdim(conj(m(:,n2+2:end,:)),2);
    m(:,:,1:n2)=flipdim(conj(m(:,:,n2+2:end)),3);

    m=fftshift(m);


    y=ifftn(m);
        %clear m;
        y=real(y);
        y=y-mean((y(:))); y=(1/std2(y))*y;

