function [ar, clm, acfslice, lag, imgacf]=acfprofile(img)
%function [ar, clm, acfslice, lag, imgacf]=acfprofile(img)
%computes radial sections of autocovariance function from input image.
%correlation length is taken to be the lag at which the correlation
%falls off to 1/e.
%outputs: ar = correlation anisotropy ratio a_max/a_min
%         clm = median correlation length (1/e fold distance)
%         acfslice = matrix with colms representing slices of acf
%                    along profiles from 0 to 180 degrees azimth.
%         lag = lag corresponding to rows of acfslice
%         imgacf = 2D acf of the input image

%written by T. Mukerji, March 2003


img=double(img);

imgft=fft2(img-mean2(img));
imgsp=imgft.*conj(imgft);
imgacf=fftshift(ifft2(imgsp));
imgacf=real(imgacf);
imgacf=imgacf./max(imgacf(:));


[nx,ny]=size(imgacf);
[ii,jj]=find(imgacf==1); %% 0-lag indices at center of image acf
lagy=[1-ii:nx-ii]; lagx=[1-jj:ny-jj];
r=min(ii,jj);

npt =2*nx; theta=[0:180]; thetar=theta*pi/180;
lag=zeros(npt,length(theta)); acfslice=zeros(npt,length(theta));
for kk=1:length(theta)
[cx,cy,c]=improfile(lagx,lagy,imgacf,[0, r*cos(thetar(kk))],[0, -r*sin(thetar(kk))],npt,'bilinear');
lag(:,kk)=sqrt(cx.^2 + cy.^2);
acfslice(:,kk)=c;
end;

cl=max(lag.*(acfslice>exp(-1)));
clm = median(cl);
ar = max(cl)/min(cl);

if nargout==0
figure; imagesc(imgacf); axis image;
figure; plot(cl);
plot(lag, acfslice,'-k'), xlabel('lag'), ylabel('ACF')
clm, ar
end
