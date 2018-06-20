function [mat]=smooth(matf,rangei,rangej,filtre);
%function [mat]=smooth(matf,rangei,rangej,filtre);
% function used by bornfilt. Usually need not be called seperately by user.
% see BORNFILT
% Used to truncate and smooth a 2-D filter.
% matf = input filter, mat=output truncated and smoothed filter
% size(matf)=[2*ni,2*nj];
%
%function [mat]=fftshift(ifft2(smooth(fft2(fftshift(matf)),rangei,rangej)));
% mat=matf(2*range,2*range) + 1/2*(1+cos(alpha*x))*(1+cos(beta*z))

figure(1);
ifi=-1;ifj=-1;
while ((ifi~=1)|(ifj~=1))
clf;
colormap(gray);
imagesc(real(fftshift(fft2(fftshift(matf)))));
title('spatial filter');
xlabel('horizontal axis (j-coordinate for a matrix (i,j))')
ylabel('vertical axis (i-coordinate for a matrix (i,j))')
hold on;
[ni,nj]=size(matf);
rangeimax=ni/2;
rangejmax=nj/2;
% special case of smoothing
rangei=ni/4;
rangej=nj/4;
% end of special case of smoothing
rangeimin=min(10,rangeimax);
rangejmin=min(10,rangejmax);
plot([0.5 (nj+0.5)],[(ni/2-rangei+1) (ni/2-rangei+1)])
plot([0.5 (0.5+nj)],[(ni/2+rangei) (ni/2+rangei)])
plot([(nj/2-rangej+1) (nj/2-rangej+1)],[0.5 (ni+0.5)])
plot([(nj/2+rangej) (nj/2+rangej)],[0.5 (0.5+ni)])
hold off
rangei
rangei0=rangei;
rangei=-1;
while((rangei<rangeimin)|(rangei>rangeimax))
rangei=input('rangei ?  : ');
if isempty(rangei)
   rangei=rangei0;
   ifi=1;
else 
   ifi=0;
end
rangei=round(rangei);
ifi=1;
end
rangej
rangej0=rangej;
rangej=-1;
while((rangej<rangejmin)|(rangej>rangejmax))
rangej=input('rangej ?  : ');
if isempty(rangej)
   rangej=rangej0;
   ifj=1;
else
   ifj=0;
end
rangej=round(rangej);
ifj=1;
end
end


if filtre ~= 5
matf=fftshift(fft2(fftshift(matf)));
alphai = pi/rangei;
alphaj = pi/rangej;
aux_mat=zeros(rangei*2,rangej*2);
for i=1:2*rangei;
   di=abs(i-rangei-0.5);
   for j=1:2*rangej;
      dj=abs(j-rangej-0.5);
      aux_mat(i,j)=0.25*(1+cos(alphai*di))*(1+cos(alphaj*dj));
   end
end
ni=ni/2;nj=nj/2;
mat=zeros(2*ni,2*nj);
aux_mat=matf((ni-rangei+1):(ni+rangei),(nj-rangej+1):(nj+rangej)).*aux_mat;
mat((ni-rangei+1):(ni+rangei),(nj-rangej+1):(nj+rangej))=aux_mat;
mat=real(fftshift(ifft2(fftshift(mat))));
else
matf=fftshift(fft(fftshift(matf(:,nj/2))));
alphai = pi/rangei;
aux_mat=zeros(rangei*2,1);
for i=1:2*rangei;
   di=abs(i-rangei-0.5);
   aux_mat(i,1)=0.5*(1+cos(alphai*di));
end
ni=ni/2;
mat=zeros(2*ni,nj);
aux_mat=matf((ni-rangei+1):(ni+rangei),1).*aux_mat;
mat((ni-rangei+1):(ni+rangei),nj/2)=aux_mat;
mat(:,nj/2)=real(fftshift(ifft(fftshift(mat(:,nj/2)))));
mat(:,(nj/2+1))=mat(:,nj/2);
end
