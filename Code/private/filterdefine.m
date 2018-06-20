function [filt,vec,err,S]=filterdefine(dataref,xx,obs,nfilt);

if ~iscell(xx) xx{1}=xx; end;
ndim=length(xx);
if length(obs)==1, obs=ones(1,ndim)*obs; end
obs=obs(:)';

correlatederror=1;

if any(obs~=0)
   
   
   S=std(dataref);
   
%S=cov(dataref(dataref(:,end)==max(dataref(:,end)),1:ndim));
%S=sqrt(diag(S))';
for i=1:ndim
dx(i)=mean(diff(xx{i}));
end

SS=S./dx;
ss=obs.*SS;
%nfilt=2*floor(ss*6/.1*11/100/2)+1;
scaled_data=dataref(:,1:ndim).*repmat(obs./dx,[length(dataref) 1]);
cv=cov(scaled_data);
if correlatederror==0,
cv=diag(diag(cv));
end
filt=fspecials('ga',nfilt,cv);


%error=4*ss.*dx;

error=obs.*S;

e1=4*ss(1)*dx(1);
for i=1:ndim
avx(i)= max(xx{i})-min(xx{i});
end;
err=repmat(error,[2 1]);
err(2,:)=err(2,:)./avx;
for i=1:ndim
   vec{i}=xx{i};
%vec{i}=extrapolate(xx{i},nfilt(i));
end

else
for i=1:ndim;
vec{i}=xx{i};
end
filt=1;
err=zeros(2,ndim);

end
