function WW=weight(W)
%function WW=weight(W)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function calculate a matrix of zero-padding and weight
%  For the absorbing boundary and the free surface boundary condition.
%  WW is a matrix of unity with decaying boundaries  
%  at the Z direction after absorbing boundaries.
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WW=W;

[N,n]=size(W);

nw=N/16; nw=4;
nw34=nw*3/4;

%for i=1:N
%for j=1:N
%---------------------------------------------------
%  The X direction boundaries
%-------------------------------------------------------------------- 

%if i <=nw
i=1:nw; 
WW(:,i)=exp(-(0.25*(nw-i(ones(N,1),:))).^2).*W(:,i);

%end

%if i >=N-nw
i=N-nw:N;
WW(:,i)=exp(-(0.25*(i(ones(N,1),:)-N+nw)).^2).*W(:,i);

%end
%end,end
%---------------------------------------------------
%  The Z direction boundaries (decay+ zero padding)
%--------------------------------------------------- 

%for i=1:N
%for j=1:N


%if j>=N/2+nw34
%j=N/2+nw34:N; 
%WW(j,:)=exp(-(0.25*(j(ones(N,1),:)'-N/2+nw34)).^2).*W(j,:);
j=N-nw:N;
WW(j,:)=exp(-(0.25*(j(ones(N,1),:)'-N+nw)).^2).*WW(j,:);

%end
%if j >=N-(nw-nw34)
%j= N-(nw-nw34):N;
%WW(j,:)=zeros(length(j),N);
%end
%end,end

