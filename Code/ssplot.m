function seisplot(x,t0,dt,nd,sf)
% SEISPLOT(X,T0,DT,ND,SF) 
% Surface seismic wiggle plot with filled positive area
% of traces forming the columns of matrix X against the time
% axis given by beginning time T0 and time interval DT.
% ND is the number of traces skipped between each plotted trace.
% Time axis positive downwards.
% If only one trace is given (X has one column), it is repeated 25 times.
% SF is an optional scale factor which controls the trace spacing. The
% default if nothing is given is 0.5

% Written by Li Teng

if nargin==4, sf=0.5; end; %default scale factor

%if only one trace is given, then plot ten traces here
[nr, nc]=size(x);
if ( (nr==1) | (nc==1) )
  [junk,x]=meshgrid(1:25,x);
end

[nr, nc]=size(x);
if nd ~= 0.
  nc=floor(nc/nd);
end

%choose the traces needed to be plotted
y=zeros(nr,nc);
i=1:nc;
  if nd ~= 0.
    y(:,i)=x(:,(i-1)*nd+1);
  else
    y(:,i)=x(:,i);
  end

%zero mean
m=mean(y);
y=y-m(ones(nr,1),:);

%normalize y by the biggest element
nrm=max(max(y));
y=y./nrm; 

%z is the positive area
z=y.*(y>0);
%zw0 ie z with zero at the two ends will garantee that the
%filled area is between the positive values and mean value
zw0=zeros(nr+2,nc); zw0(2:nr+1,:)=z;
zw0=z;

i=1:nc; scale=(i-1)*sf;
% y(:,i)=y(:,i)-(i-1.)*0.5;  % for vsp
% z(:,i)=z(:,i)-(i-1.)*0.5;  % for vsp
%  y(:,i)=y(:,i)+(i-1.)*0.5; 
%  zw0(:,i)=zw0(:,i)+(i-1.)*0.5;
  y=y+scale(ones(nr,1),:);  % vectorized version of the above two lines
  zw0=zw0+scale(ones(nr,1),:);

ty=t0:dt:t0+(nr-1)*dt;
tz(1)=t0;
tz(2:nr+1)=ty;
tz(nr+2)=tz(nr+1);
tz=tz';
tz=ty;

%fill(tz,zw0,'w'), hold on
%plot(ty,y,'-w'), hold off
%xlabel('time');
fill(zw0,tz,'w'), hold on
plot(y,ty,'-w'), hold off
set(gca,'ydir','reverse')
ylabel('time');

