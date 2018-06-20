function lplotseis2(log,sc)
%function lplotseis2(log)
%LPLOTSEIS2 Plots logs and 2 filled wiggle trace seismogram with 
%depth (or time) axis increasing downwards. Input is a matrix
%with depth (or time) in the first column, various logs in the other columns,
%and the 2 seismograms in the last 2 columns.
%See also LOGAX to change the depth axis of the log plots.
%See also LPLOTSEIS

%Written by T. Mukerji

n=size(log,2); m=n-1;

if nargin==1, sc=log(:,end); end;

clf
subplot('position',[0 0 1 1]); 
for k=1:m-2
subplot(1,m,k)
plot(log(:,k+1),log(:,1)),set(gca,'ydir','reverse'),end; 
ss=log(:,end-1)*ones(1,5); sc=sc*ones(1,5);
subplot(1,m,m-1); seisrwb(ss,log(1,1),log(2,1)-log(1,1),0,1,0,sc);
xlabel(''); ylabel('');
axis([-1 5 log(1,1),log(end,1)]); logax([log(1,1),log(end,1)]);
ss=log(:,end)*ones(1,5);
subplot(1,m,m); seisrwb(ss,log(1,1),log(2,1)-log(1,1),0,1,0,sc);
xlabel(''); ylabel('');
axis([-1 5 log(1,1),log(end,1)]); logax([log(1,1),log(end,1)]);
lzoom yon

