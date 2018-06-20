function lplotseis(log)
%function lplotseis(log)
%LPLOTSEIS Plots logs and filled wiggle trace seismogram with depth 
%(or time) axis %increasing downwards. Input is a matrix
%with depth (or time) in the first column, various logs in the other columns,
%and the seismogram in the last column.
%
%See also LOGAX and LZOOM to change the depth axis of the log plots.
%See also LPLOTSEIS2 

%Written by T. Mukerji

n=size(log,2); m=n-1;

clf
subplot('position',[0 0 1 1]); 
for k=1:m-1
subplot(1,m,k)
plot(log(:,k+1),log(:,1)),set(gca,'ydir','reverse'),end; 
ss=log(:,end)*ones(1,5);
subplot(1,m,m); seisplot(ss,log(1,1),log(2,1)-log(1,1));
xlabel(''); ylabel('');
axis([-1 5 log(1,1),log(end,1)]); logax([log(1,1),log(end,1)]);
lzoom yon

