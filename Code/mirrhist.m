function mirrhist(data1,data2,n1,n2,col1,col2)
%function mirrhist(data1,data2,n1,n2,col1,col2)
%Horizontal mirror histogram plots of two data sets around a central line.
%DATA1, DATA2, column vector of two data sets to be compared.
%N1, N2, number of bins for histograms of data1 and data2.
%COL1, COL2 (optional inputs), colors (rgb triplet or string) for the 
%two histograms.

%Written by T. Mukerji

if nargin<5, col1=[0 0.8 0.9]; col2=[0.8 0.3 0]; end;
if nargin<3, n1=10; n2=10; end;
[f1,bin1]=hist(data1,n1); [f2,bin2]=hist(data2,n2); 
f1=f1./sum(f1); f2=f2./sum(f2);

h1=barh(bin1,-f1,1); hold on; h2=barh(bin2,f2,1);
set(h1,'facecolor',col1); set(h2,'facecolor',col2);
set(gca,'xticklabel',num2str(abs(str2num(get(gca,'xticklabel')))));
ylmt=get(gca,'ylim');
plot([0 0],ylmt,'-b'); ylim(ylmt); hold off;
%set(gca,'color',[0.3 0.5 0.]);


