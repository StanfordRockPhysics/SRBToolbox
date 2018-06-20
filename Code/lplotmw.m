function lplotmw(varargin)
%function lplotmw(wlog1,defclm1,wlog2,defclm2,....,clm)
%Plots selected logs from multiple well log data matrices WLOG1, WLOG2,...
%First column in each matrix is the depth, and may be different in
%each WLOG matrix. 
%DEFCLM1, DEFCLM2,...are structure arrays defining the columns of the 
%corresponding WLOG matrices. See LPLOT for the definition of DEFCLM.
%Any number of WLOG and corresponding DEFCLM inputs may be given.
%CLM is a cell array of strings selecting which column to plot. The string
%should correspond to the field names in DEFCLM.
%
%See also LPLOT, LOGAX

%Written by T. Mukerji

clm=varargin{end}; nlog=(length(varargin)-1)/2; ymin=[]; ymax=[];
for k=1:length(clm)
figure(k)
 for kk=1:2:length(varargin)-1
   wlog=varargin{kk}; wlogclm=varargin{kk+1}; fn=fieldnames(wlogclm);
   wlog(wlog<=-999)=nan;
   [tmp,ifn]=intersect(fn,clm{k}); 
   subplot(1,nlog,fix(kk/2)+1)
   if ~isempty(ifn)
   colno=getfield(wlogclm,fn{ifn});
   plot(wlog(:,colno),wlog(:,1)),set(gca,'ydir','reverse'); 
   xlabel(upper(clm{k}));
   ylmt=get(gca,'ylim'); ymin=min([ymin,ylmt]); ymax=max([ymax,ylmt]);
   end;
 end;
set(findobj('type','axes'),'ylim',[ymin,ymax]);
lzoom yon;
end;  

