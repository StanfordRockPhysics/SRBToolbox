function lplot(wlog,defclm,clm)
%function lplot(wlog,defclm,clm)
%LPLOT Plots well logs with depth axis increasing downwards. Input is 
%either a matrix WLOG with depth in the first column, and various log
%curves in the other columns, or is a data structure with fields 
%'header' and different column names as obtained from LOADLAS.
%
%Optional input: DEFCLM is a cell array of strings defining the column names.
%This should have as many elements as the number of columns in data matrix. 
%When WLOG is a structure, an empty matrix [] may be passed in for
%DEFCLM as the column names are taken from the fieldnames of WLOG.
%Optional input: CLM is a cell array of strings selecting which columns to plot.
%The strings should correspond to the names in DEFCLM (or fieldnames of WLOG).
%When WLOG is a structure (as obtained from LOADLAS), calling LPLOT with
%single input argument, LPLOT(WLOG), brings up list dialog box to
%select curves to be plotted.
%Linked y-axis zoom is enabled by default in LPLOT.
%
%EXAMPLE: >> wlog = 1000 x 5 matrix (1st column = depth) 
%         >> lplot(wlog)
%        This will plot all 4 columns of the matrix
%
%         >> defclm={'depth','gr','por','ild','dt'};
%
%This defines gamma-ray to be column 2, porosity, column 3, resistivity,
%column 4, and sonic traveltime, column 5.
%Selected fields from defclm may be plotted as:
%
%         >> lplot(wlog,defclm,{'gr','ild'})
%This will select and plot only the gamma-ray and resistivity columns.
%It is equivalent to
%         >> lplot(wlog(:,[1 2 4]));
%
%EXAMPLE: when wlog is a struct:
%         >> wlog
%
%    header: [64x82 char]
%     depth: [16310x1 double]
%        gr: [16310x1 double]
%      rhob: [16310x1 double]
%      nphi: [16310x1 double]
%       ild: [16310x1 double]
%
%	  >> lplot(wlog) %will bring up dialog box to select curves to plot.
%Command line selection will also work:
%	  >> lplot(wlog,[],{'gr','ild','nphi'});
%
%See also LOGAX and LZOOM to change the depth axis of the log plots.
%See also LOADLAS to load LAS formatted well log data files.

%Written by T. Mukerji, I. Takahashi, 1999.
%Modified from older version.

if isstruct(wlog), 
 if isfield(wlog,'header'), wlog=rmfield(wlog,'header'); end;
 colnames=fieldnames(wlog);
   if nargin==1, 
   [sel,ok]=listdlg('liststring',colnames(2:end),... %don't list depth
                    'name','LPLOT-SRBTOOLS',...
                    'promptstring','Select curves to plot'); 
   clm=colnames(sel+1);
   end;
% log(:,1)=wlog.depth;
 log(:,1)=getfield(wlog,colnames{1});
 for k=1:length(clm),log(:,k+1)=getfield(wlog,clm{k}); end;
 log(log<=-999)=nan;
 xlbl=clm;
else                     %% when wlog is not a struct
 wlog(wlog<=-999)=nan;
 if nargin==1, 
  log=wlog; xlbl=[]; 
 elseif nargin==2
   xlbl=defclm(2:end); log=wlog;
 elseif nargin==3
  if ~isempty(defclm)
   colnos=find(ismember(defclm,clm)); xlbl=clm;
   colnos=[1,colnos(:)']; log=wlog(:,colnos);
  else
   xlbl=clm; log=wlog;
  end
 end;
end

n=size(log,2); m=n-1;

for k=1:m
subplot(1,m,k)
plot(log(:,k+1),log(:,1),'b'),set(gca,'ydir','reverse'); grid on; hold on;
if ~isempty(xlbl), 
xlabel(upper(xlbl{k}),'fontsize',10), 
if ismember(xlbl{k},{'ild','ilm'}), set(gca,'xscale','log'), end;
end;
if k>1, set(gca,'yticklabel',[]); end;
if k==1, ylabel('DEPTH','fontsize',10), end;
set(gca,'fontsize',10)
end;
lzoom yon;
