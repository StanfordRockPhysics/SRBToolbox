function h1=scatterqq(x,y,s, clim)
%        h1=scatterqq(x,y,s, clim)
%
% Makes a "quick" scatter plot.  Similar to Matlab's scatter, except that
% this version takes the data, divides it into 64 data ranges and makes
% regular plots with different colors.  Therefore, it is very much faster
% than scatter.
% This version will not work with LSELECT, because each color is a line
% with a different length.

% Alternate version SCATTERQ plots each color as a vector of length equal
% to the original data, with lots of NANs for the invisible points.

% G. Mavko, April 2004

smin = min(s);
smax = max(s);
smin = smin - .01*(abs(smax-smin));
smax = smax + .01*(abs(smax-smin));
% answer=inputdlg({'Min value of color range','Maximum value:'},'Quick Scatter Plot',1,{num2str(smin) num2str(smax)});
% smin = str2num(answer{1});
% smax = str2num(answer{2});

if nargin == 4, smin = clim(1); smax = clim(2); end

linecolor = colormap;
n = size(linecolor,1);
srange = linspace(smin, smax, n+1);
caxis([smin smax]);

% this piece tricks the colorbar by making a dummy imagesc plot, then
% deleting it.
h=imagesc(smin+(smax-smin)*rand(200,1));hold on;h1=colorbar;p=get(gca,'position');delete(h);

% now make the separate plots
for k = 1: n
     inds = srange(k) <= s & srange(k+1) > s;
     h=plot(x(inds),y(inds),'o'); 
set(h,'markeredgecolor',linecolor(k,:),'markerfacecolor',linecolor(k,:),'markersize',3); 

     hold on;
end;

% fix the axes
%%%set(gca,'position',p);
axis xy
axis auto
