function seisfillcolor(colorfill,colorline,colormp,colorfillneg,polarity)
% seisfillcolor(colorfill,colorline,colormp,colorfillneg,polarity)
%
% Changes colors of seismic plots made with seisplot.m and seisrwb.m
% This program takes an existing plot (the active figure in a Matlab
% session) and allows several options for coloring the line, fill, and
% background.  In fact, with seisplot, there is no need to have separate 
% seisplot.m and seisrwb.m programs, because seisfillcolor.m allows
% conversion from one to the other.
%
% --> It is important to use the newer version of seisplot or seisrwb, which
% give the axes the 'tag' 'seisplot'
%
% Inputs:
%    colorfill  color for positive wiggle fill.
%               = 'none' or '' for just wiggle
%               = 'r','b','g','m','c','k','y','w'
%                  or [r g b], ie., [.8 .8 .8] for light gray
%               defaults to 'k'
%    colorline  (optional) color for wiggle lines.  
%               = 'none' or '' for no line
%               = 'r','b','g','m','c','k','y','w',
%                  or [r g b], ie., [.8 .8 .8] for light gray
%               Defaults to colredge=colrfill
%    colormp    colormap to for the underlying image
%               = 'none' or '' to hide background image
%    colorfillneg similar to input colorfill, except this allows trough to
%               be filled with a separate color.  
%               = '' or 'none' will give no trough fill
%               defaults to no trough fill
%    polarity   = 1 or -1
%
% colorfill and  colorline can be specified as a string such as 'k' or as a vector [1 0 0]
% use single quotes in 'k', 'g', etc; but no quotes if [1 0 0].
% set to [] or 'none' to delete line or fill
% colormap can be set to one of the preset maps: seis, rwb, a 64x3 array,
%               or to fixed solid colors 'k', 'g', etc; or rgb colors   [1 0 0]
%              
%    use single quotes in 'k', 'g', etc;  no quotes if [1 0 0]; no quotes
%    if rwb, seis, jet, etc
%
% Example 1:  
%   >figure; seisplot(rand(200,1));
%   >seisfillcolor('k','k',seis);             
%   will yield a plot with black wiggles, black fill, and a colored
%   background having colormap seis.
% Example 2:
%   >figure; seisplot(rand(200,1));
%   >seisfillcolor('','k',rwb);             
%   will yield a plot with black wiggles, no fill, and a colored
%   background having colormap rwb.
% Example 3:
%    figure; seisplot(rand(200,1));
%    seisfillcolor('b','k','','r');             
%   will yield a plot with black wiggles, peaks filled in blue, 
%   troughs filled in red, and no colored background

% Written by Gary Mavko


% take care of line and patch defaults
if nargin==0, colorfill='k'; colorline='k'; end;
if nargin==1, colorline=colorfill;  end;
if nargin<4, colorfillneg = []; end;
if nargin<5, polarity=1; end;
coloredge=colorfill;
coloredgeneg=colorfillneg;
 
% Find the subplot in the current figure that is a seismic section
haxis=findobj(gcf,'tag','seisplot');

for i=1:length(haxis),
    axes(haxis(i));
    % Take care of background image first.
    himage= findobj(haxis(i),'type','image');           % handle to background image
    if isempty(himage), himage= findobj(haxis(i),'type','surface'); end;
    if nargin>=3,                               % change in colormap had been requested in the input
        if isempty(colormp) | strcmp(colormp,'none'),                    % specified colormap is empty, so delete image
%            cmap = white;
            set(himage,'visible','off');
        else                                    % colormap is specified, if necessary, convert to 64x3 array
            if prod(size(colormp))==1,          % see if the specified colormap is a single color string
                if ismember(colormp,{'r', 'b', 'g', 'm', 'c', 'y', 'k','w'}),
                    i = strmatch(colormp, {'r', 'b', 'g', 'm', 'c', 'y', 'k','w'}, 'exact');
                    carray = [1 0 0; 0 0 1; 0 1 0; 1 0 1; 0 1 1; 1 1 0; 0 0 0; 1 1 1];
                    cmap = repmat(carray(i,:), 64,1);   % convert letter to a colormap array
                else,
                    warndlg('Cannot understand the specified colormap');
                    return;
                end;
            elseif length(colormp)==3,       % specified as constant rgb array
                cmap = repmat(colormp, 64, 1);
            elseif length(colormp)>=4        
                if strcmp(colormp, 'none'), 
                    cmap = white;
                elseif length(colormp)>=64,
                    cmap = colormp;
%                    if polarity<0,cmap=flipud(cmap); end;
                else
                    warndlg('Cannot understand the specified colormap');
                    return;
                end;
            end;
        
            % colormap is now well defined, apply it to background image
            if ~isempty(himage)                 % existing background image exists, so  apply the colormap
                colormap(cmap);
                set(himage,'visible','on');
            elseif isempty(himage)                          % no image exists,but map is specified, so must redraw
                hpatch = findobj(haxis(i),'type','patch');     % look for filled traces (patches)
                hline=findobj(haxis(i),'type','line');         % look for the wiggles (lines)
                xdata=get(hline, 'xdata');                  % get data from wiggles for redrawing background
                if iscell(xdata),  for j=1:length(xdata), x(:,j)=xdata{length(xdata)-j+1}'; end; else, x=xdata; end;
                for j=1:size(x,2), xmean=nanmean(x(:,j)); x(:,j)=x(:,j)-xmean; end; % remove mean from wiggle data
                xl=xlim;
                x0 = xl(1);
                dx = diff(xl)/(size(x,2)+.5);
                x=x/(1.5*dx);   % remove the scale factor applied in seisplot
                xnorm = x/max(max(x));                                                % normalize to max; same as in seisplot

                ydata=get(hline, 'ydata');                  % repeat same process for ydata
                if iscell(ydata),  y=ydata{1}; else, y=ydata(1); end;
                ttop = y(1);
                dt = mean(diff(y));
            
               % save old line and patch colors before redrawing
                old_linecolor = [];
                if ~isempty(hline),  
                    old_linecolor = get(hline,'color');
                    delete(hline); 
                end;
                old_facecolor = [];
                old_edgecolor = [];
                if ~isempty(hpatch), 
                    old_facecolor = get(hpatch,'facecolor');
                    old_edgecolor = get(hpatch,'edgecolor');                
                    delete(hpatch); 
                end;
                axis manual;
                hold on;
                seisrwb(x,ttop,dt,x0,dx,0,xnorm);                  % redraw wiggles and background
                if ~isempty(old_linecolor), set(findobj(gca,'type','line'),'color',old_linecolor{1}); end;            
                if ~isempty(old_facecolor), set(findobj(gca,'type','patch'),'facecolor',old_facecolor{1}); end;
                if ~isempty(old_edgecolor), set(findobj(gca,'type','patch'),'edgecolor',old_edgecolor{1}); end;
                colormap(colormp);
            end;
        end;
    end;
    

    % Work on line color

    hline = findobj(haxis(i),'type','line');
    if ~isempty(hline), 
        if isempty(colorline) | strcmp(colorline,'none'),
            set(hline,'visible','off');
        else
            set(hline,'color',colorline); 
            set(hline,'visible','on');
        end;
    end;

    % Work on peak and trough-filling patch color

    hpatchpos = findobj(haxis(i),'tag','fillpos');
    hpatchneg = findobj(haxis(i),'tag','fillneg');
    if isempty(hpatchpos) | isempty(hpatchneg) | polarity<0,    % one patch or another does not exist, might need to redraw
%        if ~isempty(colorfill) & ~strcmp(colorfill,'none'),
            hline=findobj(haxis(i),'type','line');
            xdata=get(hline, 'xdata');
            if iscell(xdata),  for j=1:length(xdata), x(:,j)=xdata{length(xdata)-j+1}'; end; else, x=xdata; end;
            for j=1:size(x,2), xmean=nanmean(x(:,j)); x(:,j)=x(:,j)-xmean; end; % remove mean from wiggle data
            x=x*polarity;
%            polarity
            xl=xlim;
            x0 = xl(1);
            dx = diff(xl)/(size(x,2)+.5);
            x=x/(1.5*dx);   % remove the scale factor applied in seisplot
            xnorm = x/max(max(x));                                                % normalize to max; same as in seisplot

            ydata=get(hline, 'ydata');
            if iscell(ydata),  y=ydata{1}; else, y=ydata(:,1); end;
            ttop = y(1);
            dt = mean(diff(y));
            axis manual;
            delete(hline);
            delete(findobj(haxis(i),'tag','fillpos'));
            delete(findobj(haxis(i),'tag','fillneg'));
            hold on; seisplot(x,ttop,dt,x0,dx,0,xnorm);
%            hpatchpos = findobj(haxis(i),'tag','fillpos');
%            hpatchneg = findobj(haxis(i),'tag','fillneg');
%            set(hpatchpos,'facecolor',colorfill);
%            set(hpatchpos,'edgecolor',coloredge);
%            set(hpatchneg,'facecolor',colorfillneg);
%            set(hpatchneg,'edgecolor',coloredgeneg);
            set(findobj(gca,'type','line'),'color',colorline);
%        end;
    end;
    hpatchpos = findobj(haxis(i),'tag','fillpos');
    if ~isempty(hpatchpos),                                    % positive patch exists; delete it or color it
        if isempty(colorfill) | strcmp(colorfill,'none'),
            set(hpatchpos,'visible','off');
        else
            set(hpatchpos,'facecolor',colorfill);
            set(hpatchpos,'edgecolor',coloredge);
            set(hpatchpos,'visible','on');
        end;
    end;
    hpatchneg = findobj(haxis(i),'tag','fillneg');
    if ~isempty(hpatchneg),                                    % positive patch exists; delete it or color it
        if isempty(colorfillneg) | strcmp(colorfillneg,'none'),
            set(hpatchneg,'visible','off');
        else
            set(hpatchneg,'facecolor',colorfillneg);
            set(hpatchneg,'edgecolor',coloredgeneg);
            set(hpatchneg,'visible','on');
        end;
    end;
end;
    
    

    

