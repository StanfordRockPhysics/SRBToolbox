function  [x,y]=dthief(filename)
%
% [x,y]=dthief(filename)
% Allows manual digitizing of data points from a scanned image.
% Reads image file 'filename' entered in single quotes as a string.
%
% for example,
%             [x,y] = dtheif('myscannedimage.jpg')
%
% Program then plots it, and prompts the user to click various points on
% the X and Y axes, to register the orientation and units of the input
% image.  Program makes correction for rotated graphs on the input image.
% Finally, prompts the user to click on the data points to digitize.
% Left-click adds new digitized point.   Right-click deletes previous
% points;  multiple right-click successively deletes points.
% Outputs the (x,y) locations of each clicked data point, in the dimensions
% of the axes.

%written by Gary Mavko July 2003

% read the file
M=imread(filename);

% display the file
himage =figure;
imagesc(M);
hold on;
axis off;

% Flip sense of image y-axis, because it is in Matlab IJ format
ysign = -1;

% Calibrate X-axis

h = text(0,0,'click any point along the LEFT half of the X-(horizontal) axis');
set(h,'fontsize',16); set(h,'color','r');
[xclick1, yclick1]=ginput(1);
yclick1 = ysign*yclick1;
delete(h);
answer = [];
while isempty(answer),
      answer = inputdlg({'Enter the value for the X-point just 
clicked'},'',1,{''});
      xin1 = str2num(answer{1});
end;

h = text(0,0,'click any point along the RIGHT half of the 
X-(horizontal) axis');
set(h,'fontsize',16); set(h,'color','r');
[xclick2, yclick2]=ginput(1);
yclick2 = ysign*yclick2;
delete(h);
answer = [];
while isempty(answer),
      answer = inputdlg({'Enter the value for the X-point just 
clicked'},'',1,{''});
      xin2 = str2num(answer{1});
end;
alphax = atan2((yclick2-yclick1),(xclick2-xclick1));

% Calibrate Y-axis
h = text(0,0,'click any point along the LOWER half of the Y-(vertical) axis');
set(h,'fontsize',16); set(h,'color','r');
[xclick3, yclick3]=ginput(1);
yclick3 = ysign*yclick3;
delete(h);
answer = [];
while isempty(answer),
      answer = inputdlg({'Enter the value for the Y-point just 
clicked'},'',1,{''});
      yin3 = str2num(answer{1});
end;

h = text(0,0,'click any point along the UPPER half of the Y-(vertical) axis');
set(h,'fontsize',16); set(h,'color','r');
[xclick4, yclick4]=ginput(1);
yclick4 = ysign*yclick4;
delete(h);
answer = [];
while isempty(answer),
      answer = inputdlg({'Enter the value for the Y-point just 
clicked'},'',1,{''});
      yin4 = str2num(answer{1});
end;
alphay = atan2((yclick4-yclick3),(xclick4-xclick3));

% QC to see if axes are orthogonal
if( abs(alphax - alphay)-pi/2 > .02),
     warndlg('Axes not orthogonal. Start program again');
     return;
end;
alphaav = (alphax+alphay-pi/2)/2;

figure(himage);
h = text(0,0,{'Left-Click to Digitize Points';
               'Right-Click to Undo';
               'RETURN when finished'});
set(h,'fontsize',14,'color','r');
xdigt = 1;
k=0;
while ~isempty(xdigt),
     [xdigt,ydigt,button]=ginput(1);
     if ~isempty(xdigt),
         if button == 1,
             k = k+1;
             hdig(k) = plot(xdigt,ydigt,'+'); set(hdig(k),'markersize',14);
             xdig(k) = xdigt;
             ydig(k) = ysign*ydigt;
         else,
             delete(hdig(k));
             xdig(k) = [];
             ydig(k) = [];
             k = k-1;
         end;
     end;
end;
delete(h);

% Scale digitized points

% Define the intersection of the two axes, defined by the four clicks,
% as the origin, even though their values are not zero

% X-axis is a line y=a1*x+b1  in the pixel space; Y-axis is a line x = a2*y+b2;
a1 = (yclick2-yclick1)/(xclick2-xclick1);
a2 = (xclick4-xclick3)/(yclick4-yclick3);
b1 = (xclick2*yclick1-xclick1*yclick2)/(xclick2-xclick1);
b2 = (yclick4*xclick3-yclick3*xclick4)/(yclick4-yclick3);

% point of intersection of these two, in pixel space
y0 = (a1*b2+b1)/(1-a1*a2);
x0 = a2*y0 + b2;

% coordinates of the data at this defined original
xin0 = xin1 + (x0-xclick1)*(xin1-xin2)/(xclick1-xclick2);
yin0 = yin3 + (y0-yclick3)*(yin3-yin4)/(yclick3-yclick4);

% translate the digitized data to this new origin, without rotation
xdigtemp = xdig - x0;
ydigtemp = ydig - y0;

% now rotate the axes

xdig =  cos(alphaav)*xdigtemp + sin(alphaav)*ydigtemp;
ydig = -sin(alphaav)*xdigtemp + cos(alphaav)*ydigtemp;

% now scale these points to the desired data space
xfact = (xin2-xin1)/sqrt((xclick2-xclick1)^2+(yclick2-yclick1)^2);
yfact = (yin4-yin3)/sqrt((xclick4-xclick3)^2+(yclick4-yclick3)^2);

x = xdig*xfact + xin0;
y = ydig*yfact + yin0;
figure; plot(x,y,'o');


