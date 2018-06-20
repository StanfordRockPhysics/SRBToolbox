function plotxyband(in,ltype,bcol)
%function plotxyband(in,ltype,bcol)
%in=[x,y,u,l], nx4 matrix, x, y, upper and lower range for band
%optional inputs, ltype - linetype, bcolor - bandcolor

if nargin==1, ltype='-b'; bcol=[0.8 0.8 0.8]; end;

hst=ishold;
x=in(:,1); y=in(:,2); u=in(:,3); l=in(:,4);
yy=[l;flipud(u)]; xx=[x;flipud(x)];
fill(xx,yy,bcol,'edgecolor',bcol);
hold on;
h=plot(x,y,ltype,'linewidth',1);
if hst, hold on, else, hold off; end;

