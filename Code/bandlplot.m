function bandlplot(in,ltype,lcol,bcol)
%function bandlplot(in,ltype,lcol,bcol)
%in=[depth,val,u,l] nx4 matrix
%optional inputs, ltype - linetype, lcolor - linecolor, bcolor - bandcolor

if nargin==1, ltype='-'; lcol='b'; bcol=[0.8 0.8 0.8]; end;

d=in(:,1); val=in(:,2); u=in(:,3); l=in(:,4);
x=[l;flipud(u)]; y=[d;flipud(d)];
fill(x,y,bcol,'edgecolor',bcol);
hold on;
h=plot(val,d,ltype,'linewidth',1);
set(h,'color',lcol);
set(gca,'ydir','reverse');
hold off;

