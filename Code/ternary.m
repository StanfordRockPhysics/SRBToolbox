function ternary(data,varargin);
%FUNCTION TERNARY(DATA,VARARGIN)
% Ternary diagram plot. 
% 
% DATA: Three column matrix consisting of fractions (non-negative) of 
% 1st, (Top), 2nd (Lower Right), and 3rd (Lower Left) constituent. 
% Fractions will be normalized if sum of each row is not equal to 1.
%
% VARARGIN: Optional input parameter to specify attributes of display. 
%	    Ex: ternary(data,'color','red','marker','o')

%	Written by Isao Takahashi. 11/19/99 


data=data./repmat(sum(data,2),[1 3]);

nn=size(data,1);
r3=sqrt(3);
po=[0 2/3;1/r3 -1/3;-1/r3 -1/3];

hold on
ho=findobj(gca,'tag','ternary_diagram_background');

if isempty(ho);
	ho=patch(po(:,1),po(:,2),'w');
	set(ho,'tag','ternary_diagram_background','edgecolor','k');
end;

vec=[0 2/3 1/r3 -1/3 -1/r3 -1/3];
vecrep=repmat(vec,[nn 1]);

output=data(:,[1 1]).*vecrep(:,1:2)+data(:,[2 2]).*vecrep(:,3:4)+data(:,[3 3]).*vecrep(:,5:6);

h=plot(output(:,1),output(:,2),'x');
set(h,'markersize',8,varargin{:});

axis equal 
axis off
axis([-1/r3 1/r3 -1/3 2/3])

dx=.05;
ht(1)=text(0,2/3+dx,'Frac 1');
ht(2)=text(1/r3,-1/3-dx,'Frac 2');
ht(3)=text(-1/r3,-1/3-dx,'Frac 3');

set(ht,'horizontalalignment','center')
