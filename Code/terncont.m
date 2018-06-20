function [c,h,output]=terncont(data,z,varargin);
%FUNCTION TERNCONT(DATA,Z,VARARGIN)
% Ternary contour plot. 
% 
% DATA: Nx3 matrix. Three columns consisting of fractions (non-negative) of 
% 1st, (Top), 2nd (Lower Right), and 3rd (Lower Left) constituent. 
% Z: Property to be contoured at each composition. Nx1 vector 
% Fractions will be normalized if sum of each row is not equal to 1.
%
% VARARGIN: Optional input parameter to specify attributes of contour plot
%           (see CONTOUR). 
%	    Ex: terncont(data,15,'-r') draws 15 contour lines in red.
% [C, H]=TERNCONT(...)  returns contour matrix C as described in
%    CONTOURC and a column vector H of handles to LINE or PATCH
%    objects, one handle per line.  Both of these can be used as
%    input to CLABEL.
%
% See also: TERNARY, CONTOUR

%	Written by Isao Takahashi, Tapan Mukerji 11/19/99, 7/2000 

if any(data(:)<0), error('negative fractions!'); end;

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

axis equal off; 
axis([-1/r3 1/r3 -1/3 2/3])

xi=linspace(-1/r3,1/r3,100); yi=linspace(-1/3,2/3,100)';
zi=griddata(output(:,1),output(:,2),z,xi,yi);

[c,h]=contour(xi,yi,zi,varargin{:});
dx=.05;
ht(1)=text(0,2/3+dx,'Frac 1');
ht(2)=text(1/r3,-1/3-dx,'Frac 2');
ht(3)=text(-1/r3,-1/3-dx,'Frac 3');
set(ht,'horizontalalignment','center')
