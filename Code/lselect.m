function dx=lselect(onoff,colr)
%function dx=lselect(onoff,colr)
%Linked selection of data point clusters across multiple figure windows.
%Typical usage might be in conjuction with multiple logs in one window
%and one or more cross-plots in other windows.
%DX, output of index of selected points.
%Inputs: ONOFF, string variable can take values 'on' or 'off'.
%Selected points are highlighted with color (and marker) given by COLR.
%COLR is an optional input, the default being '.r'.
%lselect off deselects all previously selected and highlighted points.
%lselect on starts the mouse based selection in the current figure window.
%Use the mouse to select a polygonal region enclosing the points of interest.
%
%See also GINPUT, GETPTS, GETRECT, GETLINE

%Written by T. Mukerji

if nargin<2, colr='.r'; end;
hfig=findobj(0,'type','figure');

switch onoff
  case 'off'
for k=1:length(hfig)
hsp=findobj(hfig(k),'tag','selectedpts'); delete(hsp);
end;

 case 'on'
hofig=setdiff(hfig,gcf); hca=gca;
[x y]=getline;
haline=findobj(gca,'type','line'); 
hspline=findobj(gca,'tag','selectedpts');
hline=setdiff(haline,hspline);
for k=1:length(hline)
cxdata=get(hline(k),'xdata'); cydata=get(hline(k),'ydata');
dx=inpolygon(cxdata,cydata,x,y);
%hold on;h=plot(cxdata(dx),cydata(dx),colr);set(h,'tag','selectedpts');
hofig=hfig;
  for kk=1:length(hofig)
  hax=findobj(hofig(kk),'type','axes');
    for j=1:length(hax)
     hoaline=findobj(hax(j),'type','line');
     hospline=findobj(hax(j),'tag','selectedpts');
     holine=setdiff(hoaline,hospline);
       for jj=1:length(holine)
        ooxdata=get(holine(jj),'xdata');
        if length(ooxdata)==length(cxdata), 
         oxdata=ooxdata; oydata=get(holine(jj),'ydata'); 
         ozdata=get(holine(jj),'zdata');
        end %end if length(ooxdata)==length(cxdata)
        axes(hax(j)); hold on; 
         if isempty(ozdata)
         h=plot(oxdata(dx),oydata(dx),colr);set(h,'tag','selectedpts');
         else
         h=plot3(oxdata(dx),oydata(dx),ozdata(dx),colr);
         set(h,'tag','selectedpts');
         end  %end if isempty(ozdata)
        end; %end for jj=1:length(holine)
       end; %end for j=1:length(hax)
    end; %end for kk=1:length(hofig)
end; %end for k=1:length(hline)
axes(hca);
end; %end switch
