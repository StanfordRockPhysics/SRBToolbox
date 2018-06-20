function xslidercb
hxslider=findobj(gcf,'tag','xslider'); hxtext=findobj(gcf,'tag','xtext');
hax3d=findobj(gcf,'tag','ax3d');hximg=findobj(gcf,'tag','ximg');
hxl=findobj(gcf,'tag','sxline');
sx=get(hxslider,'value'); 

data=get(hax3d,'userdata');
x=data{1}; y=data{2}; z=data{3}; v=data{4};
sxindex=sum(sx>=x); if sxindex<=0, sxindex=1; end;
set(hxtext,'string',['x= ' num2str(x(sxindex))]);
vslice=squeeze(v(:,sxindex,:)).'; 
set(hximg,'cdata',vslice);
set(hxl,'xdata',[sx sx sx ]);
