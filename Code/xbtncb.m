function xbtncb
hxslider=findobj(gcf,'tag','xslider'); hxtext=findobj(gcf,'tag','xtext');
hax3d=findobj(gcf,'tag','ax3d');hximg=findobj(gcf,'tag','ximg');
hxl=findobj(gcf,'tag','sxline');
%sx=fix(get(hxslider,'value'));
%set(hxtext,'string',sx);

data=get(hax3d,'userdata');
x=data{1}; y=data{2}; z=data{3}; v=data{4};

for k=1:length(x)
sx=x(k); set(hxslider,'value',sx); 
set(hxtext,'string',['x= ' num2str(sx)]);
sxindex=sum(sx>=x); if sxindex<=0, sxindex=1; end;
vslice=squeeze(v(:,sxindex,:)).'; 
set(hximg,'cdata',vslice);
set(hxl,'xdata',[sx sx sx ]);
drawnow;
end;
