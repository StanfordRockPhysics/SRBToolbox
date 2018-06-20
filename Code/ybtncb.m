function ybtncb
hyslider=findobj(gcf,'tag','yslider'); hytext=findobj(gcf,'tag','ytext');
hax3d=findobj(gcf,'tag','ax3d');hyimg=findobj(gcf,'tag','yimg');
hyl=findobj(gcf,'tag','syline');
%sx=fix(get(hxslider,'value'));
%set(hxtext,'string',sx);

data=get(hax3d,'userdata');
x=data{1}; y=data{2}; z=data{3}; v=data{4};

for k=1:length(y)
sy=y(k); set(hyslider,'value',sy); 
set(hytext,'string',['y= ' num2str(sy)]);
syindex=sum(sy>=y); if syindex<=0, syindex=1; end;
vslice=squeeze(v(syindex,:,:)).'; 
set(hyimg,'cdata',vslice);
set(hyl,'ydata',[sy sy sy ]);
drawnow;
end;
