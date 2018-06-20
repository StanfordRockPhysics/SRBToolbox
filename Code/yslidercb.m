function yslidercb
hyslider=findobj(gcf,'tag','yslider'); hytext=findobj(gcf,'tag','ytext');
hax3d=findobj(gcf,'tag','ax3d');hyimg=findobj(gcf,'tag','yimg');
hyl=findobj(gcf,'tag','syline');
sy=get(hyslider,'value');
set(hytext,'string',['y= ' num2str(sy)]);

data=get(hax3d,'userdata');
x=data{1}; y=data{2}; z=data{3}; v=data{4};
syindex=sum(sy>=y); if syindex<=0, syindex=1; end;
set(hytext,'string',['y= ' num2str(y(syindex))]);
vslice=squeeze(v(syindex,:,:)).'; 
set(hyimg,'cdata',vslice);
set(hyl,'ydata',[sy sy sy ]);
