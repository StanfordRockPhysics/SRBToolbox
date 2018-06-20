function zslidercb
hzslider=findobj(gcf,'tag','zslider'); hztext=findobj(gcf,'tag','ztext');
hax3d=findobj(gcf,'tag','ax3d');hzimg=findobj(gcf,'tag','zimg');
hzl=findobj(gcf,'tag','szline');
sz=get(hzslider,'value');

data=get(hax3d,'userdata');
x=data{1}; y=data{2}; z=data{3}; v=data{4};
szindex=sum(sz>=z); if szindex<=0, szindex=1; end;
set(hztext,'string',['z= ' num2str(z(szindex))]);
vslice=v(:,:,szindex); 
set(hzimg,'cdata',vslice);
set(hzl,'zdata',[sz sz sz ]);
