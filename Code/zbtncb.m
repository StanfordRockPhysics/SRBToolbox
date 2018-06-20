function zbtncb
hzslider=findobj(gcf,'tag','zslider'); hztext=findobj(gcf,'tag','ztext');
hax3d=findobj(gcf,'tag','ax3d');hzimg=findobj(gcf,'tag','zimg');
hzl=findobj(gcf,'tag','szline');
%sx=fix(get(hxslider,'value'));
%set(hxtext,'string',sx);

data=get(hax3d,'userdata');
x=data{1}; y=data{2}; z=data{3}; v=data{4};

for k=1:length(z)
sz=z(k); set(hzslider,'value',sz); 
set(hztext,'string',['z= ' num2str(sz)]);
szindex=sum(sz>=z); if szindex<=0, szindex=1; end;
vslice=v(:,:,szindex); 
set(hzimg,'cdata',vslice);
set(hzl,'zdata',[sz sz sz ]);
drawnow;
end;
