function mscatout=edgecol2col(mscatin);

mscatout=mscatin;
for i=1:length(mscatout)
j=strmatch('edgecolor',strvcat(mscatout{i}{:}));
mscatout{i}{j}='color';
end

