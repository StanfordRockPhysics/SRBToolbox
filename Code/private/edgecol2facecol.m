function mscatout=edgecol2facecol(mscatin);

mscatout=mscatin;
for i=1:length(mscatout)
j=strmatch('edgecolor',strvcat(mscatout{i}{:}));
mscatout{i}{j}='facecolor';
end

