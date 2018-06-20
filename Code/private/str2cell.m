function out=str2cell(in);

if isstruct(in)
in=in(:);
names=fieldnames(in);
fields=struct2cell(in);
	for i=1:length(in)
	temp={names{:};fields{:,i}};
	out(i,:)=temp(:)';
	end;
else
out={};
end;
