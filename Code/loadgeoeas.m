function [out, colnames, line1]=loadgeoeas(file)
%[OUT, COLNAMES, LINE1]=LOADGEOEAS(FILE)
% Load GEO-EAS formatted file (GSLIB)
%
% FILE: Input filename of file in GEO-EAS format.
%       FILE should be character array within single quotes 'file'.
%       For details of GEO-EAS format,see page 21 of 
%       "GSLIB: Geostatistical Software Library and User's Guide",
%       by Deutsch and Journel, 1998.
%       First line is title/heading. Second line has the number 
%       of columns, ncol.
%       Following 'ncol' lines are names for the ncol columns. 
%       Finally the data, with 'ncol' numbers on each line.
% OUT: the numeric data matrix with 'ncol' columns.
% COLNAMES: names of the columns as specified in file (cell array)
% LINE1:    line 1 of file as a char array in matlab

% Written by Tapan Mukerji, Isao Takahshi, 12/4/1999


if ischar(file)
fid=fopen(file);
line1=fgets(fid);
ncol=str2num(fgets(fid));
for k=1:ncol; colnames{k}=fgets(fid); end;
out=fscanf(fid,'%f');
nrow=length(out)/ncol;
out=reshape(out,ncol,nrow).';
fclose(fid);
end;
