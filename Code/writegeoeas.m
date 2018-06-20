function writegeoeas(file,a,colnames,line1)
%WRITEGEOEAS(FILE,A,COLNAMES,LINE1)
%Write GEO-EAS formatted file (GSLIB)
%
% FILE: Filename of file to be written in GEO-EAS format.
%       FILE should be character array within single quotes 'file'.
%       For details of GEO-EAS format,see page 21 of 
%       "GSLIB: Geostatistical Software Library and User's Guide",
%       by Deutsch and Journel, 1998.
%       First line is title/heading. Second line has the number 
%       of columns, ncol.
%       Following 'ncol' lines are names for the ncol columns. 
%       Finally the data, with 'ncol' numbers on each line.
% A:    the numeric data matrix 
% COLNAMES (optional): names of the columns specified in a cell array
%       of characters. If specified, the number of columns 
%       of A should be equal to length(colnames).
% LINE1 (optional): line 1 of file as a char array in matlab
%
% E.G. >> writegeoeas('gslib.dat',a,{'x','y','por','seis'},'title')
% where a is the four column matrix:
%    1.000    2.000    3.000    4.000
%    5.000    6.000    7.000    8.000
%    ...       ...      ...     ....
%
% writes out the disk file 'gslib.dat' in GEO-EAS format.
% >> !more gslib.dat
%title
%4
%x
%y
%por
%seis
%    1.000    2.000    3.000    4.000
%    5.000    6.000    7.000    8.000
%    ...       ...      ...     ....
%

% Written by Tapan Mukerji, 12/4/1999

[nrow,ncol]=size(a);
fmt=[];
if nargin<4, line1='Line 1'; end;
if nargin<3, 
for k=1:ncol, colnames{k}=['column ', num2str(k)]; end; 
end;

if ischar(file)
fid=fopen(file,'w');
fprintf(fid,'%s\n',line1);
fprintf(fid,'%s\n',num2str(ncol));
for k=1:ncol, fprintf(fid,'%s\n',colnames{k}); fmt=strcat(fmt,' %8.3f'); end;
fmt=strcat(fmt,'\n');
fprintf(fid,fmt,a.');
fclose(fid);
end;
