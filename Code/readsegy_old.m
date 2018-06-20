function [D,H] = readsegy(FID,hw,min,max)
%READSEGY  READSEGY(FID,hw,min,max) returns the data D and
%          the trace header H. The data and headers are extracted 
%          in the  range of the hw given by min and max.
%  FID = fopen('filename','r','b');  % when running on a Linux based PC
%  FID = fopen('filename','r');      % when running on a SUN system
%  FID = fopen('filename','r','l');  % when reading PC segy file (byte-flipped)
%
%   example:    [D,H] = readsegy(FID,'cdp',500,550) will provide
%               the traces and associated headers for traces with the 
%               header word cdp goes from 500 to 550.
%
%   example:    [D,H] = readsegy(FID,'offset',250,510) like the above
%               example but now the header word 'offset' is used to
%               read traces with offsets in the range 250-500.
%
%   example:    [D,H] = readsegy(FID) reads everything until end
%               of file.
%
% Requires file to be in IEEE and NOT IBM float format. 
% The segy format has no way of indicating which floating point format
% the file is in, IBM or IEEE. It only indicates whether it is floating
% point or not. Seismic data in the wrong floating format can be deceptive
% since both formats will appear like reasonable traces, with only subtle 
% differences. An IBM float data file read as though it was IEEE float
% may be detected by ploting and zooming in on a single trace. If the 
% peaks appear flattened or clipped, then the file was IBM format and should 
% be converted to IEEE.
%
%See also FOPEN, FREAD


%  M.D.Sacchi, July 1997, Dept. of Physics, UofA.
%  sacchi@phys.ualberta.ca
% 
%  Modified by T. Mukerji


hdrblock=3600;                       % EBCDIC header (3200) + line header (400)
 

   load segy.mat                       % load the definitions of 
                                       % the header words

   load countsegy.mat                      % load the position of 
                                       % each word in the header (in bytes)

status = fseek(FID,3224,'bof'); 
fmt = fread(FID,1,'int16');            % read data format from binary header
switch fmt
 case 1				       % 4 bytes float
   segy.trace='float';
   nsfactor=1;
 case 2				       % 4 bytes integer
   segy.trace='int';
   nsfactor=1;
 case 3				       % 2 bytes integer
   segy.trace='int16';
   nsfactor=2;
 case 5				       % 1 byte integer
   segy.trace='int8'; 
   nsfactor=4;
 otherwise
   disp('Data format code unrecognised. Assuming 4 bytes float.');
   segy.trace='float';
   nsfactor=1; fmt=1;
end;

   status = fseek(FID,hdrblock+countsegy.ns,'bof'); % go to the beggining of file

   ns = fread(FID,1,segy.ns);          % read ns from first trace
                                       % ns is the number of samples per trace

   total = 60+ns/nsfactor;                      % total nuber of 4bytes words


   max_traces=9999999;                 % maximum number of traces (will
                                       % stop before). The variable status
                                       % will make the code stop when
                                       % the end of file is reached  

 if nargin>2;
   hc=eval(strcat('countsegy.',hw));       % assigned the header word required  
   hp=eval(strcat('segy.',hw));        % to extract the traces.
 j = 1;                                % counter 
 for k =1:max_traces
   position = total*(k-1)*4+hc;        % where in bytes is the header word
   status = fseek(FID,hdrblock+position,'bof'); 
  if status == 0                       % stop when end of file is reached
    w = fread(FID,1,hp);
     if  w>=min                        % pick traces with a given range
      if w<=max                        % of the desired header word
       position = total*(k-1)*4+countsegy.trace;
       status = fseek(FID,hdrblock+position,'bof'); 
       trace = fread(FID,ns,segy.trace);
       j = j + 1;  
       D(:,j-1)  = trace(:);           % load traces into D
       H(j-1)  = header(FID,ns,k,fmt);     % load each header in a structure H
%       disp(j-1)
      end
     end
   else
  return
 end
 end 

 else 

% when no hw and bounds are give, reads evrything

 for k =1:max_traces
       position = total*(k-1)*4+countsegy.trace;
       status = fseek(FID,hdrblock+position,'bof'); 
       if status == 0                 
       trace = fread(FID,ns,segy.trace);
       if length(trace) ~= ns, return, end;
       D(:,k)  = trace(:);           % load traces into D
       H(k)  = header(FID,ns,k,fmt);     % load each header in a structure H
       if rem(k,100)==0, disp(k), end;
   else
  return
 end
 end 
 end 







 [message,errnum] = ferror(FID)
 fclose(FID);


   

