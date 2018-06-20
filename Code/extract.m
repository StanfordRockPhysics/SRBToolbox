function [value] = extract(filename,hw)
%EXTRACT EXTRACT(filename,hw) return the numerical value
%        of a SEGY header word. hw is the desired header 
%        word to retrieve. 
%        filename and hw are characters.
%
%   example:    cdp = extract('data','cdp') will provide
%               the cdp numbers of each trace in 'data' 
%               To list the rest of the header words in the
%               SEGY format load segy.mat. This loads a structure
%               call segy where the header words and the associated
%               precission are listed. 
%
%   example:    D = extract('data','trace') will extract the traces
%               into a matrix called D.

%  M.D.Sacchi, July 1997, Dept. of Physics, UofA.
%        
%  sacchi@phys.ualberta.ca
%
 

   FID = fopen(filename,'r','b');      % when running on a Linux based PC

%  FID = fopen(filename,'r');          % when running on a SUN system 

   load segy.mat                       % load the definitions of 
                                       % the header words

   load countsegy.mat                      % load the position of 
                                       % each header word in the file

   status = fseek(FID,countsegy.ns,'bof'); % go to begging of file

   ns = fread(FID,1,segy.ns);          % read ns from first trace

   if strcmp(hw,'ns') == 0;            % get only ns and return
    else
     value = ns;
      return
       end

   total = 60+ns;                      % total num of 4bytes words 
                                       % in the header

   if strcmp(hw,'trace') == 0;              
    elements = 1;
     word=zeros(1,1);
      else
      elements = ns;
       word=zeros(ns,1);
        end

   hc=eval(strcat('countsegy.',hw));
   hp=eval(strcat('segy.',hw));


   max_traces = 999999;

   for j=1:max_traces;

     position = total*(j-1)*4+hc;
      status = fseek(FID,position,'bof'); 

       if status == 0                  % check status and read 
        word = fread(FID,elements,hp);
        j = j + 1  
        value(:,j-1)  = word(:,1);
      else
       return                          % return if end of file reached 
      end 
   end 

 [message,errnum] = ferror(FID)



   

