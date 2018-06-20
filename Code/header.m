function header=header(FID,ns,k,fmt)
%HEADER    HEADER(FID,ns,k) returns the complete trace header
%          into a structure called header. k is the sequential trace
%          number in the segy file.
%          ns the is number of samples per trace.
%          fmt is the data format code:
%               1  4 bytes float
%               2  4 bytes integer
%               3  2 bytes integer
%               5  1 byte integer
%
%   Example: header = header(FID,1024,10,1) gets the complete header
%                     of the trace 10 in the file FID (already opened with 
%                     FID=fopen(...)). The trace has 1024 samples. The latter
%                     can be got using the program 'extract.m' 

%
%  M.D.Sacchi, July 1997, Dept. of Physics, UofA.
%  sacchi@phys.ualberta.ca
%
% Modified by T. Mukerji

hdrblock=3600;                       % EBCDIC header (3200) + line header (400)
switch fmt
 case 1                                % 4 bytes float
   nsfactor=1;
 case 2                                % 4 bytes integer
   nsfactor=1;
 case 3                                % 2 bytes integer
   nsfactor=2;
 case 5                                % 1 byte integer
   nsfactor=4;
end;
 
   offset = (240+ns*4/nsfactor)*(k-1);
   status=fseek(FID,hdrblock+offset,'bof');
   header.tracl=fread(FID,1,'int');    
   header.tracr=fread(FID,1,'int');     
   header.fldr=fread(FID,1,'int');     
   header.tracf=fread(FID,1,'int'); 
   header.ep=fread(FID,1,'int');       
   header.cdp=fread(FID,1,'int');     
   header.cdpt=fread(FID,1,'int');    
   header.trid=fread(FID,1,'short'); 
   header.nva=fread(FID,1,'short');    
   header.nhs=fread(FID,1,'short');   
   header.duse=fread(FID,1,'short');   
   header.offset=fread(FID,1,'int');   
   header.gelev=fread(FID,1,'int');  
   header.selev=fread(FID,1,'int');    
   header.sdepth=fread(FID,1,'int');   
   header.gdel=fread(FID,1,'int');     
   header.sdel=fread(FID,1,'int');     
   header.swdep=fread(FID,1,'int');   
   header.gwdep=fread(FID,1,'int'); 
   header.scalel=fread(FID,1,'short'); 
   header.scalco=fread(FID,1,'short'); 
   header.sx=fread(FID,1,'int');       
   header.sy=fread(FID,1,'int');       
   header.gx=fread(FID,1,'int');       
   header.gy=fread(FID,1,'int');       
   header.counit=fread(FID,1,'short'); 
   header.wevel=fread(FID,1,'short');  
   header.swevel=fread(FID,1,'short'); 
   header.sut=fread(FID,1,'short');    
   header.gut=fread(FID,1,'short');    
   header.sstat=fread(FID,1,'short');  
   header.gstat=fread(FID,1,'short');  
   header.tstat=fread(FID,1,'short');  
   header.laga=fread(FID,1,'short');   
   header.lagb=fread(FID,1,'short');   
   header.delrt=fread(FID,1,'short');  
   header.muts=fread(FID,1,'short');   
   header.mute=fread(FID,1,'short');   
   header.ns=fread(FID,1,'unsigned short');  
   header.dt=fread(FID,1,'unsigned short');  
   header.gain=fread(FID,1,'short');   
   header.igc=fread(FID,1,'short');    
   header.igi=fread(FID,1,'short');   
   header.corr=fread(FID,1,'short');   
   header.sfs=fread(FID,1,'short');    
   header.sfe=fread(FID,1,'short');    
   header.slen=fread(FID,1,'short');  
   header.styp=fread(FID,1,'short');  
   header.stas=fread(FID,1,'short');   
   header.stae=fread(FID,1,'short');  
   header.tatyp=fread(FID,1,'short');  
   header.afilf=fread(FID,1,'short');  
   header.afils=fread(FID,1,'short');  
   header.nofilf=fread(FID,1,'short'); 
   header.nofils=fread(FID,1,'short'); 
   header.lcf=fread(FID,1,'short');    
   header.hcf=fread(FID,1,'short');   
   header.lcs=fread(FID,1,'short');   
   header.hcs=fread(FID,1,'short');   
   header.year=fread(FID,1,'short');   
   header.day=fread(FID,1,'short');    
   header.hour=fread(FID,1,'short');   
   header.minute=fread(FID,1,'short'); 
   header.sec=fread(FID,1,'short');    
   header.timbas=fread(FID,1,'short'); 
   header.trwf=fread(FID,1,'short');   
   header.grnors=fread(FID,1,'short'); 
   header.grnofr=fread(FID,1,'short'); 
   header.grnlof=fread(FID,1,'short'); 
   header.gaps=fread(FID,1,'short');   
   header.otrav=fread(FID,1,'short');  
   header.d1=fread(FID,1,'float');     
   header.f1=fread(FID,1,'float');     
   header.d2=fread(FID,1,'float');     
   header.f2=fread(FID,1,'float');     
   header.ungpow=fread(FID,1,'float'); 
   header.unscale=fread(FID,1,'float');
   header.ntr=fread(FID,1,'int');     
   header.mark=fread(FID,1,'short');  
   header.unass=fread(FID,15,'short');

