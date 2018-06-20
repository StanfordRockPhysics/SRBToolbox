function out = kirchadj(adj, halfdif, t0, dt, dx, vrms, in)
%KIRCHHOFF DIFFRACTION/MIGRATION
%OUT = KIRCHADJ(ADJ, HALFDIF, T0, DT, DX, VRMS, IN)
% ADJ      : 0=diffraction modeling, 1=migration
% HALFDIF  : 0=no, 1=yes (half derivative operator to correct waveform of
%                                                   kirchoff operator)
% T0,DT,DX : start time, time sampling, horizontal spacing
% VRMS     : rms velocity (scalar is OK, if const. velocity). Velocity can
%            vary with depth, in which case length(VRMS)== number of rows of IN.
%            No lateral velocity variations.
% IN,OUT   : input/output data. In modeling mode IN is reflectivity image,
%            OUT is the stack section. In migration mode IN is the
%            unmigrated stack section, OUT is the migrated image.
%
%	     e.g.
%		>> data=zeros(50,50);
%		>> data(24:26,24:26)=1;
%		>> imagesc(data)
%		>> model=kirchadj(0,0,0,5e-3,10,2300,data);
%		>> imagesc(model)
%		>> migr=kirchadj(1,0,0,5e-3,10,2300,model); 
%		>> imagesc(migr) 
%
% See also: KIRCHMIG, GAZADJ, STOLTMIG, STOLTMOD
%           KIRCHMIG is vectorized and is an order of magnitude faster.


% 1999.12. Youngseuk Keehm
% Reference: Claerbout, J., Basic Earth Imaging 

[nt, nx]=size(in);
if(length(vrms)==1) vrms=vrms*ones(nt,1); end
if(adj==1) 	
   data=in; 		modl=zeros(nt,nx); 
   if(halfdif==1)	data=halfdifa(1,data); end
else			
   modl=in;			data=zeros(nt,nx);
end

for ib = -nx:nx
   b=dx*ib;
   for iz = 2:nt
      z = t0 + dt * (iz -1);
      t = sqrt( z*z + (b*2/vrms(iz))^2 );
      it= fix( 1.5 + (t-t0)/dt ); 
      if( it > nt ) break; end;
      amp = (z/t) * sqrt(nt*dt/t);
      for ix = max([1, 1-ib]):min([nx, nx-ib])
         if(adj == 0)
            data(it, ix+ib) = data(it, ix+ib) + modl(iz, ix)*amp;
         else
            modl(iz, ix) = modl(iz, ix) + data(it, ix+ib)*amp;
         end
      end
   end
end
if(adj==0) 
   if(halfdif==1)	data=halfdifa(0,data); end
   out=data;
else		  
   out=modl;
end


function cv=halfdifa(adj, cv)
% HALF DERIVATIVE OPERATOR - Correct the waveform of Kirchhoff operator
%                            according to Huygens' 2nd source theorem.
% ADJ : 0 - for diffraction modeling
%       1 - for Kirchhoff migration
% CV  : input/output  
% 
% Due to the sign convention of Claerbout's book (Basic Earth Imaging),
% Fourier transform here == Inverse transform in MATLAB
%

% 1999.12. Youngseuk Keehm

[nt,nx]=size(cv);
cv=ifft(cv);
w=[0:nt-1]'.*2*pi/nt;
cz=sqrt( 1-exp( complex(0,w) ) );
if(adj==1) cz=conj(cz);end
czm=repmat(cz,1,nx);
cv=cv.*czm;
cv=real(fft(cv));

