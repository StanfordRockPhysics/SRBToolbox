function out = kirchmig(half, t0, dt, dx, vrms, in)
% OUT = KIRCHMIG(HALF, T0, DT, DX, VRMS, IN)
% Kirchhoff migration with nearest neighborhood.
% IN, OUT : input and output data dim(IN)=dim(OUT)=[nt,nx]
% HALF    : =1  for half derivative operator, =0, no operator
%          (corrects the waveform of Kirchhoff operator)
% TO      : start time
% DT, DX  : time and spacing intervals
% VRMS    : rms velocity, length(vrms)==nt or 1 (const. vel.)
%
%	e.g.
%               >> data=zeros(50,50);
%               >> data(24:26,24:26)=1;
%               >> imagesc(data)
%               >> model=kirchadj(0,0,0,5e-3,10,2300,data);
%               >> imagesc(model)
%		>> migr=kirchmig(0,0,5e-3,10,2300,model);  
%		>> imagesc(migr)
%
% KIRCHMIG is vectorized code and is an order of magnitude faster than
% migration with KIRCHADJ.
%
% See also: KIRCHADJ, GAZADJ, STOLTMOD, STOLTMIG


% 1999. 12.  Youngseuk Keehm
% Reference: Claerbout, Basic Earth Imaging
%
% Vectorized code design, T. Mukerji

if (half == 1)
   in=halfdifa(1,in);                    % half derivative operator
end
[nt, nx]=size(in);   nz=nt;
out=zeros(nt,nx);
if(length(vrms)==1)                      % constant velocity model
   vrms=vrms*ones(nt,1); 
end
in=[zeros(nt,nx-1),in,zeros(nt,nx-1)];   % zero padding w.r.t. x-axis
xind =[1:2*nx-1];	                        % index of x-grid
xbin =[-nx+1:nx-1]*dx;	                        

for iz=2:nz
   z=(iz-1)*dt;                               
   t=t0 + sqrt( z*z+xbin.*xbin*4/(vrms(iz)^2) );    % hyperbola equation
   it= fix(1.5 + t/dt);                             % nearest neighbors
   amp=(z./t).*sqrt(dt*nt./t);	                     % amplitude correction		
   it(it>nz)=0; id=find(it);                        % index of time
   
   index=sub2ind([nt,(3*nx-2)],it(id),xind(id))';   % set of index (it,ix) for hyperbola
   indm=repmat(index,[1,nx])+repmat(nz*[0:nx-1],[length(index),1]);     % index matrix
   out(iz,:)=sum(in(indm).*repmat(amp(id)',[1,nx]));
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
%
[nt,nx]=size(cv);
cv=ifft(cv);
w=[0:nt-1]'.*2*pi/nt;
cz=sqrt( 1-exp( complex(0,w) ) );
if(adj==1) cz=conj(cz); end
czm=repmat(cz,1,nx);
cv=cv.*czm;
cv=real(fft(cv));
