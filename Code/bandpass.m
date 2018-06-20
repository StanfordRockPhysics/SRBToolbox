function filtx=bandpass(xx,dt,flow,fhigh,n,ps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function filtx=bandpass(xx,dt,flow,fhigh,n,ps)
%
% This function filters an array of data xx
% with sampling rate of dt with a bandpass filter
% xx-data size[nt,nx] 
% flow- lowfrequency (min frequency, Hz)
% fhigh- high frequency (max frequency, Hz)
% n - number of points in filter (optional input)
% ps- phase=0 for zero phase input, (optional input)
% the filter coeficients are calculated with fir1, 
% with n point filter. default n=128, if n not specified on input.
% The filtering is done with fftfilt, this introduces a delay as it
% is not a zero-phase filter. delay is ~ 0.5*n time points. for zero phase
% filter insert the optional parameter ps=0.

% Written by Ran Bachrach 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nt,nx]=size(xx);
if nargin<5, n=128; end;
if nargin<6  ps=1; end;
fN=0.5/dt; %nyquist frequency
f1=flow/fN;
f2=fhigh/fN;
df=1/(dt*nt);
B=fir1(n,[f1,f2]);
if ps==1
filtx=fftfilt(B,xx); %introduces linear phase delay
elseif ps==0
 filtx=filtfilt(B,1,xx); %zero phase filter
end
