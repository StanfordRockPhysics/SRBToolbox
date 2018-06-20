function cbtitle(titstr,fs)
%   cbtitle(titstr,fs)
%
% Puts title on a colorbar in the current figure
%
% titstr     is the desired title string
% fs         is the fontsize

% G Mavko

if nargin<2, fs=10; end;

set(get(findobj(gcf,'tag','colorbar'),'title'),'string',titstr,'fontsize',fs);
set(get(findobj(gcf,'tag','Colorbar'),'title'),'string',titstr,'fontsize',fs);
