function logt=d2t(logd)
%logt=d2t(logd)
%
%Depth-to-time (two-way traveltime) conversion of well logs.
%Input: LOGD, logs in depth, column 1 = depth monotonically increasing, 
%column 2 = velocity, and optional additional columns.
%Output: LOGT column 1 = two-way traveltime, column 2 = velocity
%and other columns interpolated to time axis. Time interval dt for
%new time axis is 0.5*(minimum interval two-way transit time). This prevents
%aliasing. Uses INTERP1Q for fast linear interpolation.

depth=logd(:,1); vel=logd(:,2);

tt=cumsum(2*diff([0;depth])./vel); ttu=[tt(1):0.5*min(diff(tt)):tt(end)].';

logt(:,1)=ttu;

logt(:,2:size(logd,2))=interp1q(tt,logd(:,2:end),ttu);
