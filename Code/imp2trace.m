function s=imp2trace(depth, v, imp, wvlt, dt, nwvlt2)
%  s=imp2trace(depth, v, imp, wvlt, dt, nwvlt2)
%
% Takes input log of impedance, computes reflectivity, and outputs convolved
% seismic trace.   Can be use for normal incidence acoustic impedance, or
% offset trace elastic impedance.
%
% inputs:
%   depth   array of well log depths - single column vector
%   v       array of velocities - can be a vector, 2D, or 3D array
%           must have same dimensions as array of impedances
%   imp     array of impedances - can be a vector, 2D, or 3D array
%   wvlt    wavelet
%   dt      time sample interval for wavelet and output trace
%   nwvlt2  index in wavelet representing the time=zero point.  For
%           example, a Ricker wavelet of length 21 samples would have nwvlt2=11
% outputs
%   s       array of output seismic traces, with time sample interval dt

% Gary Mavko

        nx1 = size(imp,1);
        nx2 = size(imp,2);
        nx3 = size(imp,3);

        for jj=1:nx3, 
            for kk=1:nx2, 
                v(:,kk,jj)   = fillnan(v(:,kk,jj)); 
            end; 
        end;
        
        h      = diff(depth); h(end+1,:,:) = h(end,:,:);   % layer thickness from log
        hh     = repmat(h,[1,nx2,nx3]); 
        deltat = 2*hh./v;                           % 2-way time through each layer, each trace
        ttop   = 2*depth(1)./v(1,:,:);              % 2-way time to top of log, each trace
        t      = repmat(ttop,[nx1,1,1]) + cumsum(deltat);            
        tbot   = t(end,:,:);
        
        dtavg = 0.5*mean(deltat(:,1,1));       % half average travel time through the model layers
        nover = ceil(dt/dtavg);                % oversample factor to take wavelet sample to dtavg
        nover = max(1,nover) ;                 % oversample factor
        time  = [ttop(1,1,1): dt/nover: tbot(1,1,1)]';       % Vector of equally sampled time, finely spaced

        for jj = 1:nx3,
            for kk = 1:nx2,
%             Oversample to equally spaced in time, trying to have each layer sampled
                imp(:,kk,jj) = fillnan(imp(:,kk,jj));
                impt  = interp1(t(:,kk,jj),imp(:,kk,jj),time);
%             Now smooth the impedances and resample back to wavelet interval dt
                if nover > 1,
                    impt  = runavg(impt, nover);
                    impt = impt(1:nover:end);
                end;
%             Estimate reflectivity and convolve reflectivity with wavelet.
                R = diff(impt)./(impt(1:end-1)+impt(2:end)); R(end+1) = 0;
                R(isnan(R)) = 0;
                s(:,kk,jj) = conv(R, wvlt);
            end;
        end;
        s(1:nwvlt2,:,:)=[];

