function dc = covmatch(d, t)
%function dc = covmatch(d, t)
%Scales the multivariate data d to match the covariance and mean of
%target t. d and t must have same number of columns (e.g. 2 for bivariate
%distributions). 


%Written T. Mukerji 

dz = zscore(d); tz = zscore(t); tmean = mean(t); tstd = std(t);
dcov = cov(dz); tcov=cov(tz);
dc = (sqrtm(tcov)*inv(sqrtm(dcov))*dz')';

dc = dc.*repmat(tstd,length(dc),1) + repmat(tmean, length(dc),1);

if nargout==0
    plot(dc(:,1),dc(:,2),'.k'); hold on; plot(t(:,1),t(:,2),'xg'); 
end
