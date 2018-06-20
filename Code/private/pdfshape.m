function outpdf = pdfshape(inpdf,xx,thresh)

ndim=length(xx);
outpdf = inpdf;
for i=1:ndim
th=thresh(i,:);
lowflag=(xx{i}<=th(1));
if sum(lowflag)~=0
	boundary=([0;diff(lowflag)]==-1);
	jvec1=ones(1,ndim); jvec1(i)=length(xx{i});
	jvec2=size(outpdf); jvec2(i)=1;
	junk=reshape(lowflag,jvec1);
	lowflagmat=repmat(reshape(lowflag,jvec1),jvec2);
	boundarymat=repmat(reshape(boundary,jvec1),jvec2);
	residual=sum(outpdf.*lowflagmat,i);
	outpdf(boundarymat)=outpdf(boundarymat)+residual';
	outpdf(lowflagmat)=0;
end;
highflag=(xx{i}>=th(2));
if (sum(highflag)~=0)&(1>2)
	boundary=([0;diff(highflag)]==1);	
	jvec1=ones(1,ndim); jvec1(i)=length(xx{i});
	jvec2=size(outpdf); jvec2(i)=1;
	highflagmat=repmat(reshape(highflag,jvec1),jvec2);
	boundarymat=repmat(reshape(boundary,jvec1),jvec2);
	residual=sum(outpdf.*highflagmat,i);
	outpdf(boundarymat)=outpdf(boundarymat)+residual(:);
	outpdf(highflagmat)=0;
end;
end

