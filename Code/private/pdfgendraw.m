function [smcpdf,vec,handle,err]=pdfgendraw(indata,code,obs,nbin,nfilt,plotop,CONTROL_PROPERTY,SCAT_PROPERTY);

[ucode,i,j]=unique(code(:,1));
nlist=length(ucode);
for k=1:length(ucode);pfac(k)=sum(code(code(:,1)==ucode(k),2))/sum(code(:,2));
end
   
ndim=size(indata,2);
if length(obs)==1, obs=ones(ndim,1)*obs; end
if length(nfilt)==1, nfilt=ones(1,ndim)*nfilt; end
ref=indata;

thresh =repmat([-1 1],[ndim,1])*inf;
scale=ones(ndim,1);
[dataref]=indata;

if ~iscell(nbin)
	for i=1:ndim
		xx{i}=linspace(max(min(dataref(:,i)),thresh(i,1)),min(max(dataref(:,i)),thresh(i,2)),nbin)';
	end
else xx=nbin;
end
[filt,vec,err,S]=filterdefine(dataref,xx,obs,nfilt);

[ocpdf]=cpdf(indata,vec,nlist,code,pfac);

for i=1:nlist;
	tempfilt=filterdefine(indata(code(:,1)==i,:),vec,obs,nfilt);
	if min(size(tempfilt))==1, tempfilt=tempfilt(:);end;
	temppdf=convn(ocpdf(:,:,:,i),tempfilt,'same');
	tempsmcpdf=pdfshape(temppdf,vec,thresh);
	smcpdf(:,:,:,i)=tempsmcpdf./sum(tempsmcpdf(:));
end;
junkvec=size(smcpdf); if length(junkvec)==2,  junkvec(3)=1; end;
smcpdf(:,:,:,nlist+1)=sum(smcpdf.*repmat(reshape(pfac,[1 1 1 length(ucode)]),[junkvec(1:3) 1]),4);

multiscat_property=colormarkerset(CONTROL_PROPERTY,SCAT_PROPERTY,nlist); %close(gcf);

%%%hold on;
if plotop==1
switch ndim
	case 1
	multiscat_property=edgecol2col(multiscat_property);
	figure;
	for i=1:nlist;
		handle1(i,1)=plot(vec{1},smcpdf(:,:,:,i)'); hold on;
		set(handle1(i,1),multiscat_property{i}{:});
	end
	handle=handle1;   

	case 2
	p1=.1; p2=.6; p3=.1; p4=.05;
	figure;
	h1=subplot('position',[p1 p1 p2 p2]); 
	h2=subplot('position',[p1 p1+p2+p4  p2 p3]); title('Contour Plot and Projected Marginal pdfs ','fontsize',16);
	h3=subplot('position',[p1+p2+p4 p1  p3 p2]);subplot(h1)
	hold on;
	if nlist~=1, ncont=5; else ncont=10;	end
	xaxis=vec{1};
	yaxis=vec{2};
	for i=1:nlist;
		tempcpdf=squeeze(smcpdf(:,:,:,i));
	  	[c,handle{i}]=contour(xaxis,yaxis,tempcpdf',ncont); hold on;
		set(handle{i},multiscat_property{i}{:});
		handle{i}=handle{i}';
	end;

	subplot(h2)
	hold on;
	multiscat_property2=edgecol2col(multiscat_property);
	for i=1:nlist;
		handle2=plot(xaxis,sum(smcpdf(:,:,:,i),2));
		set(handle2,multiscat_property2{i}{:});
	end;
	axis auto; ax2=axis;
	
	subplot(h3)
	hold on;
	for  i=1:nlist;
		handle3=plot(sum(smcpdf(:,:,:,i),1),yaxis);
		set(handle3,multiscat_property2{i}{:});
	end;
	axis auto; ax3=axis;

	subplot(h1)
	ax1=axis;
	subplot(h2); xlim([ax1(1:2)]); ylim([0 max([ax2(4) ax3(2)])]);
	set(gca,'xticklabel',''); box on;
	subplot(h3); ylim([ax1(3:4)]); xlim([0 max([ax2(4) ax3(2)])]);
	set(gca,'yticklabel',''); box on;

	handle=[handle{:} handle2 handle3]';
	
	case 3
	p1=.1; p2=.6; p3=.1; p4=.05;
	multiscat_property3=edgecol2facecol(multiscat_property);
	xaxis=vec{1};
	yaxis=vec{2};
	zaxis=vec{3};
	[xx,yy,zz]=meshgrid(yaxis,xaxis,zaxis);
	figure('name','surface','tag','surf123','visible','off');
	for i=1:nlist;
		handle(i,1)=patch(isosurface(xx,yy,zz,smcpdf(:,:,:,i))); hold on;
		isonormals(xx,yy,zz,smcpdf(:,:,:,i),handle(i,1));
		set(handle(i,1),multiscat_property3{i}{:},'edgecolor','none');
	end
	%daspect([1 1 .5]);
	view(3)
	camlight(0,0);
	lighting phong

	%pdfs for individual attributes
	plot_property=edgecol2col(multiscat_property);
	figure('name','PDFs for Single Attributes');
	handle_pdf1=subplot(ndim,1,1);
	for i=1:nlist;
	  	tempcpdf=squeeze(sum(smcpdf(:,:,:,i),3));
		handle_curves=plot(xaxis,sum(tempcpdf,2)); hold on;
		set(handle_curves,plot_property{i}{:});
	end
	title('Attribute 1');
	handle_pdf2=subplot(ndim,1,2);
	for i=1:nlist;
	  	tempcpdf=squeeze(sum(smcpdf(:,:,:,i),3));
		handle_curves=plot(yaxis,sum(tempcpdf,1)); hold on;
		set(handle_curves,plot_property{i}{:});
	end
	title('Attribute 2');
	handle_pdf3=subplot(ndim,1,3);
	for i=1:nlist;
	  	tempcpdf=squeeze(sum(smcpdf(:,:,:,i),1));
		handle_curves=plot(zaxis,sum(tempcpdf,1)); hold on;
		set(handle_curves,plot_property{i}{:});
	end
	title('Attribute 3');

	%combinations
	for attcomb=1:3
		switch attcomb
		case 1, t='Attribute Combination: 2&3'; axisa=vec{2};axisb=vec{3};
		case 2, t='Attribute Combination: 1&3'; axisa=vec{1};axisb=vec{3};
		case 3, t='Attribute Combination: 1&2'; axisa=vec{1};axisb=vec{2};
		end
		figure('tag',['attcomb' num2str(attcomb)],'Name',t,'visible','off');
		h1=subplot('position',[p1 p1 p2 p2]); 
		h2=subplot('position',[p1 p1+p2+p4  p2 p3]); title('Contour Plot and Projected Marginal pdfs ','fontsize',16);
		h3=subplot('position',[p1+p2+p4 p1  p3 p2]);subplot(h1)
		hold on;
		if nlist~=1, ncont=5; else ncont=10;	end
		xaxis=vec{1};
		yaxis=vec{2};
        
		for i=1:nlist;
			tempcpdf=squeeze(sum(smcpdf(:,:,:,i),attcomb));
		  	[c,handlecont{i}]=contour(axisa,axisb,tempcpdf',ncont); hold on;
		    set(handlecont{i},multiscat_property{i}{:});
			handlecont{i}=handlecont{i}';
		end;

		subplot(h2)
		hold on;
		multiscat_property2=edgecol2col(multiscat_property);
		for i=1:nlist;
		  	tempcpdf=squeeze(sum(smcpdf(:,:,:,i),attcomb));
			handle2=plot(axisa,sum(tempcpdf,2));
			set(handle2,multiscat_property2{i}{:});
		end;
		axis auto; ax2=axis;

		subplot(h3)
		hold on;
		for i=1:nlist;
			tempcpdf=squeeze(sum(smcpdf(:,:,:,i),attcomb));
			handle3=plot(sum(tempcpdf,1),axisb);
			set(handle3,multiscat_property2{i}{:});
		end;
		axis auto; ax3=axis;

		subplot(h1)
		ax1=axis;
	
		subplot(h2); xlim([ax1(1:2)]); ylim([0 max([ax2(4) ax3(2)])]);
		set(gca,'xticklabel',''); box on;
		subplot(h3); ylim([ax1(3:4)]); xlim([0 max([ax2(4) ax3(2)])]);
		set(gca,'yticklabel',''); box on;

		handle=[handlecont{:} handle2 handle3]';
	end
end
hold off;
else handle=[];
end;


