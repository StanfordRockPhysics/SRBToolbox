function [bayesresult]=pdfstat(currentdata,li,pfac,attcomb,error);

tempsize=size(currentdata);
if ndims(currentdata)<4,
	currentdata=reshape(currentdata,[tempsize(1:end-1) ones(1,4-length(tempsize)) tempsize(end)]);
end
tempdims=sum(tempsize(1:end-1)~=1);
if nargin<=3, attcomb=[1:tempdims]; 
elseif prod(attcomb)==0, attcomb=[1:tempdims];
end
if nargin<=2, pfac=ones(tempsize(end)-1,1)/(tempsize(end)-1); end;
if nargin<=1, li=1:(tempsize(end)-1); end;

fid=fopen('tempstat.txt','w');
fprintf(fid,'\nBasic Statistics\n');

switch tempdims
case 1, comb={1};
case 2, comb={[1 2],	[1;2]};
case 3, comb={[1 2 3], [1 2 3;2 3 1],[1;2;3]};
end

if nargin==5
fprintf(fid,'\n');
fprintf(fid,['Measurement Error\n']);
fprintf(fid,['Parameters      :',repmat('%10i',[1 length(attcomb)]),'\n'],attcomb);
fprintf(fid,['Error           :',repmat('%10.3f',[1 length(attcomb)]),'\n'],error(1,:)); 
fprintf(fid,['Error (percent) :',repmat('%10.3f',[1 length(attcomb)]),'\n'],error(2,:)*100);
fprintf(fid,'\n');
end

for k=1:tempdims
	fprintf(fid,'\n');
	fprintf(fid,['%3i Parameter Combinations'],k);
	fprintf(fid,'\n');

	for j=1:size(comb{k},2)
		dim=comb{k}(:,j)';
		junkdim=setdiff([1 2 3],dim);
		junkdata=currentdata;

		for l=1:length(junkdim)
		junkdata=sum(junkdata,junkdim(l));
		end

		dim=attcomb(dim);

		berror=bayes(junkdata,li,pfac);
		Berror{k}{j}=bayes(junkdata,li,pfac);
		entro=centropy(junkdata,pfac);
		n=size(berror,1);

		TotalError{k}{j}=1-sum(diag(berror));
		SuccessRate{k}{j}=sum(diag(berror));
		fprintf(fid,'\n');
		fprintf(fid,['Parameters :',repmat('%3i',[1 length(dim)]),'\n'],dim);
		fprintf(fid,'\n');
		fprintf(fid,'Bayes Errors\n');
		fprintf(fid,'\n');
		fprintf(fid,'Total: %5.3f, Success Rate: %5.3f\n',1-sum(diag(berror)),sum(diag(berror)));
		fprintf(fid,'\n');
		fprintf(fid,'Unconditioned Probabilities (p{prediction,true})\n');
		fprintf(fid,['True facies     :    :Code:',repmat('%7i',[1 n]),'\n'],li');
		for i=1:n
		fprintf(fid,['Predicted facies: #%2i: %3i:',repmat('  %5.3f',[1 n]),'\n'],i,li(i),berror(i,:));
		end
		fprintf(fid,'\n');
		fprintf(fid,'Conditioned Probabilities (p{true|prediction})\n');
		fprintf(fid,['True facies     :    :Code:',repmat('%7i',[1 n]),'\n'],li');
		for i=1:n
		junkberror=sum(berror(i,:));
			if junkberror==0, tempberror=berror(i,:)*0; condberror(i,:)=berror(i,:)*0;
			else, tempberror=berror(i,:)./junkberror; 
   			condberror(i,:)=berror(i,:)./junkberror; %this is defined for results windows
			end;
		fprintf(fid,['Predicted facies: #%2i: %3i:',repmat('  %5.3f',[1 size(berror,2)]),'\n'],i,li(i),tempberror);
		end
		Condberror{k}{j}=condberror;
		H{k}{j}=entro(4)';
		Hcond{k}{j}=entro(5)';
		Icond{k}{j}=entro(3)';
		HcondH{k}{j}=entro(6)';
		fprintf(fid,'\n');
		fprintf(fid,'\n');
		fprintf(fid,'Shannon''s Information\n');
		fprintf(fid,'\n');
		fprintf(fid,'       H(facies)              : %5.2f\n',entro(4)');
		fprintf(fid,'       H(facies|attribute)    : %5.2f\n',entro(5)');
		fprintf(fid,'       I(facies|attribute)    : %5.2f\n',entro(3)');
		fprintf(fid,'H(facies|attribute)/H(facies) : %5.2f\n',entro(6)');
		fprintf(fid,'\n');
		fprintf(fid,'\n');
		fprintf(fid,'\n');
	end
end
fclose(fid);

bayesresult.total=TotalError;
bayesresult.successrate=SuccessRate;
bayesresult.prob=Berror;
bayesresult.condprob=Condberror;
bayesresult.entro=H;
bayesresult.condentro=Hcond;
bayesresult.condinfo=Icond;
bayesresult.condentronorm=HcondH;
 
