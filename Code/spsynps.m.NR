function y=spsynps(lambda,n)
% Y = SPSYNPS(LAMBDA,N) generates a 1-D Poisson series Y of length N and 
% Poisson parameter LAMBDA

% written by Philippe Rio

beta=lambda;
         val1=2000;
         val2=4000;
         aa = poissrnd(beta,1,2*n);
         a = zeros(1,n);
         i=1;
         ii=0;
         val=0;
         while ii<n
            if aa(i)~=0
               for j=1:aa(i);
                   a(j+ii)=val;
               end
               ii=ii+aa(i);
               if val==0
                  val=1;
               else
                  val=0;
               end
            else
               if val==0
                  val=1;
               else
                  val=0;
               end
            end
            i=i+1;
         end
         a=a(1:n);
         y=a;
         y=[];
         for i=1:n;
             y=[y;a(i)*ones(1,n)];
         end
         y=val1+(val2-val1)*y;
         imagesc(y,[min(min(y))*0.8 max(max(y))*1.2]);
