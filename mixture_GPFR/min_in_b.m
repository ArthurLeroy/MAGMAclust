function Bmin = min_in_b(X, alltraindata, B0, Vm, Um, para)
%X --> \theta and \beta
%Minimise in B-spline coefficient given previous one as starting value and \theta and \beta
%

%------------------------------------------------------------------------
%Method 2: directly minimising the function fun_of_b in B given \theta and
%\beta and B0 as starting value
%------------------------------------------------------------------------

comnum = para(1);
gppno = para(2);
nord = para(3);
knots = para(4:end); 
nbas = length(knots) + nord - 2;
udim = size(Um,1);
m = size(alltraindata,1);
gpparasize = ( size(alltraindata{1},2) - 2 )*gppno + 2;

%Compute the weights
pimk = zeros(m,comnum);
infonum = size(Vm,2);                %the number of infomations of patients
beta = [X( comnum*gpparasize+1 : end ); zeros(infonum,1)];
for i = 1:m
    for j = 1:comnum
        betak = beta( (j-1)*infonum+1 : j*infonum );
        pimk(i,j) = exp( Vm(i,:)*betak );
    end
    pimk(i,:) = pimk(i,:) / sum(pimk(i,:)) ;
end

%Compute invCmk, which is independent of B
invCmk = cell(m,comnum);
covarfun = ['covfun0', num2str(gppno), '(Xk, input)'];
for k=1:comnum
    Xk = X((k-1)*gpparasize+1 : k*gpparasize);
    for i=1:m
        input = alltraindata{i}(:,2:end-1);
        %Cmk = covfun( Xk, input);
        Cmk = eval(covarfun);
        invCmk{i,k} = inv(Cmk);
    end
end

%Minimise the function
Bv = [];
for j=1:comnum
    Bv = [Bv; reshape(B0{j}, nbas*udim, 1)];
end
otherpara = cell(3,1);
otherpara{1,1} = pimk;
otherpara{2,1} = Um;
otherpara{3,1} = invCmk;


% %Approach 1:
% [Bvmin fxp ip] = minimize(Bv, 'fun_in_b', 10, alltraindata, X, otherpara, para);
% Bmin = cell(comnum,1);
% for j=1:comnum
%     Bmin{j} = reshape(Bvmin( (j-1)*nbas*udim+1 : j*nbas*udim ), nbas, udim);
% end

%Approach 2:
s = 0; smax = 10;  % smax=0 corresponds to 1-step iteration
while 1
    s = s + 1    
	Bmin = computeb(X, alltraindata, B0, otherpara, para);	
    db = 0;
    for j=1:comnum
        db = db + norm(Bmin{j}-B0{j},'fro');
    end
	if (db < 1e-6) | (s>smax)
        break;
    else
        B0 = Bmin;
    end
end


% %-----------------------------------------------------------------------
% %Method 3: express the observations by B-spline, then directly minimising
% %fun_of_b_bspline in B given \theta and \beta and B0 as starting value
% %-----------------------------------------------------------------------
% 
% comnum = para(1);
% nord = para(2);
% knots = para(3:end); 
% nbas = length(knots) + nord - 2;
% udim = size(Um,1);
% m = size(alltraindata,1);
% 
% % Express Y_m by B-spline: Y_m = PHI_m*A_m and compute some constants
% A = cell(m,1);
% PinvCP = cell(m,comnum);
% COEF = cell(m,comnum);
% for j=1:m
%     tdata = alltraindata{j}(:,1);
%     input = alltraindata{j}(:,2:end-1);
%     target = alltraindata{j}(:,end);   
%     gpparasize = size(input,2)*2 + 2;    %w, a, v1, v0
% 
%     for jj=1:size(knots,2)
%         if knots(1,jj)>=tdata(end)
%             mknots = knots(1,1:jj);            
%             break;
%         end
%     end
%     
%     PHI = bsplineM(tdata, mknots, nord);    
%     A{j} = inv(PHI'*PHI)*PHI'*target;
%     for k=1:comnum
%         Xk = X((k-1)*gpparasize+1 : k*gpparasize);
%         Cmk = covfun( Xk, input);
%         invCmk = inv(Cmk);
%         logdetCmk = 2*sum(log(diag(chol(Cmk))));
%         
%         COEF{j,k} = exp( -0.5*size(tdata,1)*log(2*pi) - 0.5*logdetCmk );
%         PinvCP{j,k} = PHI'*invCmk*PHI;
%     end
% end
% 
% %Compute the weights
% pimk = zeros(m,comnum);
% infonum = size(Vm,2);                %the number of infomations of patients
% beta = [X( comnum*gpparasize+1 : end ); zeros(infonum,1)];
% for i = 1:m
%     for j = 1:comnum
%         betak = beta( (j-1)*infonum+1 : j*infonum );
%         pimk(i,j) = exp( Vm(i,:)*betak );
%     end
%     pimk(i,:) = pimk(i,:) / sum(pimk(i,:)) ;
% end
% 
% 
% %Maximise the function
% Bv0 = [];
% for j=1:comnum
%     Bv0 = [Bv0; reshape(B0{j}, nbas*udim, 1)];
% end
% 
% otherpara = cell(6,1);
% otherpara{1} = A;
% otherpara{2} = pimk;
% otherpara{3} = Um;
% otherpara{4} = COEF;
% otherpara{5} = PinvCP;
% sp = [comnum;nbas;udim];
% otherpara{6} = sp;
% 
% options = optimset('Display','iter','LargeScale','off','MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-5,'TolX',1e-6);
% [Bvmin fval] = fminunc(@(Bv) fun_in_b_bspline(Bv,otherpara), Bv0, options);
% 
% Bmin = cell(comnum,1);
% for j=1:comnum
%     Bmin{j} = reshape(Bvmin( (j-1)*nbas*udim+1 : j*nbas*udim ), nbas, udim);
% end


