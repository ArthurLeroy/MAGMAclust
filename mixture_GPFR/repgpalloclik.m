function [fX, dfX] = repgpalloclik(X, alltraindata, B, Vm, Um, sp)
%
%Used by GP allocation model
%comnum: number of components in mixture
%


comnum = sp(1);
gppno = sp(2);
nord = sp(3);
knots = sp(4:end); 

m = size(alltraindata,1);  %number of all training curves

fXt = zeros(m,comnum);  %
dfXt = cell(m,comnum);  %
gplikfun = ['gpm0', num2str(gppno), 'lik(Xk, input, target, mu)'];

for i = 1:m    
    tdata = alltraindata{i}(:,1);
    input = alltraindata{i}(:,2:end-1);
    target = alltraindata{i}(:,end);
    tdsize = size(input,1);
    gpparasize = size(input,2)*gppno + 2;    %w, (a), v1, v0
    
    PHI = bsplineM(tdata, knots, nord);
    
    for k=1:comnum
        Xk = X((k-1)*gpparasize+1 : k*gpparasize);        
        mu = PHI*B{k}*Um(:,i);
        %[a, da] = gpm01lik(Xk, input, target, mu);
        [a, da] = eval(gplikfun);
        fXt(i,k) = a;
        dfXt{i,k} = da;        
    end
end


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

    
TempSum = zeros(m,1);
for i=1:m
    TempSum(i) = exp(-fXt(i,:))* pimk(i,:)';
end


fX = sum(log(TempSum));


dfX = zeros(comnum*gpparasize + (comnum-1)*infonum,1);
for k=1:comnum
    for i=1:m
        dfX((k-1)*gpparasize+1 : k*gpparasize) = dfX((k-1)*gpparasize+1 : k*gpparasize)...
            - pimk(i,k)*exp(-fXt(i,k))*dfXt{i,k}/TempSum(i);
    end
end
if comnum >1        % To allow the case with only one component
    for k=1:comnum-1
        for i=1:m
            dfX(comnum*gpparasize + (k-1)*infonum+1 : comnum*gpparasize + k*infonum)...
                = dfX(comnum*gpparasize + (k-1)*infonum+1 : comnum*gpparasize + k*infonum)...
                + exp(-fXt(i,k))*pimk(i,k)*Vm(i,:)' / TempSum(i) - pimk(i,k)*Vm(i,:)';
        end
    end
end

fX = -fX;
dfX = -dfX;
