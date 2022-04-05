function [fX, dfX] = fun_in_b(Bv, alltraindata, X, otherpara, sp)
%
%

pimk = otherpara{1};
Um = otherpara{2};
invCmk = otherpara{3};

comnum = sp(1);
gppno = sp(2);
nord = sp(3);
knots = sp(4:end); 
nbas = length(knots) + nord - 2;
m = size(alltraindata,1);  %number of all training curves
udim = size(Um,1);

B = cell(comnum,1);
for j=1:comnum
    B{j} = reshape(Bv( (j-1)*nbas*udim+1 : j*nbas*udim ), nbas, udim);
end

fXt = zeros(m,comnum);  %
dfXt = cell(m,comnum);  %
gplikfun = ['gpm0', num2str(gppno), 'lik(Xk, input, target, mu)'];

for i = 1:m    
    tdata = alltraindata{i}(:,1);
    input = alltraindata{i}(:,2:end-1);
    target = alltraindata{i}(:,end);
    tdsize = size(input,1);
    gpparasize = size(input,2)*gppno + 2;    %w, (a,) v1, v0
    
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


TempSum = zeros(m,1);
for i=1:m
    TempSum(i) = exp(-fXt(i,:))* pimk(i,:)';
end


fX = sum(log(TempSum));


dfX = zeros(size(Bv,1),1);
for k=1:comnum
    for i=1:m
        tdata = alltraindata{i}(:,1);
        target = alltraindata{i}(:,end);

        PHI = bsplineM(tdata, knots, nord);
        kUP = kron(Um(:,i),PHI');

        dfX((k-1)*nbas*udim+1 : k*nbas*udim) = dfX((k-1)*nbas*udim+1 : k*nbas*udim)...
            + pimk(i,k)*exp(-fXt(i,k))*kUP*invCmk{i,k}*(target-PHI*B{k}*Um(:,i))/TempSum(i);
    end
end


fX = -fX;
dfX = -dfX;
