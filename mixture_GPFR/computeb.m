function Bmin = computeb(X, alltraindata, B, otherpara, sp)
%X --> \theta and \beta
%Compute B-spline coefficient given previous one and \theta and \beta
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


R = zeros(m,comnum);
for i=1:m
    fXt(i,:) = fXt(i,:) + repmat(max(-fXt(i,:),[],2),1,comnum); %to avoid 'Not a Number'
    R(i,:) = exp(-fXt(i,:)).* pimk(i,:)/(exp(-fXt(i,:))* pimk(i,:)');
end


Bmin = cell(comnum,1);
for k=1:comnum   
    denMat = zeros(nbas*udim, nbas*udim);
    numVec = zeros(nbas*udim,1);
    
    for i=1:m
        tdata = alltraindata{i}(:,1);
        target = alltraindata{i}(:,end);

        PHI = bsplineM(tdata, knots, nord);
        kUP = kron(Um(:,i),PHI');

        denMat = denMat + R(i,k)*kron( kUP*invCmk{i,k}*PHI, Um(:,i)' );
        numVec = numVec + R(i,k)*kUP*invCmk{i,k}*target;             
    end
    
    %indvec = sum(denMat ~= zeros(nbas*udim,nbas*udim));
    indvec = abs(numVec)/max(abs(numVec)) > 1e-4;
    induse = find(full(indvec),1,'last');
    if induse == nbas*udim
        Bvec = inv(denMat)*numVec;
        Bmin{k} = reshape(Bvec, nbas, udim);
    else
        denMat1 = denMat(1:induse,1:induse);
        numVec1 = numVec(1:induse);
        Bvec1 = inv(denMat1)*numVec1;
        tempBmin = reshape(Bvec1, induse/udim, udim);
        Bmin{k} = [tempBmin; zeros(nbas-induse/udim,udim)];
    end      
end


