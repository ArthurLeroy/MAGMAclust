function [ypred, s2] = gpallocpred_new(traindata, U, V, testdata, UT, VT, B, xp, iuu, comnum, gppno)

if ~iscell(B)             %mixture of gp
    Bgp = cell(comnum,1);
    for j=1:comnum
        Bgp{j} = zeros(20,1);
    end
    B = Bgp;
end
    

m = size(traindata,1);    %number of all training curves
infonum = size(VT,2);                %the number of infomations of patients
if length(xp) == infonum*(comnum-1)
    gpparasize = 0;
    flag = 'f';     %fda mixture
else
    gpparasize = (size(traindata{1},2)-2)*gppno + 2;    %w, (a), v1, v0
    flag = 'a';    %gpalloc
end



% Use B-spline with order 4 and basis number 20
tmin = 0; tmax = 0; 
for j=1:m
    tdata = traindata{j}(:,1);
    tmin0 = min(tdata);
    tmax0 = max(tdata);
    if tmin0 < tmin, tmin = tmin0; end
    if tmax0 > tmax, tmax = tmax0; end
end        
nord = 4;
nbas = 20;
nkno = nbas + 2 - nord;
knots = linspace(tmin, tmax, nkno);


gppredfun = ['gpm0', num2str(gppno), 'pred(Xk, input, target, muin, xte, mute)'];
tte = testdata{1}(:,1);        %time points 
xte = testdata{1}(:,2:end-1);  %input points 
pte = length(tte);
allmu = zeros(pte,m);
alls = zeros(pte,m);
for i=1:m
    tdata = traindata{i}(:,1);
    input = traindata{i}(:,2:end-1);
    target = traindata{i}(:,end);
    
    PHIin = bsplineM(tdata, knots, nord);
    PHIte = bsplineM(tte, knots, nord);    
    
    mixmu = zeros(pte,comnum);
    mixs = zeros(pte,comnum);
    for k=1:comnum
        if flag == 'a'
            Xk = xp((k-1)*gpparasize+1 : k*gpparasize);        
            muin = PHIin*B{k}*U(:,i);   
            mute = PHIte*B{k}*UT;  
            [mu, s2] = eval(gppredfun);
            mixmu(:,k) = mu;
            mixs(:,k) = s2*(1+UT'*iuu*UT);
        end
        if flag == 'f'
            mixmu(:,k) = PHIte*B{k}*UT;
            sig2 = eval_fd(tte,iuu);
            mixs(:,k) = sig2*(1+UT'*inv(U*U')*UT);    
        end
    end

    %Compute the posterior mean of mixing weights
    inputdata = cell(1,1); inputdata{1} = traindata{i};
    [ppikTe pikTe] = postpi(traindata, inputdata, B, xp, VT, UT, comnum, gppno);

    %Overall prediction for m-th curve
    allmu(:,i) = mixmu*pikTe';
	alls(:,i) = mixs*pikTe' + mixmu.^2*pikTe' - allmu(:,i).^2; 
end


%Weighted average over all predicting curves

pimkTr = zeros(m, comnum);
pikTe = zeros(1,comnum);
beta = [xp( comnum*gpparasize+1 : end ); zeros(infonum,1)];

for i = 1:m
    for j = 1:comnum
        betak = beta( (j-1)*infonum+1 : j*infonum );
        pimkTr(i,j) = exp( V(i,:)*betak );
    end
    pimkTr(i,:) = pimkTr(i,:) / sum(pimkTr(i,:)) ;
end

for j = 1:comnum
    betak = beta( (j-1)*infonum+1 : j*infonum );
    pikTe(1,j) = exp( VT(1,:)*betak );
end
pikTe(1,:) = pikTe(1,:) / sum(pikTe(1,:)) ;


kl = zeros(m,1);
for i = 1:m
    kl(i) = pikTe*(log(pikTe./pimkTr(i,:)))';
end
if sum(kl == zeros(m,1)) > 0
    index = 1:m;
    induse = index(kl==0);
    curveno = size(induse,2);
    w = ones(curveno,1)/curveno;    
else
    induse = 1:m;
    curveno = size(induse,2);
    wun = 1./kl;
    w = wun/sum(wun);
end

ypred = allmu(:,induse) * w;
s2 = alls(:,induse)*w + allmu(:,induse).^2*w - ypred.^2;



