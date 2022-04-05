function [pp pik] = postpi(traindata, allinputdata, B, xp, Vm, Um, comnum, gppno)

infonum = size(Vm,2);                %the number of infomations of patients
inputnum = size(traindata{1},2) - 2; 
if length(xp) == infonum*(comnum-1)
    inittheta = zeros(inputnum*gppno+2,1);
    inittheta(end) = 1;
    xp = [repmat(inittheta,comnum,1); xp];
end


% Use B-spline with order 4 and basis number 20
nord = 4;
nbas = 20;
nkno = nbas + 2 - nord;
tmin = 0; tmax = 0; 
for j=1:size(traindata,1)
    tdata = traindata{j}(:,1);
    tmin0 = min(tdata);
    tmax0 = max(tdata);
    if tmin0 < tmin, tmin = tmin0; end
    if tmax0 > tmax, tmax = tmax0; end
end      
knots = linspace(tmin, tmax, nkno);


m = size(allinputdata,1);  %number of all test curves

fXt = zeros(m,comnum);  %
dfXt = cell(m,comnum);  %
gplikfun = ['gpm0', num2str(gppno), 'lik(Xk, input, target, mu)'];

for i = 1:m    
    tdata = allinputdata{i}(:,1);
    input = allinputdata{i}(:,2:end-1);
    target = allinputdata{i}(:,end);
    tdsize = size(input,1);
    gpparasize = size(input,2)*gppno + 2;    %w, (a), v1, v0
    
    PHI = bsplineM(tdata, knots, nord);
    
    for k=1:comnum
        Xk = xp((k-1)*gpparasize+1 : k*gpparasize);        
        mu = PHI*B{k}*Um(:,i);
        [a, da] = eval(gplikfun);
        fXt(i,k) = a;
        dfXt{i,k} = da;        
    end
end


pik = zeros(1,comnum);
beta = [xp( comnum*gpparasize+1 : end ); zeros(infonum,1)];
for j = 1:comnum
    betak = beta( (j-1)*infonum+1 : j*infonum );
    pik(1,j) = exp( Vm(1,:)*betak );
end
pik(1,:) = pik(1,:) / sum(pik(1,:)) ;

    
fXt = fXt + repmat(max(-fXt,[],2),1,comnum);      %to avoid 'Not a Number'
fX = sum(fXt,1);
pp = exp(-fX(1,:)).* pik(1,:)/(exp(-fX(1,:))* pik(1,:)');


