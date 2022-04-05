function [ypred, s2] = gpallocpred(traindata, inputdata, U, V, testdata, UT, VT, B, xp, iuu, comnum, gppno)

if ~iscell(B)             %mixture of gp
    Bgp = cell(comnum,1);
    for j=1:comnum
        Bgp{j} = zeros(20,1);
    end
    B = Bgp;
end
    

infonum = size(VT,2);                %the number of infomations of patients
if length(xp) == infonum*(comnum-1)
    gpparasize = 0;
    flag = 'f';     %fda mixture
else
    gpparasize = (size(inputdata{1},2)-2)*gppno + 2;    %w, (a), v1, v0
    flag = 'a';    %gpalloc
end

%Compute the posterior mean of mixing weights
pikTe = postpi(traindata, inputdata, B, xp, VT, UT, comnum, gppno);


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
gppredfun = ['gpm0', num2str(gppno), 'pred(Xk, input, target, muin, xte, mute)'];
tdata = inputdata{1}(:,1);
input = inputdata{1}(:,2:end-1);
target = inputdata{1}(:,end);
tte = testdata{1}(:,1);        %time points 
xte = testdata{1}(:,2:end-1);  %input points 
pte = length(tte);
PHIin = bsplineM(tdata, knots, nord);
PHIte = bsplineM(tte, knots, nord);  

mixmu = zeros(pte,comnum);
mixs = zeros(pte,comnum);
for k=1:comnum
    if flag == 'a'
        Xk = xp((k-1)*gpparasize+1 : k*gpparasize);        
        muin = PHIin*B{k}*UT;   
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
ypred = mixmu*pikTe';
s2 = mixs*pikTe' + mixmu.^2*pikTe' - ypred.^2; 

