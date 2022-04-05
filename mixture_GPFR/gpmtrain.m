function xp = gpmtrain(traindata, Um, Vm, mixnum, gppno)
%
% GP mixture training 
%

global initvalue1;

cnum = size(traindata,1);
inputnum = size(traindata{1},2) - 2; 


% Set B-spline basis and set the mean coefficient as 0, equivalent to using GP only
tmin = 0; tmax = 0; 
for j=1:cnum
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
B = cell(mixnum,1);
for j=1:mixnum
    B{j} = zeros(nbas,1);    
end

% Initialise Gaussian process parameters
infonum = size(Vm,2);                       %the number of infomations of batches
%xp = [unifrnd(-10,0,(inputnum*gppno+2)*mixnum,1); ones(infonum*(mixnum-1),1)];
xp = initvalue1;

%Mminimise the negative log-likelihood
para = [mixnum, gppno, nord, knots];
disp('Minimising the negative loglikelihood of mixture of GP...')       
[xp1 fxp ip] = minimize(xp, 'repgpalloclik', 100, traindata, B, Vm, Um, para);
xp = xp1;


