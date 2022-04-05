function [ypred, s2] = gpfrpred_new(traindata, U, testdata, UT, yregfd, xpgf, izzgf, gppno)


m = size(traindata,1);    %number of all training curves
gpparasize = (size(traindata{1},2)-2)*gppno + 2;    %w, (a), v1, v0

gppredfun = ['gpm0', num2str(gppno), 'pred(Xk, input, target, muin, xte, mute)'];
tte = testdata{1}(:,1);        %time points 
xte = testdata{1}(:,2:end-1);  %input points 
pte = length(tte);
allmu = zeros(pte,m);
alls = zeros(pte,m);
Xk = xpgf;
for i=1:m
    tdata = traindata{i}(:,1);
    input = traindata{i}(:,2:end-1);
    target = traindata{i}(:,end);
       
    muin = eval_fd(tdata,yregfd)*U(:,i); 
    mute = eval_fd(tte,yregfd)*UT;  
    [mu, s2] = eval(gppredfun);

	allmu(:,i) = mu;
	alls(:,i) = s2*(1+UT'*izzgf*UT); 
end


%Average over all predicting curves

ypred = sum( allmu, 2 ) / m ;
s2 = sum( alls, 2 )/m + sum( allmu.^2, 2 )/m - ypred.^2;



