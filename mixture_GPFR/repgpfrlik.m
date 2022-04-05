function [fX, dfX] = repgpfrlik(X, alltraindata, yregfd, gppno)
%
%Used by GPFR
%

m = size(alltraindata,1);  %number of all training curves
gplikfun = ['gpm0', num2str(gppno), 'lik(X, input, target, mu)'];

fX = 0;
dfX = 0;
for i = 1:m    
    tdata = alltraindata{i}(:,1); 
    input = alltraindata{i}(:,2:end-1);
    target = alltraindata{i}(:,end);
      
    mu = eval_fd(tdata,yregfd);
    [fXt, dfXt] = eval(gplikfun);
    fX = fX + fXt;
    dfX = dfX + dfXt;
end