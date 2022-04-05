function [ypred, s2] = gpfrpred(inputdata, U, testdata, UT, yregfd, xpgf, izzgf, gppno)


gpparasize = (size(inputdata{1},2)-2)*gppno + 2;    %w, (a), v1, v0

gppredfun = ['gpm0', num2str(gppno), 'pred(Xk, input, target, muin, xte, mute)'];
tte = testdata{1}(:,1);        %time points 
xte = testdata{1}(:,2:end-1);  %input points 

Xk = xpgf;

tdata = inputdata{1}(:,1);
input = inputdata{1}(:,2:end-1);
target = inputdata{1}(:,end);

muin = eval_fd(tdata,yregfd)*UT; 
mute = eval_fd(tte,yregfd)*UT;  
[mu, s2] = eval(gppredfun);

ypred = mu;
s2 = s2*(1+UT'*izzgf*UT); 

