function [B, xp, iuu] = gpalloctrain(traindata, Um, Vm, mixnum, gppno)
%
% GP Allocation training 
%
global initvalue;
global initvalue1;
global nbas;

cnum = size(traindata,1);
inputnum = size(traindata{1},2) - 2; 


tmin = 0; tmax = 0; 
for j=1:cnum
    tdata = traindata{j}(:,1);
    tmin0 = min(tdata);
    tmax0 = max(tdata);
    if tmin0 < tmin, tmin = tmin0; end
    if tmax0 > tmax, tmax = tmax0; end
end        

% Use B-spline with order 4 and basis number 20
nord = 4;
%nbas = 20;
nkno = nbas + 2 - nord;
knots = linspace(tmin, tmax, nkno);

% Initialise the B-spline coefficients
B = initvalue;
xp = initvalue1;


t = 0;
para = [mixnum, gppno, nord, knots];
while 1
    t = t + 1
	disp('Minimising the minus loglikelihood of allocation mixture of GPFR...')       
    [xp1 fxp ip] = minimize(xp, 'repgpalloclik', 10, traindata, B, Vm, Um, para);
	B1 = min_in_b(xp1, traindata, B, Vm, Um, para);   
    
    dB = 0;
    for j=1:mixnum
        dB = dB + norm(B1{j}-B{j},'fro');
    end
    dB
    dx = (xp1-xp)'*(xp1-xp)
    
    B = B1; xp = xp1;
	
    if ((dx < 1e-6)&(dB < 1e-6)) | (t>=30)
    %if ((dx < 1e-3)&(dB < 1e-3)) | (t>=30)    
        break;
    end
end

tval = tmin:0.1:tmax;
PHI = bsplineM(tval, knots, nord);
% figure;
% co = ['b'; 'g'; 'r'; 'c'; 'm'; 'y'; 'k'; 'k'];
% for j=1:mixnum
%     meanfun = PHI*B{j}*Um;
%     %meanfun = PHI*B{j};
%     plot(tval, meanfun, co(j));
%     hold on;
% end


% For computing the confidence interval
iuu = inv(Um*Um');

