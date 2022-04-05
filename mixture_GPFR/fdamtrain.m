function [B, beta, sigfd] = fdamtrain(traindata, Um, Vm, mixnum, gppno)
%
% FDA mixture training 
%

global initvalue;
global nbas;

cnum = size(traindata,1);
inputnum = size(traindata{1},2) - 2; 

tmin = 0; tmax = 0; tlengmax = 0; 
for j=1:cnum
    tdata = traindata{j}(:,1);
    tmin0 = min(tdata);
    tmax0 = max(tdata);
    if tmin0 < tmin, tmin = tmin0; end
    if tmax0 > tmax, tmax = tmax0; end
    if tlengmax < size(tdata,1), tlengmax = size(tdata,1); end
end     

% Use B-spline with order 4 and basis number 20
nord = 4;
%nbas = 20;
nkno = nbas + 2 - nord;
knots = linspace(tmin, tmax, nkno);

% Initialise the B-spline coefficients
B = initvalue;
inittheta = zeros(inputnum*gppno+2,1);
inittheta(end) = 1;             %set the parameters of covariance being 0 except noise
infonum = size(Vm,2);           %the number of infomations of batches
initbeta = ones(infonum,1);
beta = repmat(initbeta,mixnum-1,1);


t = 0;
para = [mixnum, gppno, nord, knots];
while 1
    t = t + 1
    
	disp('Minimising the minus loglikelihood of mixture of FDA...')       
    [beta1 fxp ip] = minimize(beta, 'repfdamlik', 10, traindata, B, Vm, Um, para);
    xp = [repmat(inittheta,mixnum,1); beta1];
	B1 = min_in_b(xp, traindata, B, Vm, Um, para);   
    
    dB = 0;
    for j=1:mixnum
        dB = dB + norm(B1{j}-B{j},'fro');
    end
    dB
    dbeta = (beta1-beta)'*(beta1-beta)
    
    B = B1; beta = beta1;
	
    if ((dbeta < 1e-6)&(dB < 1e-6)) | (t>=30)
        break;
    end
end

tval = tmin:0.1:tmax;
PHI = bsplineM(tval, knots, nord);
% figure;
% for j=1:mixnum
%     meanfun = PHI*B{j}*Um;
%     plot(tval, meanfun, 'r');
%     %pause;
%     hold on;
% end


% For computing the confidence interval
beta1 = [beta;zeros(infonum,1)];
pimk = zeros(cnum, mixnum);
for i = 1:cnum
    for j = 1:mixnum
        betak = beta1( (j-1)*infonum+1 : j*infonum );
        pimk(i,j) = exp( Vm(i,:)*betak );
    end
    pimk(i,:) = pimk(i,:) / sum(pimk(i,:)) ;
end
allmuhat2 = cell(cnum,1);
for i=1:cnum
    tdata = traindata{i}(:,1);
    target = traindata{i}(:,end);
    PHI = bsplineM(tdata, knots, nord);
    mixmu = zeros(length(tdata),mixnum);
    for j=1:mixnum, mixmu(:,j) = PHI*B{j}*Um(:,i); end
    allmu = mixmu*pimk(i,:)';
    allmuhat2{i} = [tdata, (allmu - target).^2];    
end
sigfd = comp_mean_curve(allmuhat2);
tgrid = (linspace(tmin, tmax, tlengmax))';
sig = eval_fd(tgrid,sigfd);
sig = sqrt(sig.^2);
bsbasis = create_bspline_basis([tmin,tmax], 20);
sigfd = data2fd(sig, tgrid, bsbasis);

