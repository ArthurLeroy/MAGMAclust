%
%unif_3c.m: 
%

addpath('..\..\..\..\help\fda');

global initvalue;
global initvalue1;
global nbas;

%Sampling training data
sn1 = 30;                 % number of curves for 1st classes
sn2 = 30;                 % number of curves for 2nd classes
sn3 = 30;
sn = sn1+sn2+sn3;             % number of curves for all classes
V = ones(sn,1);          % covariate for logistic model
V(1:sn1,1) = 1;          % for first component
V(sn1+1:sn1+sn2,1) = -1;  % for second component
V(sn1+sn2+1:sn,1) = 0.01;      
U = ones(1,sn);          % covariate for mean functions

%Set up true values
para0 = [1.5; 1e-20; 0.1; 25e-4; 1.0; 1e-20; 0.05; 25e-4; 0.5; 1e-20; 0.07; 25e-4; 2.0; 1.0];
[traindata,rann] = sampling_3c(V, U, para0);

%Plot the sample curves
n = size(traindata{1},1);
ydata = zeros(n,sn);
xdata = zeros(n,sn);
tdata = traindata{1}(:,1);
for i=1:sn
    ydata(:,i) = traindata{i}(:,3);
    xdata(:,i) = traindata{i}(:,2);
end
%plot(tdata,ydata(:,1:sn1),'Color',[0.5 0.5 0.5]);
plot(tdata,ydata(:,1:sn1),'g');
hold on;
plot(tdata,ydata(:,sn1+1:sn1+sn2),'b');
hold on;
plot(tdata,ydata(:,sn1+sn2+1:sn),'r');
figure;
plot(tdata,xdata(:,1:sn1),'r', tdata,xdata(:,sn1+1:sn1+sn2),'g', tdata,xdata(:,sn1+sn2+1:sn),'b');
%Plot differnet components in different gray lavels
ind1 = find(rann==1); ind2 = find(rann==2); ind3 = find(rann==3);
plot(tdata,ydata(:,ind1),'Color',[0.5 0.5 0.5]);
hold on;
plot(tdata,ydata(:,ind2),'k');
hold on;
plot(tdata,ydata(:,ind3),'Color',[0.8 0.8 0.8]);
%Plot differnet components in different colors
ind1 = find(rann==1); ind2 = find(rann==2); ind3 = find(rann==3);
plot(tdata,ydata(:,ind1),'r');
hold on;
plot(tdata,ydata(:,ind2),'b');
hold on;
plot(tdata,ydata(:,ind3),'g');

%Training

%Set up some initial values for training
mixnum = 3;   %number of mixture components
gppno = 2;    %type of GP covariance: 
              % 1 --> gpm01lik & gpm01pred & covfun01 : no linear part
              % 2 --> gpm02lik & gpm02pred & covfun02 : with linear part
udim = size(U,1);
nbas = 20;
initvalue = cell(mixnum,1);
for j=1:mixnum
    initvalue{j} = zeros(nbas,udim);    
end
initvalue1 = [log(para0(1:12)); para0(13:14)];

%Training
[B, xp, iuu] = gpalloctrain(traindata, U, V, mixnum, gppno);      % GPALLOC method
[yregfd, xpgf, izzgf] = gpfrtrain(traindata, U, gppno);           % GPFR method
[Bfda, xpfda, iuufda] = fdamtrain(traindata, U, V, mixnum, gppno);% FDA mixutre method
xpgpm = gpmtrain(traindata, U, V, mixnum, gppno);                 % GP method


%Simulation study

%Prediction a completely new curve
simtimes = 1;                %numbers of simulation running
simrcc = cell(simtimes,1);
savefig = 1;                  %not save the figures of results
i = 1;
rcc = zeros(4, simtimes, 2);
while i <= simtimes
	VT = 1;     
    UT = 1;
    testdata = sampling_3c(VT, UT, para0);
    
    %predict by gpalloc
    [predc, preds2] = gpallocpred_new(traindata, U, V, testdata, UT, VT, B, xp, iuu, mixnum, gppno);
    rcc(1,i,:) = assess(predc,preds2,testdata, savefig);

    %predict by gpfr
    [predc, preds2] = gpfrpred_new(traindata, U, testdata, UT, yregfd, xpgf, izzgf, gppno);
    rcc(2,i,:) = assess(predc,preds2,testdata, savefig);

    %predict by fda mixture
    [predc, preds2] = gpallocpred_new(traindata, U, V, testdata, UT, VT, Bfda, xpfda, iuufda, mixnum, gppno);
    rcc(3,i,:) = assess(predc,preds2,testdata, savefig);
    
    %predict by gp mixture
    [predc, preds2] = gpallocpred_new(traindata, U, V, testdata, UT, VT, 0, xpgpm, iuu, mixnum, gppno);
    rcc(4,i,:) = assess(predc,preds2,testdata, savefig);
       
    i = i + 1;
end
avrcc1 = sum(rcc,2)/simtimes


%Prediction randomly selected half curve
simtimes = 1;                %numbers of simulation running
simrcc = cell(simtimes,1);
savefig = 1;                  %not save the figures of results
i = 1;
rcc = zeros(4, simtimes, 2);
while i <= simtimes
	VT = 1;     
    UT = 1;
    testdata = sampling_3c(VT, UT, para0);
    lt = length(testdata{1}(:,1));
    [indinput,indtest] = srswor(lt,fix(lt/2));
    tempmat1 = testdata{1}(indinput,:);
    tempmat2 = testdata{1}(indtest,:);
    inputdata = cell(1,1);
    inputdata{1} = tempmat1;
    testdata{1} = tempmat2;
    
    %predict by gpalloc
    [predc, preds2] = gpallocpred(traindata, inputdata, U, V, testdata, UT, VT, B, xp, iuu, mixnum, gppno);
    rcc(1,i,:) = assess(predc,preds2,testdata, savefig);

    %predict by gpfr
    [predc, preds2] = gpfrpred(inputdata, U, testdata, UT, yregfd, xpgf, izzgf, gppno);
    rcc(2,i,:) = assess(predc,preds2,testdata, savefig);

    %predict by fda mixture   
    [predc, preds2] = gpallocpred(traindata, inputdata, U, V, testdata, UT, VT, Bfda, xpfda, iuufda, mixnum, gppno);
    rcc(3,i,:) = assess(predc,preds2,testdata, savefig);
    
    %predict by gp mixture
    [predc, preds2] = gpallocpred(traindata, inputdata, U, V, testdata, UT, VT, 0, xpgpm, iuu, mixnum, gppno);
    rcc(4,i,:) = assess(predc,preds2,testdata, savefig);
       
    i = i + 1;
end
avrcc2 = sum(rcc,2)/simtimes


%Prediction the second half curve
simtimes = 1;                %numbers of simulation running
simrcc = cell(simtimes,1);
savefig = 1;                  %not save the figures of results
i = 1;
rcc = zeros(4, simtimes, 2);
while i <= simtimes
	VT = 1;     
    UT = 1;
    testdata = sampling_3c(VT, UT, para0);
    lt = length(testdata{1}(:,1));
    indinput = 1:fix(lt/2); indtest = fix(lt/2)+1:lt;
    tempmat1 = testdata{1}(indinput,:);
    tempmat2 = testdata{1}(indtest,:);
    inputdata = cell(1,1);
    inputdata{1} = tempmat1;
    testdata{1} = tempmat2;
    
    %predict by gpalloc
    [predc, preds2] = gpallocpred(traindata, inputdata, U, V, testdata, UT, VT, B, xp, iuu, mixnum, gppno);
    rcc(1,i,:) = assess(predc,preds2,testdata, savefig);

    %predict by gpfr
    [predc, preds2] = gpfrpred(inputdata, U, testdata, UT, yregfd, xpgf, izzgf, gppno);
    rcc(2,i,:) = assess(predc,preds2,testdata, savefig);

    %predict by fda mixture   
    [predc, preds2] = gpallocpred(traindata, inputdata, U, V, testdata, UT, VT, Bfda, xpfda, iuufda, mixnum, gppno);
    rcc(3,i,:) = assess(predc,preds2,testdata, savefig);
    
    %predict by gp mixture
    [predc, preds2] = gpallocpred(traindata, inputdata, U, V, testdata, UT, VT, 0, xpgpm, iuu, mixnum, gppno);
    rcc(4,i,:) = assess(predc,preds2,testdata, savefig);
       
    i = i + 1;
end
avrcc3 = sum(rcc,2)/simtimes


