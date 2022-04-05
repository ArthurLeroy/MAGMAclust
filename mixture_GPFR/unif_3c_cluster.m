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
sn = sn1+sn2+sn3;         % number of curves for all classes
V = ones(sn,1);           % covariate for logistic model
V(1:sn1,1) = 2;           % for first component
V(sn1+1:sn1+sn2,1) = -1;  % for second component
V(sn1+sn2+1:sn,1) = -0.8;          
U = ones(1,sn);           % covariate for mean functions

%Set up true values
para0 = [1.5; 1e-20; 0.1; 64e-4; 1.0; 1e-20; 0.05; 64e-4; 0.5; 1e-20; 0.07; 64e-4; 2.0; 1.0];
[traindata,rann] = sampling_3c_revise(V, U, para0);

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
gppno = 2;    
udim = size(U,1);
nbas = 10;
initvalue = cell(mixnum,1);
for j=1:mixnum
    initvalue{j} = zeros(nbas,udim);    
end
initvalue1 = [log(para0(1:12)); para0(13:14)];

%Training
[B, xp, iuu] = gpalloctrain(traindata, U, V, mixnum, gppno);      % GPALLOC method
[Bfda, xpfda, iuufda] = fdamtrain(traindata, U, V, mixnum, gppno);% FDA mixutre method
xpgpm = gpmtrain(traindata, U, V, mixnum, gppno);                 % GP method


%Simulation study for prediction (deleted because not in use, see unif_3c.m)

%Clustering

testdata = cell(1,1);
gpa_cluster = zeros(sn,1);
fda_cluster = zeros(sn,1);
gp_cluster = zeros(sn,1);
for i = 1:sn
    VT = V(i);     
    UT = U(i);
    testdata{1} = traindata{i};
    
    pp = postpi(traindata, testdata, B, xp, VT, UT, mixnum, gppno);
    [ppm ppi] = max(pp); gpa_cluster(i) = ppi;   
    
    newxpfda = [0;0;0;1;0;0;0;1;0;0;0;1;xpfda];
    pp = postpi(traindata, testdata, Bfda, newxpfda, VT, UT, mixnum, gppno);
    [ppm ppi] = max(pp); fda_cluster(i) = ppi;

    newB = B; newB{1} = B{1}.*0; newB{2} = B{2}.*0; newB{3} = B{3}.*0;
    pp = postpi(traindata, testdata, newB, xpgpm, VT, UT, mixnum, gppno);
    [ppm ppi] = max(pp); gp_cluster(i) = ppi;
end


%Simulation study for model selection

bv = zeros(5,1);
BB = cell(5,1);
xpB = cell(5,1);
for mixno = 1:5
    udim = size(U,1);
    initvalue = cell(mixno,1);
    for j=1:mixno
        initvalue{j} = zeros(nbas,udim);    
    end  
    initvalue1 = [unifrnd(-5,0,4*mixno,1); unifrnd(0,2,mixno-1,1)];
    if mixno == 3
        Bms = B; xpms = xp; iuums = iuu;
    else
        [Bms, xpms, iuums] = gpalloctrain(traindata, U, V, mixno, gppno);
    end
    bv(mixno) = bicvalue(Bms, xpms, traindata, U, V, mixno, gppno);
    BB{mixno} = Bms;
    xpB{mixno} = xpms;
end
figure;
plot(1:5,bv,'-o');


