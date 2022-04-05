function [yregfd, xpgf, izz] = gpfrtrain(traindata, zmat, gppno);
%
%Single GPFR 
%

global initvalue1;

cnum = size(traindata,1);
inputnum = size(traindata{1},2) - 2; 

% tdata = traindata{1}(:,1);
% ydata = zeros(length(tdata),cnum);
% for j=1:cnum
%     tdata0 = traindata{j}(:,1);
%     if length(tdata0)~=length(tdata) | sum(~(tdata0==tdata))~=0 
%         error('For single GPFR method, the time of all the training data must be the same!');
%     else
%         ydata(:,j) = traindata{j}(:,end);
%     end
% end        

tmin = 0; tmax = 0; tlengmax = 0; tlengmin = 1000;
for j=1:cnum
    tdata = traindata{j}(:,1);
    tmin0 = min(tdata);
    tmax0 = max(tdata);
    if tmin0 < tmin, tmin = tmin0; end
    if tmax0 > tmax, tmax = tmax0; end
    if tlengmax < size(tdata,1), tlengmax = size(tdata,1); end
    if tlengmin > size(tdata,1), tlengmin = size(tdata,1); end
end        

tdata = zeros(tlengmax, cnum);
ydata = zeros(tlengmax, cnum);  
for i=1:cnum
    ttemp = traindata{i}(:,1);
    ytemp = traindata{i}(:,end);
    tleng = length(ytemp);
    if tleng < tlengmax
        %tdata(:,i) = [ttemp; NaN*ones(tlengmax-tleng, 1)];
        %ydata(:,i) = [ytemp; NaN*ones(tlengmax-tleng, 1)];
        tdata(:,i) = [ttemp(1:end-1); (linspace(ttemp(end), tmax, tlengmax-tleng+1))'];
        ydata(:,i) = [ytemp; ytemp(end)*ones(tlengmax-tleng, 1)];
    else
        if ttemp(end) == tmax
            tdata(:,i) = ttemp;
            ydata(:,i) = ytemp;
        else
            tmp = round(tlengmax*ttemp(end)/tmax);
            tind = srswor(tleng-1, tmp);
            tdata(:,i) = [ttemp(tind); (linspace(ttemp(end), tmax, tlengmax-tmp))'];
            ydata(:,i) = [ytemp(tind); ytemp(end)*ones(tlengmax-tmp, 1)];
        end
    end
end


%Use FDA to estimate the mean
zmat = zmat';
bsbasis = create_bspline_basis([tmin,tmax], 20);
yfd = data2fd(ydata, tdata, bsbasis);
linmodstr = linmod(zmat, yfd);
yregfd = linmodstr.reg;


%Use GP to estimate the covariance
disp('Minimising the minus loglikelihood of single GPFR...')
inittheta = unifrnd(-10,0,inputnum*gppno+2,1);  %theta contains: w, (a,) v1, v0
%inittheta = initvalue1(1:(inputnum*gppno+2));
[xpgf fxpgf ipgf] = minimize(inittheta, 'repgpfrlik', 300, traindata, yregfd, gppno);


%For computing CI
izz = inv(zmat'*zmat); 



