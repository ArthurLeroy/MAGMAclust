function meanfd = comp_mean_curve(traindata);
%
%Compute the mean curve  
%

cnum = size(traindata,1);

tmin = 0; tmax = 0; tlengmax = 0; 
for j=1:cnum
    tdata = traindata{j}(:,1);
    tmin0 = min(tdata);
    tmax0 = max(tdata);
    if tmin0 < tmin, tmin = tmin0; end
    if tmax0 > tmax, tmax = tmax0; end
    if tlengmax < size(tdata,1), tlengmax = size(tdata,1); end
end        

tgrid = (linspace(tmin, tmax, tlengmax))';
ymean = zeros(tlengmax,1);
nymean = zeros(tlengmax,1);
for i=1:cnum
    ttemp = traindata{i}(:,1);
    ytemp = traindata{i}(:,end);    
    bsbasis = create_bspline_basis([ttemp(1),ttemp(end)], ceil((tmax-tmin)/20*(ttemp(end)-ttemp(1))) );
    yfd = data2fd(ytemp, ttemp, bsbasis);
    
    ind = find(tgrid<=ttemp(end),1,'last');
    cy = eval_fd(tgrid(1:ind),yfd);
    ymean(1:ind) = ymean(1:ind) + cy;
    nymean(1:ind) = nymean(1:ind) + 1;
end
    
ymean = ymean./nymean;
bsbasis = create_bspline_basis([tmin,tmax], 20);
meanfd = data2fd(ymean, tgrid, bsbasis);

