function [bv nllik2 pen] = bicvalue(B, xp, traindata, Um, Vm, mixnum, gppno)
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

para = [mixnum, gppno, nord, knots];
[nllik dnllik] = repgpalloclik(xp, traindata, B, Vm, Um, para);

nB = size(B);
nnB = size(B{1});
nxp = size(xp);
nD = 0;
for j=1:cnum
    nD = nD + size(traindata{j},1);
end

nllik2 = 2*nllik;
pen = (prod(nB)*prod(nnB)+prod(nxp))*log(nD);
bv = nllik2 + pen;
