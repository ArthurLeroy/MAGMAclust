function [traindata, rann] = sampling_3c(Vm, Um, para0)

n = 51;                % number of points
L = -4; H = 4;         % limits of points are "L" and "H"
t = (L:(H-L)/(n-1):H)'; % sample points: n x 1

sn = size(Vm,1);
y = zeros(n,sn);

beta = [para0(13); para0(14); 0];
pimk = zeros(1,3);
traindata = cell(sn,1);
rann = zeros(sn,1);

for i=1:sn
    for j = 1:3
        betak = beta(j);
        pimk(1,j) = exp( Vm(i,:)*betak );
    end
    pimk = pimk(1,:) / sum(pimk(1,:));
    
    ran = multinomial(pimk, [1 2 3], 1, 1);
        
    
    if ran == 1
        %y0 = Um(i)*0.5*sin((0.5*t).^3);  
        y0 = - Um(i)*0.8*cos(0.8*t + 4.5) - 0.2;  
        w1 = para0(1); a1 = para0(2); v1 = para0(3); v0 = para0(4);
    elseif ran == 2
        y0 = Um(i)*exp(t/5) - 1.5;      
        w1 = para0(9); a1 = para0(10); v1 = para0(11); v0 = para0(12);
    else
        y0 = Um(i)*0.8*atan(t);
        w1 = para0(5); a1 = para0(6); v1 = para0(7); v0 = para0(8);
    end    

    %x = t;
    x = samplex2(t);

    Q = v1*exp(-0.5*w1*(x*ones(1,n)-ones(n,1)*x').^2)...
        + a1*(x*(ones(1,n)).*(ones(n,1)*x')) + v0*eye(n);          
    
    y = ones(1,1)*y0 + (randn(1,n)*chol(Q))'; 
    
    traindata{i} = [t,x,y];
    rann(i) = ran; 
end


function x=samplegp(t)

n = size(t,1);
w1 = 0.5; a1 = 1e-20; v1 = 2.0; v0 = 9e-10;
%w1 = 0.3; a1 = 0.3; v1 = 1.0; v0 = 9e-10;
Q = v1*exp(-0.5*w1*(t*ones(1,n)-ones(n,1)*t').^2)...
    + a1*(t*(ones(1,n)).*(ones(n,1)*t')) + v0*eye(n);     
x = 4*exp(t/5)-5 + (randn(1,n)*chol(Q))'; 


function x = samplex1(t)

n = size(t,1);

t0 = unifrnd(0.1,1,1,1);
x = normpdf((t-max(t))/2.5,0,t0);
x = x*8/max(x) - 4;


function x = samplex2(t)

n = size(t,1);

t0 = unifrnd(0.5,1,1,1);
t1 = unifrnd(-0.5,0.5,1,1);
x = 3.7*atan(t/(t0*5)+t1) - 0.3;

