function [mu, sigma2] = gpm02pred(X, input, target, muin, test, mute);

%GP1PRED Compute predictions based on hyperparameters, training inputs and 
%    targets and test inputs.
%
%    USAGE: [mu sigma2] = GP1PRED(X, input, target, test) 
%
%    where: 
%
%    X      is a (column) vector (of size 2*D+2) of hyperparameters
%    input  is a n by D matrix of training inputs
%    target is a (column) vector (of size n) of targets
%    test   is a nn by D matrix of test inputs
%    mu     is a (column) vector (of size nn) of prediced means
%    sigma2 is a (column) vector (of size nn) of predicted variances
%
%    The form of the covariance function is
%
%    C(x^p,x^q) = sum_{d=1..D} a_d * x^p_d * x^q_d
%               + v1 * exp( -0.5 * sum_{d=1..D} w_d * (x^p_d - x^q_d)^2 )
%               + v0 * delta_{p,q}
%
%    where the first term is linear, the second non-linear and the final term
%    with the kronecker delta is the noise contribution. In this function, the
%    hyperparameters w_i, a_i, v1 and v0 are collected in the vector X as
%    follows:
%
%    X = [ log(w_1)
%           .
%          log(w_D)
%          log(a_1)
%           .
%          log(a_D)
%          log(v1)
%          log(v0) ]
%
%    which corresponds exactly to the covariance function used by GP1.
%
%    See also: GP1
%              CONJGRAD
%      
%    (C) Copyright 1999, Carl Edward Rasmussen (15-01-1999).

[n, D] = size(input);   % number of training cases and dimension of input space
[nn, D] = size(test);       % number of test cases and dimension of input space
expX = exp(X);              % exponentiate the hyperparameters once and for all


% first, we write out the covariance matrix Q for the training inputs ...

[Q, Z] = covfun02(X, input);


% ... then we compute the covariance between training and test inputs ...

a = zeros(n, nn);                                       % create and zero space
for d = 1:D                                           % non-linear contribution
  a = a + (input(:,d)*ones(1,nn)-ones(n,1)*test(:,d)').^2*expX(d);
end
a = expX(2*D+1)*exp(-0.5*a);
for d = 1:D                                               % linear contribution
  a = a + expX(D+d)*input(:,d)*test(:,d)';
end


% ... and covariance between the test input and themselves 

b = expX(2*D+1) + test.^2*expX(D+1:2*D);


% Now, write back mean prediction and variance

invQ = inv(Q);
mu = mute + a'*invQ*(target - muin);
sigma2 = b - sum(a.*(invQ*a))';

