function [fX, dfX] = gp01lik(X, input, target);

%GP1 Compute minus log likelihood and its derivatives with respect to
%    hyperparameters for a Gaussian Process for regression.
%
%    GP1 is an implementation of Gaussian Process for Regression, using a
%    simlpe covariance function containing a non-linear Gaussian term, a linear
%    term and a noise term. These terms are controlled by hyperparameters.
%
%    USAGE: [fX dfX] = GP1(X, input, target)
%
%    where:
%
%    X      is a (column) vector (of size 2*D+2) of hyperparameters
%    input  is a n by D matrix of training inputs
%    target is a (column) vector (of size n) of targets
%    fX     is the returned value of minus log likelihood
%    dfX    is a (column) vector (of size 2*D+2) of partial derivatives
%           of minus the log likelihood wrt each of the hyperparameters
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
%    Note: the reason why the log og the parameters are used in X is that we
%    may then used uncontrained optimisation, while the hyperparameters
%    themselves must, naturally, always be positive.
%
%    This function can conveniently be used with the CONJGRAD function to train
%    a Gaussian process:
%
%    [M, FM, IT, NB] = CONJGRAD(X, 'gp1', TOL, input, target)
%
%    See also: GP1PRED
%              CONJGRAD
%      
%    (C) Copyright 1999, Carl Edward Rasmussen (15-01-1999).

[n, D] = size(input);         % number of examples and dimension of input space
expX = exp(X);              % exponentiate the hyperparameters once and for all


% first, we write out the covariance matrix Q

Z = zeros(n,n); Q = zeros(n,n);                         % create and zero space
for d = 1:D                                           % non-linear contribution
  Z = Z + (input(:,d)*ones(1,n)-ones(n,1)*input(:,d)').^2*expX(d);
end
Z = expX(2*D+1)*exp(-0.5*Z);
for d = 1:D                                               % linear contribution
  Q = Q + expX(D+d)*input(:,d)*input(:,d)';
end
Q = Q + Z + expX(2*D+2)*eye(n);        % linear term, non-linear term and noise


% then, we compute the negative log likelihood ...

invQ = inv(Q);
logdetQ = 2*sum(log(diag(chol(Q))));            % don't compute det(Q) directly
fX = 0.5*logdetQ + 0.5*target'*invQ*target + 0.5*n*log(2*pi);


% ... and its partial derivatives

dfX = zeros(2*D+2,1);                   % set the size of the derivative vector

invQt = invQ*target;
for d = 1:D
  V = -0.5*expX(d)*(input(:,d)*ones(1,n)-ones(n,1)*input(:,d)').^2.*Z;
  dfX(d) = 0.5*sum(sum(invQ.*V)) - 0.5*invQt'*V*invQt;
end 
for d = 1:D
  V = expX(D+d)*input(:,d)*input(:,d)';
  dfX(D+d) = 0.5*sum(sum(invQ.*V)) - 0.5*invQt'*V*invQt;
end
dfX(2*D+1) = 0.5*sum(sum(invQ.*Z)) - 0.5*invQt'*Z*invQt;
dfX(2*D+2) = 0.5*trace(invQ)*expX(2*D+2) - 0.5*invQt'*invQt*expX(2*D+2);
