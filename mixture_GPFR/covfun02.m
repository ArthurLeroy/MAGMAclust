function [Q, Z] = covfun(X, input)
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


[n, D] = size(input);         % number of examples and dimension of input space
expX = exp(X);              % exponentiate the hyperparameters once and for all


Z = zeros(n,n); Q = zeros(n,n);                         % create and zero space
for d = 1:D                                           % non-linear contribution
  Z = Z + (input(:,d)*ones(1,n)-ones(n,1)*input(:,d)').^2*expX(d);
end
Z = expX(2*D+1)*exp(-0.5*Z);
for d = 1:D                                               % linear contribution
  Q = Q + expX(D+d)*input(:,d)*input(:,d)';
end
Q = Q + Z + expX(2*D+2)*eye(n);        % linear term, non-linear term and noise