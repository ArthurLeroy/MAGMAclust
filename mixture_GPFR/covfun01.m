function Q = trycovfun(X, input)
% The form of the covariance function is
%
% C(x^p,x^q) = v1 * exp( -0.5 * sum_{d=1..D} w_d * (x^p_d - x^q_d)^2 )
%            + v0 * delta_{p,q}
%
% where the first term is Gaussian and the second term with the kronecker
% delta is the noise contribution. In this function, the hyperparameters w_i,
% v1 and v0 are collected in the vector X as follows:
%
% X = [ log(w_1)
%       log(w_2) 
%        .
%       log(w_D)
%       log(v1)
%       log(v0) ]


[n, D] = size(input);         % number of examples and dimension of input space
expX = exp(X);              % exponentiate the hyperparameters once and for all


% first, we write out the covariance matrix Q

Z = zeros(n,n);
for d = 1:D
  Z = Z + (repmat(input(:,d),1,n)-repmat(input(:,d)',n,1)).^2*expX(d);
end
Z = expX(D+1)*exp(-0.5*Z);
Q = Z + expX(D+2)*eye(n);                             % Gaussian term and noise
