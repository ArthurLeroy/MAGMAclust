function [output] = multinomial(ps,outcomes,nr,nc)
%MULTIN Create random numbers from a multinomial distribution.
%
% USE:
% [output] = multinomial(ps, ks, nr, nc)
%
% ps = probability of outcomes in a vector,
% ks = the assigned values for each probabilities,
% nr = the number of rows for the output, and
% nc = the number of columns for the output.

% Tomo Eguchi
% 27 January 2000

if (nargin == 1),
   error('You need at least two input arguments.'),
elseif (nargin == 2),
   nr = 1;
   nc = 1;
elseif (nargin == 3),
   nc = 1;
elseif (nargin > 4),
   error('Incorrect number of input arguments.  Require four or less.'), 
end

[nr1, nc1] = size(ps);
if (nr1 > 1),
   error('Probabilities have to be in a vector.'); end

errorps = 1 - sum(ps);
if (abs(errorps) > eps),
   error('Probabilities have to add up to one.'), end

[nr2, nc2] = size(outcomes);
if (nr1 ~= nr2 | nc1 ~= nc2),
   error('Dimensions of input vectors did not match.'), end

output = zeros(nr, nc);    % define the output matrix.
randomn = rand(nr, nc);   % create a random number matrix.
endpoints = [0,cumsum(ps)];   % find the endpoints for each category in prob.

for c1 = 1:nc1,
   maxp = endpoints(c1+1);
   minp = endpoints(c1);
   index = find(randomn < maxp & randomn >= minp);
   output(index) = outcomes(c1);
end

