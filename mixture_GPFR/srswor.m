function [rn, restn] = srswor(a,c,b)
%If two inputs, randomly select c numbers from 1 to a
%rn: selected numbers
%restn: the rest a-c numbers
%If three inputs, 
%rn = [1:b, randomly select c numbers from b+1 to a]
%restn = the rest a-b-c numbers of [1, a]

if nargin == 2
    r0 = randperm(a);
    r1 = r0(1,1:c);
    r2 = r0(1,c+1:a);

    rn = sort(r1);
    restn = sort(r2);
elseif nargin == 3
    if b+c>a
        error('Selected numbers exceed the total number!');
    end

    r0 = randperm(a-b);
    r1 = r0(1,1:c);
    r2 = r0(1,c+1:a-b);

    rn = [1:b, b + sort(r1)];
    restn = b + sort(r2);
else
    error('Incorrect inputs!');
end
