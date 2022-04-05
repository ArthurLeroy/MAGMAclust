function y = rmse(t,a)
%compute the root mean squar error between two vectors

y = sqrt(sum((a-t).^2)/length(t));
