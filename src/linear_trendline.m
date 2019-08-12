function [m,b,Rsq] = linear_trendline(x,y)
X = [ones(length(x),1) x];
a = X\y;
b = a(1);
m = a(2);
if nargout > 2
    Rsq = 1 - sum((y - (m*x + b)).^2)/sum((y - mean(y)).^2);
end
end

