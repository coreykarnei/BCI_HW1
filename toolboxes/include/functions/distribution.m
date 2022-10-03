function [w value] = distribution(x, num, xmin, xmax)

if nargin < 2 | nargin == 3
    return
elseif nargin == 2
    xmin = min(x);
    xmax = max(x);
end


value = xmin : (xmax - xmin) / num : xmax;
step = (xmax - xmin) / num;
value = value + (step)/2;
value(end) = '';

w = zeros(1,num);

for i = 1:length(x)
    if x(i) ~= xmax
        index = floor((x(i)-xmin)/step) + 1;
        w(index) = w(index) + 1;
    end
end