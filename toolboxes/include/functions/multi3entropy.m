function res = multi3entropy(A)

if(nargin < 1)
    error('No input data for multientropy!');
    return;
end

try
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            res(i,j) = entropy(A(i,j,:));
        end
    end
catch
    error('Errors in entropy computation!');
    return;
end