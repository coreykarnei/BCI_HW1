function res = rmbaseline(sig)

flag = 0;
if size(sig,1) < size(sig,2)
    sig = sig';
    flag = 1;
end

for i = 1:size(sig,2)
    sig(:,i) = sig(:,i) - mean(sig(:,i));
end

if flag == 1
    sig = sig';
end

res = sig;