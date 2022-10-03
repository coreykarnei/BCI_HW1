function [CDSP CDSP_feature DP] = cva(class1,class2)

if(nargin < 2)
    error('cva has no input !!!');
    return;
end

class{1} = class1;
class{2} = class2;

% initialization of between-classes dispersion matrix (B) and
% within-classes dispersion matrix (W)
B = zeros(size(class{1},2),size(class{1},2));
W = zeros(size(class{1},2),size(class{1},2));


% compute average for each class and overall dataset
avg_all = zeros(1,size(class{1},2));
num_all = 0;
for i = 1:length(class) 
    avg_class{i} = mean(class{i},1);
    avg_all = avg_all + avg_class{i}*size(class{i},1);
    num_all = num_all + size(class{i},1);
end
avg_all = avg_all / num_all;

% update matrix B and matrix W for each class
for i = 1:length(class)
    
   % update matrix B
   B = B + size(class{i},1)*(avg_class{i} - avg_all)'*(avg_class{i} - avg_all);
   
   % update matrix W
   for j = 1:size(class{i},1)
       W = W + (class{i}(j,:) - avg_class{i})'*(class{i}(j,:) - avg_class{i});
   end
   
end

% eigenvector and eigenvalues of inv(W)*B
[A,lambda] = eig(pinv(W)*B);
lambda = sum((abs(lambda)),1);

% get first (classnum - 1) eigenvectors and eigenvalues
A = A(:,1:length(class)-1);
lambda = lambda(1:length(class)-1);

% canonical discriminant spatial pattern (CDSP) for each class
CDSP = A;

for i = 1:length(class)
    for j = size(A,2)
        CDSP_feature{i}(j,:) = class{i}*A(:,j);
    end
end

% Discriminant Power (DP) for each electrode
DP = zeros(size(A,1),1);
for j = size(A,2)
    DP = DP + A(:,j).^2;
end
DP = DP / sum(DP);
