function d = discriminant(train_feature, train_d, test_feature, n, lambda, gamma)
% discriminant was used for discriminant analysis
%
% d = discriminant(training, testing, n, lq, lambda, gamma)
%
% training: training set, training{1:n}, size(training{i}) = feature number x sample number
% testing: testing set, testing(feature number x sample number)
% n: number of classes, integer
% lq: linear or quadratic discriminant analysis
% lambda, gamma: regularization parameters
% d: output, testing decision set

%%
if nargin < 4
    error('Not enough arguments!');
    return;
end

if nargin > 6
    error('Too many arguments!');
    return;
end

if nargin <= 5
    gamma = 0;
end

if nargin <= 4
    % typical LDA
    lambda = 1;
end

% compute covariance across all training data
covmx_all = cov(train_feature');

for i = 1:n
    % compute mean and covariance for each class
    training_class = train_feature(:,train_d == i);
    avg(i,:) = mean(training_class,2);
    covmx(i,:,:) = cov(training_class');
    prior(i) = size(training_class,2);
end
% get prior probability
prior = prior ./ size(train_feature,2);


% compute regularized covariance matrix
for i = 1:n
    covmx_lambda(i,:,:) = squeeze((1-lambda)*covmx(i,:,:)) + lambda*covmx_all;
    covmx_lambda_gamma(i,:,:) = squeeze((1-gamma)*covmx_lambda(i,:,:)) + ...
        gamma*trace(covmx_lambda(i))/n*eye(size(covmx_all,1),size(covmx_all,1));
end

% terms for output
for i = 1:n
    inv_conmx(i,:,:) = pinv(squeeze(covmx_lambda_gamma(i,:,:)));
    log_det_conmx(i) = log(det(squeeze(covmx_lambda_gamma(i,:,:))));
    log_prior(i) = log(prior(i));
end

% output of testing set
for j = 1:size(test_feature,2)
    for i = 1:n
        w(j,i) = (test_feature(:,j)' - avg(i,:)) * squeeze(inv_conmx(i,:,:)) * (test_feature(:,j)' - avg(i,:))'...
            + log_det_conmx(i) - 2*log_prior(i);
    end
    [min_value min_index] = min(w(j,:));
    d(j) = min_index;
end
