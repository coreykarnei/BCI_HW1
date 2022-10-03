function [ar,e] = mmvar(data,order)
% mvar for multi-trials
% data: multi-trials data, channel number x data length x trial number
% trial_num: trial number
% ar: auto regression coefficients
% e: covariance matrix of noise

if nargin < 1
    error('mmvar has no input !!!');
    return;
end

if nargin == 1
    order = 0;
end

if length(size(data)) == 3
    trial_num = size(data,3);
else
    trial_num = 1;
end

p = size(data,1);

R = zeros(p,p,order+1);
for i = 1:trial_num
    % autovariance matrix
    trial = squeeze(data(:,:,i));
    for j = 0:order
        R(:,:,j+1) = R(:,:,j+1) + avm(trial,j);
    end
end
% average of autovariance matrix across trials
R = R./trial_num;

% matrix Y and X from R
% Y = X * A, where size Y = mp x p, size X = mp x mp
X = [];
Y = [];
for i = 1:order
    Xtmp = [];
    for j = 1:order
        Xtmp = [Xtmp R(:,:,abs(j-i)+1)];
    end
    X = [X;Xtmp];
    Y = [Y;R(:,:,i+1)];
end

% auto regression coefficient
A = pinv(X)*Y;
for i = 1:order
    ar(:,:,i) = A((i-1)*p+1:i*p,:);
end

% noise covariance matrix
e = squeeze(R(:,:,1));
for i = 1:order
    e = e + squeeze(ar(:,:,i)) * squeeze(R(:,:,i+1));
end


function Rn = avm(X,n)
% autovariance matrix for p channel N points data
% X: data with p channel and m points
% n: delay of autovariance matrix
p = size(X,1);
N = size(X,2);

Rn = zeros(p,p);
for i = 1:N-n
    Rn = Rn + 1/(N-n) * X(:,i) * X(:,i+n)';
end