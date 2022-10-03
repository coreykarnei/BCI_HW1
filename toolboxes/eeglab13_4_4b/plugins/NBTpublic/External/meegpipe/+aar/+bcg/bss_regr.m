function myNode = bss_regr(varargin)
% BSS_REGR - BCG correction using simple BSS+Regression

import meegpipe.node.*;
import misc.process_arguments;
import misc.split_arguments;

opt.EventSelector = physioset.event.class_selector('Type', 'QRS');
opt.Correction    = 50; % In percentage, 100% harshest correction
[~, opt] = process_arguments(opt, varargin);

[~, varargin] = split_arguments(opt, varargin);


%% Regression filter
myRegrFilter = filter.mlag_regr('Order', 10);

%% Component selection criterion

myFeat = spt.feature.erp(...
    'EventSelector', opt.EventSelector, ...
    'Offset',        -0.5, ...
    'Duration',      1 ...
    );
    

myCrit = spt.criterion.threshold(myFeat, ...
    'Max',      1-opt.Correction/100, ...
    'MinCard',  5); 

%% Post-processing component filter
myPCA = spt.pca('MaxCard', @(lambda) ceil(floor(0.4*numel(lambda))));
tpcaFilter = filter.tpca('Order', @(sr) min(50, ceil(sr/10)), 'PCA', myPCA);
tpcaFilter = filter.sliding_window(tpcaFilter, ...
    'WindowLength', @(sr) 30*sr, ...
    'WindowOverlap', 50);

%% PCA
myPCA = spt.pca(...
    'LearningFilter',   @(sr) filter.lpfilt('fc', 20/(sr/2)), ...
    'RetainedVar',      99.75, ...
    'MaxCard',          40,   ...
    'MinCard',          @(lambda) max(3, ceil(0.05*numel(lambda))) ...
    );

%% Putting it all together: building the node
mySel1 = pset.selector.sensor_class('Class', 'EEG');
mySel2 = pset.selector.good_data;
mySel  = pset.selector.cascade(mySel1, mySel2);

myNode = bss.new(...
    'DataSelector', mySel, ...
    'BSS',          spt.bss.multicombi, ...
    'PCA',          myPCA, ...
    'Filter',       tpcaFilter, ...
    'RegrFilter',   myRegrFilter, ...
    'Criterion',    myCrit, ...
    varargin{:});

end