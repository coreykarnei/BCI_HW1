function [popt sbc]=arorder(v, pmin, pmax, selector)
%arorder	
%
%  popt=arorder(v, pmin, pmax, selector)
%
%  where p lies between pmin and pmax and is chosen as the optimizer
%  of Schwarz's Bayesian Criterion. The input matrix v must contain
%  the time series data, with columns (first dimension) of v(:,l)
%  representing variables l=1,...,m and rows v(k,:) (second
%  dimension) of v representing observations at different (equally
%  spaced) times k=1,..,n. Optionally, v can have a third
%  dimension, in which case the matrices v(:,:, itr) represent
%  the realizations (e.g., measurement trials) itr=1,...,ntr of the
%  time series.
%
%  As order selection criteria, this function computes approximations to
%  Schwarz's Bayesian Criterion and to the logarithm of Akaike's Final
%  Prediction Error. The order selection criteria for models of order
%  pmin:pmax are returned as the vectors SBC and FPE.
%
%
%  If the optional argument SELECTOR is included in the function call,
%  as in ARFIT(v,pmin,pmax,SELECTOR), SELECTOR is used as the order
%  selection criterion in determining the optimum model order. The
%  three letter string SELECTOR must have one of the two values 'sbc'
%  or 'fpe'. (By default, Schwarz's criterion SBC is used.) If the
%  bounds pmin and pmax coincide, the order of the estimated model
%  is p=pmin=pmax.

% n:   number of time steps (per realization)
% m:   number of variables (dimension of state vectors)
% ntr: number of realizations (trials)
[n,m,ntr]   = size(v);

if (pmin ~= round(pmin) | pmax ~= round(pmax))
    error('Order must be integer.');
end
if (pmax < pmin)
    error('PMAX must be greater than or equal to PMIN.')
end

% set defaults and check for optional arguments
if (nargin == 3)              % no optional arguments => set default values
    mcor       = 1;           % fit intercept vector
    selector   = 'sbc';	      % use SBC as order selection criterion
elseif (nargin == 4)          % one optional argument
    if strcmp(selector, 'zero')
        mcor     = 0;               % no intercept vector to be fitted
        selector = 'sbc';	          % default order selection
    else
        mcor     = 1; 		  % fit intercept vector
    end
end

ne  	= ntr*(n-pmax);         % number of block equations of size m
npmax	= m*pmax+mcor;          % maximum number of parameter vectors of length m

if (ne <= npmax)
    error('Time series too short.')
end

% compute QR factorization for model of order pmax
[R, scale]   = arqr(v, pmax, mcor);

% compute approximate order selection criteria for models
% of order pmin:pmax
[sbc, fpe]   = arord(R, m, mcor, ne, pmin, pmax);

% get index iopt of order that minimizes the order selection
% criterion specified by the variable selector
[val, iopt]  = min(eval(selector));

% select order of model
popt         = pmin + iopt-1; % estimated optimum order
np           = m*popt + mcor; % number of parameter vectors of length m