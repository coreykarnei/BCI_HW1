function warning2error(cmd, varargin)
% WARNING2ERROR - Translate warnings into errors
%
% warning2error(cmd)
%
% Where 
%
% CMD is a string containing a MATLAB statement. Such statement will be
% evaluated by this function and any generated warning will be translated
% into a warning having the same ID and message.
%
% ## Notes:
%
% An alternative solution is to use the following undocumented MATLAB
% feature:
%
% http://undocumentedmatlab.com/blog/trapping-warnings-efficiently/
%
% However, the solution above does not work always. In particular it does
% not seem to catch some warning generated by some closed-source MATLAB
% built-ins. This function should be more robust.
%
% See also: warning


warning('off', 'dummy:dummy');
warning('dummy:dummy', 'dummy');

s = warning('query', 'all');

if nargin > 1,
    
    for i = 1:numel(varargin),
        warning('off', varargin{i});
    end
    
else
    
    warning('off', 'all');    
    
end

try
    
    evalin('caller', cmd);
    [msg, id] = lastwarn;
    if ~strcmp(id, 'dummy:dummy'),
        ME = MException(id, msg);
        throw(ME);
    end
    warning(s);    
    
catch ME    
        
    warning(s); 
    rethrow(ME);
    
end

end