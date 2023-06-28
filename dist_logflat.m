%DIST_LOGFLAT   Construct a log-flat (uniform) probability distribution.
%   The distribution is specified by the minimum and maximum values.
%
%   dist_fn = DIST_LOGFLAT(min_exp, max_exp);
%
%   MIN_EXP is the minimum value for the exponent of the log distribution.
%
%   MAX_EXP is the maximum value for the exponent of the log distribution.
%
%   DIST_FN is a function that wraps the given log-uniform distribution.
%
%See also dist_wrapper
function [dist_fn] = dist_logflat(min_exp, max_exp)
    fn = @(x) (10^(min_exp+(max_exp-min_exp)*x)) .* ones(1, length(x));
    dist_details = struct( ...
        'name', 'logflat', 'minimum', min_exp, 'maximum', max_exp ...
        );
    dist_fn = dist_wrapper(fn, dist_details);
end
