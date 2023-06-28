%DIST_CONST   Construct a constant probability distribution.
%   The distribution is specified by a single value.
%
%   dist_fn = DIST_CONST(value);
%
%   VALUE is the fixed value returned by the constant distribution.
%
%   DIST_FN is a function that wraps the given constant distribution.
%
%See also dist_wrapper
function [dist_fn] = dist_const(value)
    fn = @(x) value;
    dist_details = struct( ...
        'name', 'constant', 'value', value ...
        );
    dist_fn = dist_wrapper(fn, dist_details);
end
