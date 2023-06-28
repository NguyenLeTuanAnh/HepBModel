%DIST_FLATINT   Construct a flat (uniform) probability distribution on
%   an integer range.
%   The distribution is specified by the minimum and maximum values.
%
%   dist_fn = DIST_FLATINT(min_val, max_val);
%
%   MIN_VAL is the minimum value for the uniform integer distribution.
%
%   MAX_VAL is the maximum value for the uniform integer distribution.
%
%   DIST_FN is a function that wraps the given uniform distribution.
%
%See also dist_wrapper
function [dist_fn] = dist_flatint(min_val, max_val)
    fn = @(x) round((min_val-0.5 + (max_val+1 - min_val) * x));
    dist_details = struct( ...
        'name', 'flatint', 'minimum', min_val, 'maximum', max_val ...
        );
    dist_fn = dist_wrapper(fn, dist_details);
end
