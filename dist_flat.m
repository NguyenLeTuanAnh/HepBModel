%DIST_FLAT   Construct a flat (uniform) probability distribution.
%   The distribution is specified by the minimum and maximum values.
%
%   dist_fn = DIST_FLAT(min_val, max_val);
%
%   MIN_VAL is the minimum value for the uniform distribution.
%
%   MAX_VAL is the maximum value for the uniform distribution.
%
%   DIST_FN is a function that wraps the given uniform distribution.
%
%See also dist_wrapper
function [dist_fn] = dist_flat(min_val, max_val)
    fn = @(x) (min_val + (max_val - min_val) * x) .* ones(1, length(x));
    dist_details = struct( ...
        'name', 'flat', 'minimum', min_val, 'maximum', max_val ...
        );
    dist_fn = dist_wrapper(fn, dist_details);
end
