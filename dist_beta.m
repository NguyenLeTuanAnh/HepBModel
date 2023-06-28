%DIST_BETA   Construct a beta probability distribution.
%   The distribution is specified by the mean, variance and minimum/maximum
%   values.
%
%   dist_fn = DIST_BETA(mean, variance, min_v, max_v);
%
%   MEAN is the mean value for the beta distribution, over the (pre-scaling)
%   interval [0, 1].
%
%   VARIANCE is the variance for the beta distribution, over the (pre-scaling)
%   interval [0, 1].
%
%   MIN_V is the minimum value for the scaled beta distribution.
%
%   MAX_V is the maximum value for the scaled beta distribution.
%
%   DIST_FN is a function that wraps the given beta distribution.
%
%See also dist_wrapper
function [dist_fn] = dist_beta(mean, variance, min_v, max_v)
    % Select a random value with a beta distribution,
    % where linear_val is in the unit interval
    function [v] = choose_beta(alpha_1, alpha_2, min_val, max_val, linear_val)
        range = max_val - min_val;
        v = min_val + range * betainv(linear_val, alpha_1, alpha_2);
    end

    % scale the mean and variance to the [0,1] beta distribution
    range = max_v - min_v;
    scaled_mean = (mean - min_v) / range;
    scaled_var = variance / range;

    % calculate the distribution parameters alpha_1 and alpha_2
    alpha_1 = (scaled_mean/scaled_var)^2 * (1 - scaled_mean) - scaled_mean;
    alpha_2 = alpha_1 * (1 - scaled_mean) / scaled_mean;

    % return the beta distribution function
    fn = @(val) choose_beta(alpha_1, alpha_2, min_v, max_v, val);
    dist_details = struct( ...
        'name', 'beta', 'mean', mean, 'variance', variance, ...
        'scaled_mean', scaled_mean, 'scaled_var', scaled_var, ...
        'range', range, ...
        'minimum', min_v, 'maximum', max_v, ...
        'alpha_1', alpha_1, 'alpha_2', alpha_2 ...
        );

    dist_fn = dist_wrapper(fn, dist_details);
end
