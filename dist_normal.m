%DIST_NORMAL   Construct a normal (Gaussian) probability distribution.
%   The distribution is specified by the mean and standard deviation.
%
%   dist_fn = DIST_NORMAL(mean, stddev);
%
%   MEAN is the mean value for the normal distribution.
%
%   STDDEV is the standard deviation for the normal distribution.
%
%   DIST_FN is a function that wraps the given normal distribution.
%
%See also dist_wrapper
function [dist_fn] = dist_normal(mean, stddev)
    % Select a random value with a Gaussian distribution,
    % where linear_val is in the unit interval
    function v = choose_gaussian(mean, stddev, linear_val)
        v = norminv(linear_val, mean, stddev);
    end

    fn = @(val) choose_gaussian(mean, stddev, val);
    dist_details = struct( ...
        'name', 'normal', 'mean', mean, 'stddev', stddev ...
        );
    dist_fn = dist_wrapper(fn, dist_details);
end
