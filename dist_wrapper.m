%DIST_WRAPPER   Combine a distribution with distribution parameter values.
%   This function is used to wrap a probability distribution inside a separate
%   function that allows other MATLAB code to randomly select values from the
%   distribution, but also to retrieve the distribution type (eg, beta) and
%   distribution parameters (eg, mean).
%
%   dist_fn = DIST_WRAPPER(fn_handle, dist_details);
%
%   FN_HANDLE is a handle to a function that maps the interval [0, 1] to some
%   range of values, according to a probability distribution (eg, see
%   DIST_BETA).
%
%   DIST_DETAILS is a record (structure array) that contains parameters
%   specific to the given probability distribution. The type of distribution is
%   uniquely identified by the contents of the 'name' field.
%
%   DIST_FN is a function that wraps the given probability distribution. When
%   called with a scalar argument in the interval [0, 1], it returns the
%   corresponding value from the probability distribution. When called with no
%   arguments, it returns a record (structure array) that contains the
%   distribution type and the parameter values.
%
%   Example:
%      % Map the [0,1] interval to some probability distribution
%      >> fn = @some_function;
%      % Record the parameters of this distribution
%      >> dist_details = struct('name', 'some_distribution', 'param1', value1);
%      % Create (and wrap) the probability distribution
%      >> dist_fn = dist_wrapper(fn, dist_details);
%      % Select a random value from the distribution
%      >> random_value = dist_fn( rand() );
%      % Retrieve the distribution parameters
%      >> dist_details = dist_fn();
%
%See also dist_beta, dist_const, dist_flat, dist_normal, dist_generic
function [dist_fn] = dist_wrapper(fn_handle, dist_details)
    function [rand_val] = call_dist(linear_val)
        if nargin == 1
            if ((linear_val < 0) || (linear_val > 1))
                error('Argument must be in [0, 1]');
            else
                rand_val = fn_handle(linear_val);
            end
        else
            rand_val = dist_details;
        end
    end
    dist_fn = @call_dist;
end
