%DIST_GENERIC   Construct a probability distribution defined by some CDF.
%   This function is used to construct probability distributions beyond those
%   provided explicitly by the LHS framework. Such distributions are defined by
%   providing an inverse cummulative distribution function, and any parameters
%   for the inverse CDF.
%
%   NOTE: Some inverse CDFs are not completely general (eg, betainv only
%   returns values in the [0,1] interval) and any scaling or shifting must be
%   performed manually; the code for DIST_BETA demonstrates how to do this.
%
%   dist_fn = DIST_GENERIC(name, inv_fn, params);
%
%   NAME is a string that identifies the distribution family when distribution
%   details are retrieved (see DIST_WRAPPER).
%
%   INV_FN is the inverse cummulative distribution function of the probability
%   distribution being constructed. The probability must be the first parameter
%   to this function. If the inverse CDF does not meet this condition, an
%   anonymous function must be used instead (see the example below).
%
%   PARAMS is a cell array of parameter names and values for the inverse CDF,
%   which must be provided in the same order as used by the inverse CDF. The
%   field names are purely descriptive, and are returned as part of the
%   distribution details (see DIST_WRAPPER).
%
%   DIST_FN is a function that wraps the given normal distribution.
%
%   Example:
%      % Define a normal distribution (mean of 5, standard deviation of 3)
%      >> dist_fn = dist_generic('normal', @norminv, {'mean' 5 'stddev' 3});
%      % Use an inverse CDF that takes parameters in a different order
%      >> inv_fn = @(probability, x1, x2) someinv(x1, x2, probability);
%      % Create a distribution with parameters x1=2 and x2=4
%      >> dist_fn = dist_generic('somedist', inv_fn, {'x1' 2 'x2' 4});
%
%See also dist_wrapper, dist_beta
function [dist_fn] = dist_generic(name, inv_fn, params)
    % extract the parameter names and values from the third argument
    param_names = params(1:2:end-1);
    param_values = params(2:2:end);
    % construct the function to return random values from the distribution
    fn = @(v) inv_fn(v, param_values{:});
    % construct the record that contains the details of this distribution
    details = struct();
    for i = 1:2:length(params)
        details.(params{i}) = params{i+1};
    end
    details.name = name;
    % order the fields so that 'name' is first, followed by the each of the
    % parameters, in the order that they are given.
    field_order = [{'name'} param_names];
    details = orderfields(details, field_order);
    % return the generic distribution function
    dist_fn = dist_wrapper(fn, details);
end
