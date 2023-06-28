%DIST_MANY   Combine multiple probability distributions.
%   The distribution is a vector of probability distributions, typically of use
%   when a single state in the model is divided into multiple strata.
%
%   dist_fn = DIST_MANY(dists);
%
%   DISTS is a cell array of the probability distributions to combine together.
%
%   DIST_FN is a function that wraps the given distributions.
%
%See also dist_wrapper
function [dist_fn] = dist_many(dists)
    dist_count = length(dists);
    if dist_count < 1
        error('dist_many:dists', 'At least one distribution is required');
    elseif dist_count == 1
        warning('dist_many:dists', 'Only one distribution is given');
    end

    function v = choose_values(val, dists)
        % choose random values from an N-dimensional matrix of distributions
        if isscalar(val)
            values = cell(size(val));
            values{:} = val;
        elseif size(val) ~= size(dists)
            error('dist_many:choose', ...
                'Size mismatch between values and distributions')
        elseif isnumeric(val)
            values = num2cell(val);
        elseif iscell(val)
            values = val;
        else
            error('dist_many:choose', 'Invalid input');
        end
        cell_v = cellfun(@(d, v) d(v), dists, values, 'UniformOutput', 0);
        v = cell2mat(cell_v);
    end
    
    % collect the details of each distribution and store in a cell array
    dist_details = cellfun(@(dist) dist(), dists, 'UniformOutput', 0);

    fn = @(val) choose_values(val, dists);
    dist_fn = dist_wrapper(fn, dist_details);
end
