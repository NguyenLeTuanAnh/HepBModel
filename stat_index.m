%STAT_INDEX   Get the scalar index of a single statistic.
%   When the statistics for a simulation are non-trivial (ie, one or more
%   statistics are non-scalar), STAT_INDEX should be used to index a single
%   element of a statistic. Given a statistic name and coordinates into the
%   matrix associated with that statistic, STAT_INDEX returns the corresponding
%   offset into the lhs_data.stats matrix. This is of use when using
%   PLOT_STATISTIC (see the example below).
%
%   index = STAT_INDEX(layout, stat_name, coords);
%
%   LAYOUT is either a layout structure (as returned by LAYOUT_OF) or the
%   results of an LHS experiment, where the layout of the statistics must be
%   stored in the 'stat_layout' field.
%
%   STAT_NAME is the name of the statistic in question.
%
%   COORDS is the coordinates of the single element from this statistic whose
%   index is to be returned. Coordinates must be a vector whose elements
%   specify a valid index into the matrix associated with the statistic.
%
%   INDEX is the scalar index of the statistic, which can be passed to
%   PLOT_STATISTIC.
%
%   Example:
%      % Load the results of some LHS experiment
%      >> results = most_recent('some_experiment');
%      % Get the index of the statistic I_max(1,2)
%      >> index = stat_index(results.lhs_data, 'I_max', [1 2]);
%      % Plot the maximum value of this statistic over the parameter space
%      >> plot_statistic(index, @max, results.lhs_data);
%
%See also vec_index, layout_of, plot_statistic
function index = stat_index(layout, stat_name, coords)
    % accept an experimental result structure and extract the stat layout
    % otherwise, the first argument must be a valid layout
    if isfield(layout, 'stat_layout')
        layout = layout.stat_layout;
    end

    offset = 1;
    % ignore the 'size' field
    stat_names = fieldnames(rmfield(layout, 'size'));
    for i = 1:length(stat_names)
        name = stat_names{i};
        if strcmp(name, stat_name)
            % we have found the offset to this statistic
            break;
        else
            % increase the offset by the number of elements associated with the
            % statistic
            stat_details = layout.(name);
            offset = offset + stat_details(1);
        end
    end
    dims = layout.(stat_name);
    index = vec_index(dims(2:end), coords, offset);
end
