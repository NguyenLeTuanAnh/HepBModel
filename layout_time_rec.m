%LAYOUT_TIME_REC   Translate a matrix of time series into a record.
%   Given a pre-calculated layout of how to convert the data, LAYOUT_TIME_REC
%   translates time-series data into a record of matrices, where the first
%   dimension of each matrix is time.
%
%   This function is used by the LHS framework to convert time-series data into
%   a record format, for the purposes of calculating statistics from a single
%   model simulation.
%
%   rec = LAYOUT_TIME_REC(layout, time_vec);
%
%   LAYOUT is the result of a previous call to LAYOUT_OF.
%
%   TIME_VEC is a matrix containing the time series to be converted.
%
%   REC is a record of arbitrary matrices.
%
%See also layout_of, layout_vec, layout_rec, layout_merge
function rec = layout_time_rec(layout, time_vec)
    % NOTE: The first dimension of time_vec is time
    time_steps = size(time_vec, 1);

    % Build the state record from the state vector and state layout
    rec = struct();
    position = 1;
    statenames = fieldnames(rmfield(layout, 'size'));
    for i = 1:length(statenames)
        state = statenames{i};
        state_layout = layout.(state);
        state_length = state_layout(1);
        state_values = time_vec(:, position:position + state_length - 1);
        state_dims = state_layout(2:end);
        if state_dims(1) == 1
            state_dims = state_dims(2:end);
        end
        rec.(state) = reshape(state_values, [time_steps state_dims]);
        position = position + state_length;
    end
end
