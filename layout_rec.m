%LAYOUT_REC   Translate a column vector of values into a record.
%   Given a pre-calculated layout of how to convert the data, LAYOUT_REC
%   translates the data into a record.
%
%   rec = LAYOUT_REC(layout, vec);
%
%   LAYOUT is the result of a previous call to LAYOUT_OF.
%
%   VEC is a column vector containing the data to be converted.
%
%   REC is a record of arbitrary matrices.
%
%See also layout_of, layout_vec, layout_time_rec, layout_merge
function rec = layout_rec(layout, vec)
    % Build the state record from the state vector and state layout
    rec = struct();
    position = 1;
    statenames = fieldnames(rmfield(layout, 'size'));
    for i = 1:length(statenames)
        state = statenames{i};
        state_layout = layout.(state);
        state_length = state_layout(1);
        state_values = vec(position:position + state_length - 1);
        rec.(state) = reshape(state_values, state_layout(2:end));
        position = position + state_length;
    end
end
