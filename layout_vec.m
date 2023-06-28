%LAYOUT_VEC   Translate a record of values into a column vector.
%   Given a pre-calculated layout of how to convert the data, LAYOUT_VEC
%   translates the data into a column vector.
%
%   vec = LAYOUT_VEC(layout, rec);
%
%   LAYOUT is the result of a previous call to LAYOUT_OF.
%
%   REC is a record containing the data to be converted.
%
%   VEC is a column vector of scalar values.
%
%See also layout_of, layout_rec, layout_time_rec, layout_merge
function vec = layout_vec(layout, rec)
    vec = zeros([layout.size 1]);
    position = 1;
    statenames = fieldnames(rec);
    for i = 1:length(statenames)
        state = statenames{i};
        state_values = rec.(state);
        state_layout = layout.(state);
        state_length = state_layout(1);
        vec(position:position + state_length - 1) = state_values;
        position = position + state_length;
    end
end
