%LAYOUT_OF   Determine how to translate between a record and a column vector.
%   A layout defines how a record whose fields contain numerical arrays can be
%   converted into a column vector and back into an identical record. This is
%   used by the LHS framework to allow models to treat the model states and
%   parameters as non-scalars, without having to manually convert these values
%   into a single for the Matlab equation solvers.
%
%   A layout is stored as a record (struct) with fieldnames identical to the
%   original record. The contents of each field is a vector that contains the
%   number of elements in the matrix associated with the field, and the
%   dimensions of this matrix. For example, if record.S is a 3x2 matrix. then
%   layout.S is [6 3 2], since 6 = 3 x 2.
%
%   A layout also has an additional field, 'size', which records the total
%   number of elements in all of the matrices. Therefore, there must not be a
%   field called 'size' in the original record (this will raise an exception).
%
%   layout = LAYOUT_OF(record);
%
%   RECORD is a structure where each field contains an arbitrary matrix.
%
%   LAYOUT is a structure that defines how to convert data between the format
%   in which it is stored in RECORD and a column vector.
%
%   Example:
%      >> s1 = struct('a', [1 2 3], 'b', [5 6; 8 9]);
%      >> layout = layout_of(s1);
%      % Convert s1 to a vector
%      >> vec = layout_vec(layout, s1);
%      % Convert from a vector to a record
%      >> s2 = layout_rec(layout, vec);
%      >> isequal(s1, s2)
%      ans = 1
%
%See also layout_vec, layout_rec, layout_time_rec, layout_merge
function layout = layout_of(record)
    % Build the layout record and vector size
    layout = struct(record);
    if isfield(layout, 'size')
        error('layout_of:size', ...
            'The field "size" is reserved and cannot be used');
    end
    vec_size = 0;
    statenames = fieldnames(record);
    for i = 1:length(statenames)
        state = statenames{i};
        state_size = size(record.(state));
        state_len = prod(state_size);
        layout.(state) = [state_len state_size];
        vec_size = vec_size + state_len;
    end
    layout.size = vec_size;
end
