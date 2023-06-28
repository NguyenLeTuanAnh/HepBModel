%LAYOUT_MERGE   Merge two layouts together.
%   This function returns a layout that is the combination of two existing
%   layouts. This is used by the LHS framework when dealing with sampling the
%   independent and dependent parameters, as these two parameter sets must be
%   treated both separately and together in various parts of the LHS code.
%
%   layout = LAYOUT_MERGE(layout1, layout2);
%
%   LAYOUT1, LAYOUT2 are two layouts returned by previous calls to LAYOUT_OF.
%   They must not share any field names (apart from 'size') or an exception
%   will be raised.
%
%   LAYOUT is the merged layout.
%
%See also layout_of, layout_vec, layout_rec, layout_time_rec
function layout = layout_merge(layout1, layout2)
    total_size = layout1.size + layout2.size;
    layout1 = rmfield(layout1, 'size');
    layout2 = rmfield(layout2, 'size');

    M = [fieldnames(layout1)'  fieldnames(layout2)'; ...
         struct2cell(layout1)' struct2cell(layout2)'];
    % this will throw an exception if layout1 and
    % layout2 have any parameter names in common
    layout = struct( M{:} );
    layout.size = total_size;
end
