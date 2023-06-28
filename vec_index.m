%VEC_INDEX   Get the scalar index of an element in an N-dimensional matrix.
%   This function converts N-dimensional coordinates into the corresponding
%   scalar offset into the single-column vector of matrix entries.
%
%   index = VEC_INDEX(dims, coords, offset);
%
%   DIMS is a vector that defines the matrix dimensions, as returned by SIZE().
%
%   COORDS is a vector of coordinates the specify an entry in the matrix.
%
%   OFFSET is an initial position in the single-column vector from which the
%   scalar index should be calculated. This is an optional argument and has a
%   default value of 1.
%
%   INDEX is the position of the matrix entry in the single-column vector.
%
%   Example:
%      >> M(1, 2, 3) = 99;
%      >> V = M(:);
%      >> i = vec_index(size(M), [1, 2, 3]);
%      >> V(i)
%      ans = 99
function index = vec_index(dims, coords, offset)
    if nargin < 3
        offset = 1;
    end
    index = offset + sum( (coords - 1) .* [1 cumprod(dims(1:end-1))]);
end
