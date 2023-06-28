%FILENAME_TIMESTAMPED   Adds a timestamp to filenames.
%   This function takes a basename and extension for a filename and returns a
%   filename that includes the current date and time.
%
%   filename = FILENAME_TIMESTAMPED(basename, extension);
%
%   BASENAME is the name of the file, excluding any extension.
%
%   EXTENSION is the extension of the file, including the leading period. This
%   argument is optional and defaults to ''.
%
%   FILENAME is the new filename that includes the current date and time.
%
%   Example:
%      >> filename_timestamped('test_run', '.eps')
%      'test_run_2009-09-16_13-20.eps'
%
%See also filename, most_recent
function [name] = filename_timestamped(basename, extension)
    c = clock;
    year = sprintf('%02d', c(1));
    mon = sprintf('%02d', c(2));
    day = sprintf('%02d', c(3));
    hr = sprintf('%02d', c(4));
    mins = sprintf('%02d', c(5));

    if nargin < 2
        extension = '';
    end

    name = strcat(basename, '_', year, '-', mon, '-', day, ...
        '_', hr, '-', mins, extension);
end
