%FILENAME   Generates a filename for saving experimental results and plots.
%   This function generates filenames relative to the location of the directory
%   'experiments', which MUST be included in the MATLAB path. Filenames can be
%   timestamped and organised into a hierarchy of subdirectories. This function
%   can also be used to generate paths, by only providing the first argument.
%
%   [filename filepath] = FILENAME(subdirs, basename, extension, timestamp);
%
%   SUBDIRS is a string that names a subdirectory, or a cell array of nested
%   subdirectories.
%
%   BASENAME is the name of the file, excluding any extension.
%
%   EXTENSION is the extension of the file, including the leading period.
%
%   TIMESTAMP is a logical value that specifies whether the current date and
%   time should be incorporated into the filename. This is of particular use
%   when saving simulation results.
%
%   FILENAME is the (absolute) filename.
%
%   FILEPATH is the (absolute) path of the directory that contains the file.
%
%   Example:
%      >> filename('data', 'some_experiment', '.mat', 1);
%
%See also filename_timestamped, most_recent
function [filename filepath] = filename(subdir, filename, ext, timestamp)
    if nargin < 1
        subdir = ''; % default to no subdirectory
    end
    if nargin < 2
        filename = ''; % default to no filename
    end
    if nargin < 3
        ext = ''; % default to no file extension
    elseif ~ strncmp(ext, '.', 1)
        % ensure the extension begins with a period
        warning('filename:extension', ...
            'Leading period added to filename extension "%s"', ext);
        ext = ['.' ext];
    end
    if nargin < 4
        timestamp = 0; % default to no timestamp in the filename
    end

    % find the path to the directory that contains the model experiments
    % NOTE: all paths are constructed relative to this directory
    s = what('lhsexperiments');
    if length(s) < 1
        error('filename:path', 'Unable to locate the "experiments" directory');
    elseif length(s) > 1
        warning('filename:path', ...
            'Multiple "experiments" directories in the MATLAB path');
    end

    % cater for nested subdirectories
    if iscellstr(subdir)
        subdir = fullfile(subdir{:});
    end

    % generate the absolute filename, including the subdirectory
    if strcmp(subdir, '')
        filepath = s(1).path;
    else
        filepath = fullfile(s(1).path, subdir);
        % create the subdirectory if it does not exist
        if ~ isdir(filepath)
            status = mkdir(filepath);
            if ~ status
                error('filename:mkdir', ...
                    'Unable to create directory "%s"', filepath);
            end
        end
    end
    file_base = fullfile(filepath, filename);

    % return a time-stamped filename
    if timestamp
        filename = filename_timestamped(file_base, ext);
    else
        filename = [file_base ext];
    end
end
