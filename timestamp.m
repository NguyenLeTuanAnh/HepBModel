% TIMESTAMP   Generates a timestamp to be appended to filenames.
%   This function returns a stamp with the current date and time.
%
%   stamp = TIMESTAMP;
%
%   
%   Example:
%      >> timestamp
%      '_2009-09-16_13-20'
%
%
function [stamp] = timestamp
    c = clock;
    year = sprintf('%02d', c(1));
    mon = sprintf('%02d', c(2));
    day = sprintf('%02d', c(3));
    hr = sprintf('%02d', c(4));
    mins = sprintf('%02d', c(5));

    stamp = strcat('_', year, '-', mon, '-', day, ...
        '_', hr, '-', mins);
end
