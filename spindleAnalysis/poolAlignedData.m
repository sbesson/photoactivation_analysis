function [alignedTimes, alignedData] = poolAlignedData(s)
% Pool aligned times and data from a series of measurements
%
% SYNOPSIS
%
% This function first reads the dimenesions
%
% INPUT
%     s - an array of structure with 2  fields: alignedTimes and
%     alignedData
%
% OUTPUT
%     alignedTimes - a nTimePointsMax x 1 array of pooled aligned times
%
%     alignedData - a nTimePointsMax x nSamples matrix of pooled aligned

% Read number of series and dimenesions
nSeries = numel(s);
[nTimePoints, nSamples] = cellfun(@size, {s.alignedData});
nSamples = [0 cumsum(nSamples)];

% Get the largest alignedTimes series
[nTimePoints_max, imax] = max(nTimePoints);
alignedTimes = s(imax).alignedTimes;

% Initialize alignedData matrix
alignedData = NaN(nTimePoints_max, nSamples(end));

for i = 1 : nSeries
    % Find the intersection with the times
    [~, ind] = intersect(alignedTimes, s(i).alignedTimes);
    assert(numel(ind) == numel(s(i).alignedTimes),...
        'Aligned times do not correspond');
    
    % Fill the alignedData matrix
    alignedData(ind, nSamples(i)+1:nSamples(i+1)) = s(i).alignedData;
end