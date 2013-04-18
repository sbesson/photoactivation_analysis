function alignedData = alignData(data, index)
% Align series of data with respect to an index array
%
% SYNOPSIS alignedData = alignData(data, index)
%
% This function align input time-series with respect to a specified index
% series
%
% INPUT
%       data - a nTimepoints x nSamples matrix containing the time-series
%
%       index - a 1 x nSamples array of indexes to align th data
%
% OUTPUT
%      aligneData - a 2nTimepoints x nSamples matrix


% Sebastien Besson Nov 2012 (last modified Apr 2012)

% Input check
ip = inputParser;
ip.addRequired('data',@isnumeric);
ip.addRequired('index',@(x) isvector(x) && numel(x) == size(data, 2));
ip.parse(data, index);


%% Data alignment
nTimePoints = size(data, 1);
validIndex = find(~isnan(index));
nSamples = numel(validIndex);

% Create matrix of size 2n x m
alignedData = NaN(2*nTimePoints, nSamples);
for i = 1 : nSamples
    tmin = nTimePoints - index(validIndex(i)) + 1;
    tmax = 2 * nTimePoints - index(validIndex(i));
    alignedData(tmin : tmax, i) =  data(:,i);
end