function alignedData = alignData(data, index)
% Align series of data with respect to an index array
%
% SYNOPSIS alignedData = alignData(data, index)
%
% This function align input time-series with respect to a specified index
% series
%
% INPUT
%       data - a n x m matrix containing the time-series
%
%       index - a 1 x m
%
% OUTPUT
%      aligneData - a 2n x m matrix


% Sebastien Besson Nov 2012

% Input check
ip = inputParser;
ip.addRequired('data',@isnumeric);
ip.addRequired('index',@(x) isvector(x) && numel(x) == size(data, 2));
ip.parse(data, index);


%% Data alignment
nTimePoints = size(data, 1);
nSamples = size(data, 2);

% Create matrix of size 2n x m
alignedData = NaN(2*nTimePoints, nSamples);
for i=1:nSamples
    alignedData(nTimePoints - index(i) + 1: 2* nTimePoints - index(i),i) =  data(:,i);
end