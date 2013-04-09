function plotBarGraphs(data, names, ylabel, varargin)

% define small and large fonts for graphical output
default_tfont = {'FontName', 'Helvetica', 'FontSize', 14, 'FontAngle', 'italic'};
default_sfont = {'FontName', 'Helvetica', 'FontSize', 18};
default_lfont = {'FontName', 'Helvetica', 'FontSize', 22};

% Input checl
ip = inputParser;
ip.addRequired('data',@iscell);
ip.addRequired('names',@iscell);
ip.addRequired('ylabel',@ischar);
ip.addParamValue('tfont', default_tfont, @iscell);
ip.addParamValue('sfont', default_sfont, @iscell);
ip.addParamValue('lfont', default_lfont, @iscell);
ip.parse(data, names, ylabel, varargin{:});

% Plot indidividual events

means = cellfun(@nanmean, data);
std = cellfun(@(x)nanstd(x,[],2), data);
nPoints = cellfun(@numel, data);
ste = std./sqrt(nPoints);

barplot2(means, ste, 'XTickLabel',names, 'YLabel', ylabel);
