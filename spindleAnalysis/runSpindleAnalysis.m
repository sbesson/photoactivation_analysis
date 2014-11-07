clear
clc
close all

% Define constants
%dataPath = uigetdir('Select the main directory');
dataPath = '/Users/sebastien/Desktop/Ammended data - accounts for EMCCD swap';
if isequal(dataPath, 0), return; end

% Create output folder for aligned data
alignedDataPath = fullfile(dataPath, 'AlignedData');

% Read conditions from folder names
directories  = dir(dataPath);
directories = directories([directories.isdir]);
conditionsDirs = {directories.name};
dotDirs = strcmp(conditionsDirs, '.')  | strcmp(conditionsDirs, '..');
conditionsDirs(dotDirs) = [];
conditionsDirs(strcmp(conditionsDirs, 'AlignedData')) = [];

% Initialize output structure array
conditions =  cell2struct(conditionsDirs, 'name');
nConditions = numel(conditions);
fileLog = cell(nConditions, 1);

% Create output folder for aligned data
if ~isdir(alignedDataPath), mkdir(alignedDataPath); end

% Alignment parameters
minSpindleLength = 5; % minimum spindle length (in microns)
sigma = 3; % Sigma to exclude outliers (using final spindle length)

% Define mapping for bar graphs
graph_map{1,1} = 'elongationRate';
graph_map{1,2} = 'Elongation rate (microns/min)';
graph_map{2,1} = 'finalLengths';
graph_map{2,2} = 'Final spindle length (microns)';
nGraphs = size(graph_map, 1);

for iCondition = 1:nConditions
    
    % Find Excel file in condition folder
    conditionPath = fullfile(dataPath, conditions(iCondition).name);
    conditions(iCondition).path = conditionPath;
    xlsfiles = dir(fullfile(conditionPath, '*.xls*'));
    assert(numel(xlsfiles) == 1, 'More than one XLS file in folder %s',...
        conditionPath);
    xlsfile = fullfile(conditionPath, xlsfiles(1).name);
    
    %% Read data
    fileLog{iCondition} = sprintf('Reading file %s...\n', xlsfile);
    
    % Read number of worksheets
    [status, series] = xlsfinfo(xlsfile);
    nSeries = numel(series);
    fileLog{iCondition} = [fileLog{iCondition}...
        sprintf('Found %g worskheets\n', nSeries)];
    seriesLog = cell(nSeries, 1);
    
    for iSeries = 1 : nSeries
        %% Read raw data
        [data, str] = xlsread(xlsfile, series{iSeries});
        if isempty(data), continue; end
        conditions(iCondition).series(iSeries).name = series{iSeries};
        seriesLog{iSeries} = sprintf('\nReading series %g/%g: %s\n',...
            iSeries, nSeries, series{iSeries});
        
        data(1,:) = []; % Remove first row
        
        conditions(iCondition).series(iSeries).originaldata = data; % Log original data
        
        % Extract first column (time points)
        times = data(:,1);
        dt = unique(diff(times));
        assert(isscalar(dt), 'File %s, worskheet %s: Times not evenly spaces', ...
            xlsfile, series{iSeries});
        conditions(iCondition).series(iSeries).times = times;
        data(:,1)=[];
        data(:,all(isnan(data),1))=[];
        
        %% Outlier detection
        
        % Compute final spindle lengths
        finalLengths = NaN(1, size(data,2));
        for i = 1:size( data,2)
            finalLengths(i) = data(find(~isnan( data(:,i)), 1, 'last'), i);
        end
        log = sprintf('Final spindle length: %g +/- %g\n', mean(finalLengths),...
            std(finalLengths));
        seriesLog{iSeries} = [seriesLog{iSeries} log];
        
        % Detect outliers using final spindle length
        outlierIndex = detectOutliers(finalLengths, sigma);
        conditions(iCondition).series(iSeries).outlierIndex = outlierIndex;
        
        % Remove outliers
        if ~isempty(outlierIndex)
            seriesLog{iSeries} = [seriesLog{iSeries}...
                'Removing outliers: ' sprintf('%g ', outlierIndex) sprintf('\n')];
            data(:, outlierIndex)=[];
            finalLengths(outlierIndex) = [];
        end
        
        conditions(iCondition).series(iSeries).finalLengths = finalLengths;
        conditions(iCondition).series(iSeries).data = data;
        
        % Log number of time points and samples
        nTimePoints = size(data,1);
        nSamples = size(data,2);
        log1 = sprintf('Number of timepoints: %g\n', nTimePoints);
        log2 = sprintf('Number of samples: %g\n', nSamples);
        seriesLog{iSeries} = [seriesLog{iSeries} log1 log2];
        
        %% Detect events
        [events, eventlog] = detectEvents(data, times, minSpindleLength);
        seriesLog{iSeries} = [seriesLog{iSeries} eventlog sprintf('\n')];
        
        conditions(iCondition).series(iSeries).events = events;
        
        
        
        %% Measure elongation rates
        elongationRate = @(x) (x(3).values- x(1).values)./(x(3).times- x(1).times);
        conditions(iCondition).series(iSeries).elongationRate = ...
            elongationRate(conditions(iCondition).series(iSeries).events);
        
        %% Align data
        seriesLog{iSeries} = [seriesLog{iSeries} ...
            'Aligning data using ' events(1).name sprintf('\n')];
        alignedData = alignData(data, events(1).index);
        alignedTimes = (-nTimePoints+1:nTimePoints)' * dt;
        conditions(iCondition).series(iSeries).alignedData = alignedData;
        conditions(iCondition).series(iSeries).alignedTimes  = alignedTimes;
        
        
    end
    
    %% Create comparative graphs per condition
    
    for iGraph = 1 : nGraphs
        % Aggregate scalar data
        conditions(iCondition).(graph_map{iGraph, 1}) = ...
            horzcat(conditions(iCondition).series.(graph_map{iGraph, 1}));
    end
    
    % Pool aligned data
    [alignedTimes, alignedData] = ...
        poolAlignedData(conditions(iCondition).series);
    conditions(iCondition).alignedTimes = alignedTimes;
    conditions(iCondition).alignedData = alignedData;
    
    % Export aligned series into TXT format
    dlmwrite(fullfile(alignedDataPath, [conditions(iCondition).name '-alignedData-raw.txt']),...
        [conditions(iCondition).alignedTimes...
        conditions(iCondition).alignedData], '\t');
    
    % Export aligned series into TXT format
    alignedmean = nanmean(conditions(iCondition).alignedData,2);
    alignedstd = nanstd(conditions(iCondition).alignedData, [], 2);
    nPoints  =sum(~isnan(conditions(iCondition).alignedData),2);
    alignedste = alignedstd./sqrt(nPoints);
    dlmwrite(fullfile(alignedDataPath, [conditions(iCondition).name '-alignedData-MeanStdSte.txt']),...
        [conditions(iCondition).alignedTimes alignedmean alignedstd...
        alignedste], '\t');
    
    % Display and save full log
    disp([fileLog{iCondition} seriesLog{:}])
    fid = fopen(fullfile(conditionPath, 'Summary-log.txt'), 'w');
    fprintf(fid, [fileLog{iCondition} seriesLog{:}]);
    fclose(fid);
end

% Save output
save(fullfile(dataPath, 'analysis.mat'), 'conditions');


%% Graphic properties
conditionNames = {conditions.name};
% Conditions to be plotted, set to 1 : numel(conditions) to plot all graphs
conditions2plot = 1:numel(conditions);
% Maximum value of the spindle length
ymax = 18;
disp(sprintf('Generating graph for %s\n', conditionNames{conditions2plot}))

%% Plot individual series/event
rawFig = figure;
for iCondition = conditions2plot
    
    for iSeries = 1 : numel(conditions(iCondition).series)
        % Plot aligned data
        plotIndividualSeries(conditions(iCondition).series(iSeries).alignedTimes,...
            conditions(iCondition).series(iSeries).alignedData);
        if ~isempty(ymax), ylim([0 ymax]); end
        exportFigure(rawFig, conditions(iCondition).path,...
            [conditions(iCondition).name '-'...
            conditions(iCondition).series(iSeries).name '-aligneddata']);
        clf(rawFig)
        
        % Plot indidividual events
        plotIndividualSeries(conditions(iCondition).series(iSeries).times,...
            conditions(iCondition).series(iSeries).data,...
            conditions(iCondition).series(iSeries).events);
        if ~isempty(ymax), ylim([0 ymax]); end
        exportFigure(rawFig, conditions(iCondition).path,....
            [conditions(iCondition).name '-'...
            conditions(iCondition).series(iSeries).name '-events']);
        clf(rawFig)
    end
    
    % Plot aligned data per condition
    plotIndividualSeries(conditions(iCondition).alignedTimes,...
        conditions(iCondition).alignedData);
    if ~isempty(ymax), ylim([0 ymax]); end
    exportFigure(rawFig, conditions(iCondition).path,...
        [conditions(iCondition).name '-aligneddata']);
    clf(rawFig)
    
    % Set summary figure options
    plotMeanAlignedData(conditions(iCondition).series)
    xlim = get(gca, 'XLim');
    plot([xlim(1) xlim(2)], [minSpindleLength minSpindleLength] ,'--k');
    if ~isempty(ymax), ylim([0 ymax]); end
    exportFigure(rawFig, conditions(iCondition).path, 'AlignedData');
    clf(rawFig)
end
close(rawFig)

%% Generate comparative graphs between conditions
conditionFig = figure;
for iGraph = 1 : nGraphs
    for iCondition = conditions2plot
        % Plot summary bar graphs
        plotBoxPlots({conditions(iCondition).series.(graph_map{iGraph, 1})},...
            {conditions(iCondition).series.name}, graph_map{iGraph, 2});
        exportFigure(conditionFig, conditions(iCondition).path, graph_map{iGraph, 1});
        clf(conditionFig);
    end
    
    plotBoxPlots({conditions(conditions2plot).(graph_map{iGraph, 1})},...
        {conditions(conditions2plot).name}, graph_map{iGraph, 2});
    exportFigure(conditionFig, alignedDataPath, graph_map{iGraph, 1});
    clf(conditionFig);
end

% Generated all conditions comparison
plotMeanAlignedData(conditions(conditions2plot))
exportFigure(conditionFig, alignedDataPath, 'All conditions');
clf(conditionFig);

% % Generate 2 by 2 graph comparisons
% for i = conditions2plot
%     for j = 1 : i -1
%         plotMeanAlignedData(conditions([i,j]))
%         exportFigure(conditionFig, dataPath, [conditions(i).name ...
%             '-' conditions(j).name]);
%         clf(conditionFig);
%     end
% end
% close(conditionFig);

%%
for iCondition = 1:numel(conditions)
    prealignement =  conditions(iCondition).alignedData(conditions(iCondition).alignedTimes < 0,:);
    fprintf(1, '%s\n', conditions(iCondition).name);
    fprintf(1, '%g +/- %g um\n', nanmean(prealignement(:)), nanstd(prealignement(:), 1));
end
