clear
clc
close all

% Define constants
dataPath = uigetdir('Select the main directory');
% dataPath = '/Users/sebastien/Desktop/Bipolar spindle assembly assay';
if isequal(dataPath, 0), return; end

% Read conditions from folder names
directories  = dir(dataPath);
directories = directories([directories.isdir]);
conditionsDirs = {directories.name};
dotDirs = strcmp(conditionsDirs, '.')  | strcmp(conditionsDirs, '..');
conditionsDirs(dotDirs) = [];

% Initialize output structure array
conditions =  cell2struct(conditionsDirs, 'name');
nConditions = numel(conditions);
fileLog = cell(nConditions, 1);

% Alignment parameters
minSpindleLength = 6; % minimum spindle length (in microns)
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
        
        % Plot indidividual events
        figure;
        plotIndividualSeries(times, data, events);
        print(gcf,'-dtiff',fullfile(conditionPath,...
            [conditions(iCondition).name '-' series{iSeries} '-events.tif']));
        close(gcf)
        
        %% Measure elongation rates
        elongationRate = @(x) diff(vertcat(x.events.values),1)./diff(vertcat(x.events.times),1);
        conditions(iCondition).series(iSeries).elongationRate = ...
            elongationRate(conditions(iCondition).series(iSeries));
        
        %% Align data
        seriesLog{iSeries} = [seriesLog{iSeries} ...
            'Aligning data using ' events(1).name sprintf('\n')];
        alignedData = alignData(data, events(1).index);
        alignedTimes = (-nTimePoints+1:nTimePoints)' * dt;
        conditions(iCondition).series(iSeries).alignedData = alignedData;
        conditions(iCondition).series(iSeries).alignedTimes  = alignedTimes;
        
        % Plot aligned data
        figure;
        plotIndividualSeries(alignedTimes, alignedData);
        print(gcf,'-dtiff',fullfile(conditionPath,...
            [conditions(iCondition).name '-' series{iSeries} '-aligneddata.tif']));
        close(gcf)
    end
    
    %% Create comparative graphs per condition
    
    for iGraph = 1 : nGraphs
        % Aggregate scalar data
        conditions(iCondition).(graph_map{iGraph, 1}) = ...
            horzcat(conditions(iCondition).series.(graph_map{iGraph, 1}));
        
        % Plot summary bar graphs
        figure;
        plotBarGraphs({conditions(iCondition).series.(graph_map{iGraph, 1})},...
            {conditions(iCondition).series.name}, graph_map{iGraph, 2});
        print(gcf,'-dtiff',fullfile(conditionPath, [graph_map{iGraph, 1} '.tif']));
        close(gcf);
    end
    
    % Pool aligned data
    [alignedTimes, alignedData] = ...
        poolAlignedData(conditions(iCondition).series);
    conditions(iCondition).alignedTimes = alignedTimes;
    conditions(iCondition).alignedData = alignedData;
    
    % Set summary figure options
    conditionFig = figure();
    plotMeanAlignedData(conditions(iCondition).series)
    xlim = get(gca, 'XLim');
    plot([xlim(1) xlim(2)], [minSpindleLength minSpindleLength] ,'--k');
    print(gcf,'-dtiff',fullfile(conditionPath, 'AlignedData.tif'));
    close(conditionFig);   

    % Display and save full log
    disp([fileLog{iCondition} seriesLog{:}])
    fid = fopen(fullfile(conditionPath, 'Summary-log.txt'), 'w');
    fprintf(fid, [fileLog{iCondition} seriesLog{:}]);
    fclose(fid);
end

% Save output
save(fullfile(dataPath, 'analysis.mat'), 'conditions');

%% Generate comparative graphs between conditions

for iGraph = 1 : nGraphs
    figure;
    plotBarGraphs({conditions.(graph_map{iGraph, 1})},...
        {conditions.name}, graph_map{iGraph, 2});
    print(gcf,'-dtiff',fullfile(dataPath, [graph_map{iGraph, 1} '.tif']));
    close(gcf);
end

% Generated all conditions comparison
allConditionsFig = figure();
plotMeanAlignedData(conditions)
print(gcf,'-dtiff',fullfile(dataPath, 'All conditions.tif'));

% Generate 2 by 2 graph comparisons
for i = 1 : nConditions
    for j = 1 : i -1
        conditionFig = figure();
        plotMeanAlignedData(conditions([i,j]))
        print(gcf,'-dtiff',fullfile(dataPath, [conditions(i).name ...
            '-' conditions(j).name '.tif']));
        close(conditionFig);
    end
end
