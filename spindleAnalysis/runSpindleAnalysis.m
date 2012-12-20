clear
clc
close all

% Define constants
%dataPath = uigetdir('Select the main directory');
dataPath = '/Users/sebastien/Desktop/Bipolar spindle assembly assay';
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

for iCondition = 1:nConditions
    sheetFig = figure;
    
    % Find Excel file in condition folder
    conditionPath = fullfile(dataPath, conditions(iCondition).name);
    xlsfiles = dir(fullfile(conditionPath, '*.xls*'));
    assert(numel(xlsfiles) == 1, 'More than one XLS file in folder %s',...
        conditionPath);
    xlsfile = fullfile(conditionPath, xlsfiles(1).name);
    
    %%
    fileLog{iCondition} = sprintf('Reading file %s...\n', xlsfile);
    
    % Read number of worksheets
    [status, series] = xlsfinfo(xlsfile);
    nSeries = numel(series);
    fileLog{iCondition} = [fileLog{iCondition}...
        sprintf('Found %g worskheets\n', nSeries)];
    seriesLog = cell(nSeries, 1);
    pooledData = [];
    
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
        conditions(iCondition).series(iSeries).times = times;
        data(:,1)=[];
        data(:,all(isnan(data),1))=[];
        
        %% Outlier detection
        
        % Compute final spindle lengths
        finalLengths = NaN(size(data,2), 1);
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
        end
        
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
        
        %% Align data
        seriesLog{iSeries} = [seriesLog{iSeries} ...
            'Aligning data using ' events(1).name sprintf('\n')];
        alignedData = alignData(data, events(1).index);
        alignedTimes = (-nTimePoints+1:nTimePoints) * 10;
        conditions(iCondition).series(iSeries).alignedData = alignedData;
        conditions(iCondition).series(iSeries).alignedTimes  = (-nTimePoints+1:nTimePoints) * 10;
        
        % Plot aligned data
        figure;
        plotIndividualSeries(alignedTimes, alignedData);
        print(gcf,'-dtiff',fullfile(conditionPath,...
            [conditions(iCondition).name '-' series{iSeries} '-aligneddata.tif']));
        close(gcf)
        
        %% Aggregate data for the condition
        if isempty(pooledData)
            pooledData = alignedData;
            pooledTimes = alignedTimes;
        else
            if size(pooledData, 1) > size(alignedData, 1)
                nRows = (size(pooledData, 1) - size(alignedData, 1))/2;
                addRow = NaN(nRows, size(alignedData,2));
                extendedata =  [addRow; alignedData; addRow];
                pooledData = [pooledData extendedata];
                
            elseif size(pooledData, 1) < size(alignedData, 1)
                nRows = -(size(pooledData, 1) - size(alignedData, 1))/2;
                addRow = NaN(nRows, size(pooledData,2));
                extendedata =  [addRow; pooledData; addRow];
                pooledData = [extendedata alignedData];
                pooledTimes = alignedTimes;
            else
                pooledData = [pooledData alignedData];
            end
        end
    end
    
    % Save aligned times and data
    conditions(iCondition).alignedTimes = pooledTimes;
    conditions(iCondition).alignedData = pooledData;
    
    % Display and save full log
    disp([fileLog{iCondition} seriesLog{:}])
    fid = fopen(fullfile(conditionPath, 'Summary-log.txt'), 'w');
    fprintf(fid, [fileLog{iCondition} seriesLog{:}]);
    fclose(fid);
    
    % Set summary figure options
    conditionFig = figure();
    plotMeanAlignedData(conditions(iCondition).series)
    xlim = get(gca, 'XLim');
    plot([xlim(1) xlim(2)], [minSpindleLength minSpindleLength] ,'--k');
    print(gcf,'-dtiff',fullfile(conditionPath, 'AlignedData.tif'));
    close(conditionFig);
end

% Save output
save(fullfile(dataPath, 'analysis.mat'), 'conditions');

%% Generate conditions graph

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
