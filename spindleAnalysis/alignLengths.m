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
        %% Control
        seriesLog{iSeries} = sprintf('\nReading %s...\n',series{iSeries});
        [data, str] = xlsread(xlsfile, series{iSeries});
        if isempty(data), continue; end
        
        data(1,:) = []; % Remove first row
        
        % Align data
        [s, alignmentlog] = alignData(data, minSpindleLength, sigma);
        seriesLog{iSeries} = [seriesLog{iSeries} alignmentlog sprintf('\n')];
        
        % Add aligned data to output structure
        s.name = series{iSeries};
        conditions(iCondition).series(iSeries) = s;
        
        % Plot indidividual events
        figure;
        plotIndividualEvents(conditions(iCondition).series(iSeries));
        print(gcf,'-dtiff',fullfile(conditionPath,...
            [conditions(iCondition).name '-' series{iSeries} '-events.tif']));
        close(gcf)
        
        % Plot aligned data
        figure;
        plotAlignedData(conditions(iCondition).series(iSeries));        
        print(gcf,'-dtiff',fullfile(conditionPath,...
            [conditions(iCondition).name '-' series{iSeries} '-aligneddata.tif']));
        close(gcf)
                
        %% Aggregate data for the condition
        if isempty(pooledData)
            pooledData = s.alignedData;
            pooledTimes = s.alignedTimes;
        else
            if size(pooledData, 1) > size(s.alignedData, 1)
                nRows = (size(pooledData, 1) - size(s.alignedData, 1))/2;
                addRow = NaN(nRows, size(s.alignedData,2));
                extendedata =  [addRow; s.alignedData; addRow];
                pooledData = [pooledData extendedata];
                
            elseif size(pooledData, 1) < size(s.alignedData, 1)
                nRows = -(size(pooledData, 1) - size(s.alignedData, 1))/2;
                addRow = NaN(nRows, size(pooledData,2));
                extendedata =  [addRow; pooledData; addRow];
                pooledData = [extendedata s.alignedData];
                pooledTimes = s.alignedTimes;
            else
                pooledData = [pooledData s.alignedData];
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
