clear
clc
close all

% Define constants
%dataPath = uigetdir('Select the main directory');
dataPath = 'Z:\Sarah\Bipolar spindle assembly assay';
if isequal(dataPath, 0), return; end

% Read conditions from folder names
directories  = dir(dataPath);
directories = directories([directories.isdir]);
conditions = {directories.name};
dotDirs = strcmp(conditions, '.')  | strcmp(conditions, '..');
conditions(dotDirs) = [];
nConditions = numel(conditions);

% Alignment parameters
minSpindleLength = 6; % minimum spindle length (in microns)
sigma = 3; % Sigma to exclude outliers (using final spindle length)

% define small and large fonts for graphical output
tfont = {'FontName', 'Helvetica', 'FontSize', 14, 'FontAngle', 'italic'};
sfont = {'FontName', 'Helvetica', 'FontSize', 18};
lfont = {'FontName', 'Helvetica', 'FontSize', 22};

conditionscolors = hsv(nConditions);
conditionFig = figure();
tmin_cond = Inf;
tmax_cond = -Inf;
    
for iCondition = 1:nConditions
    sheetFig = figure;
    
    % Find Excel file in condition folder
    conditionPath = fullfile(dataPath, conditions{iCondition});
    xlsfiles = dir(fullfile(conditionPath, '*.xls*'));
    assert(numel(xlsfiles) == 1, 'More than one XLS file in folder %s',...
       conditionPath);
    xlsfile = fullfile(conditionPath, xlsfiles(1).name);
    
    %%
    tmin = Inf;
    tmax = -Inf;
    
    fileLog = sprintf('Reading file %s...\n', xlsfile);
    
    % Read number of worksheets
    [status, sheets] = xlsfinfo(xlsfile);
    nSheets = numel(sheets);
    fileLog = [fileLog sprintf('Found %g worskheets\n', nSheets)];
    logs = cell(nSheets, 1);
    sheetcolors = hsv(nSheets);
    pooledData = [];
    
    for iSheet = 1 : nSheets
        %% Control
        logs{iSheet} = sprintf('Reading %s...\n',sheets{iSheet});
        [data, str] = xlsread(xlsfile, sheets{iSheet});
        if isempty(data), continue; end
        
        data(1,:) = []; % Remove first row
        
        % Align data
        [s, alignmentlog] = alignData(data, minSpindleLength, sigma);
        logs{iSheet} = [logs{iSheet} alignmentlog];
        
        % Plot indidividual events
        figure;
        plot(s.times, s.data)
        hold on;
        for i = 1 : s.nSamples
            plot(s.eventTimes(i), s.eventValues(i), 'ok', 'MarkerFaceColor', 'k');
        end
        set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
        xlabel('Time (min)', lfont{:});
        ylabel('Spindle length (\mum)', lfont{:});
        print(gcf,'-dtiff',fullfile(conditionPath,...
            [conditions{iCondition} '-' sheets{iSheet} '-events.tif']));
        close(gcf)
        
        % Plot aligned data
        figure;
        alignedTimes  = (-s.nTimePoints+1:s.nTimePoints) * 10;
        plot(alignedTimes, s.aligneddata)
        
        % Calculate mean and std values of aligned data
        alignedmean = nanmean(s.aligneddata,2);
        alignedstd = nanstd(s.aligneddata, [], 2);
        nPoints  =sum(~isnan(s.aligneddata),2);
        alignedste = alignedstd./sqrt(nPoints);
        
        % Calculate mininimum and maximum times
        imin = max(1,find(~isnan(alignedmean),1,'first'));
        tmin = min(tmin, alignedTimes(imin));
        imax = min(find(~isnan(alignedmean(imin:end)),1,'last')+imin, ...
            2*s.nTimePoints);
        tmax = max(tmax, alignedTimes(imax));
        set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
        xlabel('Time (min)', lfont{:});
        ylabel('Spindle length (\mum)', lfont{:});
        xlim([alignedTimes(imin) alignedTimes(imax)]);
        
        print(gcf,'-dtiff',fullfile(conditionPath,...
            [conditions{iCondition} '-' sheets{iSheet} '-aligneddata.tif']));
        close(gcf)
        
        % Update summary plot
        figure(sheetFig);
        hold on;
        errorbar(alignedTimes, alignedmean, alignedste,...
            '-','Color', sheetcolors(iSheet,:), 'LineWidth', 1);
        
        logs{iSheet} = [logs{iSheet} sprintf('\n\n')];
        
        if isempty(pooledData)
            pooledData = s.aligneddata;
            pooledTimes = alignedTimes;
        else
            if size(pooledData, 1) > size(s.aligneddata, 1)
                nRows = (size(pooledData, 1) - size(s.aligneddata, 1))/2;
                addRow = NaN(nRows, size(s.aligneddata,2));
                extendedata =  [addRow; s.aligneddata; addRow];
                pooledData = [pooledData extendedata];
                
            elseif size(pooledData, 1) < size(s.aligneddata, 1)
                nRows = -(size(pooledData, 1) - size(s.aligneddata, 1))/2;
                addRow = NaN(nRows, size(pooledData,2));
                extendedata =  [addRow; pooledData; addRow];
                pooledData = [extendedata s.aligneddata];
                pooledTimes = alignedTimes;
            else
                pooledData = [pooledData s.aligneddata];
            end
        end
    end
    
    % Display and save full log
    disp([fileLog logs{:}])
    fid = fopen(fullfile(conditionPath, 'Summary-log.txt'), 'w');
    fprintf(fid, [fileLog logs{:}]);
    fclose(fid);
    
    % Set summary figure options
    plot([tmin tmax], [minSpindleLength minSpindleLength] ,'--k');
    box on
    set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
    xlabel('Time (min)', lfont{:});
    ylabel('Spindle length (\mum)', lfont{:});
    xlim([tmin tmax]);
    legend(sheets,'Location', 'NorthWest','Interpreter','None')
    set(gca,'LooseInset',get(gca,'TightInset'));
    
    print(gcf,'-dtiff',fullfile(conditionPath, 'AlignedData.tif'));
    
    % Update summary plot
    figure(conditionFig);

    % Calculate mean and std values of aligned data
    alignedmean = nanmean(pooledData,2);
    alignedstd = nanstd(pooledData, [], 2);
    nPoints  =sum(~isnan(pooledData),2);
    alignedste = alignedstd./sqrt(nPoints);
    
    % Calculate mininimum and maximum times
    imin = max(1,find(~isnan(alignedmean),1,'first'));
    tmin_cond = min(tmin_cond, pooledTimes(imin));
    imax = min(find(~isnan(alignedmean(imin:end)),1,'last'), ...
        2*s.nTimePoints);
    tmax_cond = max(tmax_cond, pooledTimes(imax));
    
    hold on;
    errorbar(pooledTimes', alignedmean, alignedste,...
        '-','Color', conditionscolors(iCondition,:), 'LineWidth', 1);
    
    
end


set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
xlabel('Time (min)', lfont{:});
ylabel('Spindle length (\mum)', lfont{:});
xlim([tmin_cond tmax_cond]);
legend(conditions,'Location', 'NorthWest','Interpreter','None')
set(gca,'LooseInset',get(gca,'TightInset'));
print(gcf,'-dtiff',fullfile(dataPath, 'Comparison.tif'));