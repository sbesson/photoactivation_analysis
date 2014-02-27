function data = combineDatasetPAAnalysis(session, datasetId)


% Define main namespace
ns = 'photoactivation';

mainOutputDir = fullfile(getenv('HOME'), 'omero');
outputDir = fullfile(mainOutputDir, num2str(datasetId));
s = load(fullfile(outputDir, 'analysis.mat'));
data = s.data;
validData = arrayfun(@(x) ~isempty(x.position), data);

% Save positions table
fprintf(1, 'Writing positions for dataset %g\n', datasetId);
positionsPath = fullfile(outputDir, ['Positions_' num2str(datasetId) '.txt']);
fid = fopen(positionsPath ,'w');

fprintf(fid, 'Times\t');
for i = find(validData)',
    fprintf(fid, '\t%s', data(i).name);
end
fprintf(fid, '\n');
for t = 1: numel(data(1).times),
    fprintf(fid, '%g', data(1).times(t));
    for j = find(validData)',
        fprintf(fid, '\t%g', data(j).position(t));
    end
    fprintf(fid, '\n');
end
fclose(fid);

% Create results table
uploadFileResults(session, positionsPath, 'dataset', datasetId, [ns '.positions']);


% Save intensity table
fprintf(1, 'Writing intensities for dataset %g\n', datasetId);
intensitiesPath = fullfile(outputDir, ['Intensities_' num2str(datasetId) '.txt']);
fid = fopen(intensitiesPath ,'w');

fprintf(fid, 'Times\t');
for i = find(validData)',
    fprintf(fid, '\t%s', data(i).name);
end
fprintf(fid, '\n');
for t = 1: numel(data(1).times),
    fprintf(fid, '%g', data(1).times(t));
    for j = find(validData)',
        fprintf(fid, '\t%g', data(j).dInorm(t));
    end
    fprintf(fid, '\n');
end
fclose(fid);

% Create results table
uploadFileResults(session, intensitiesPath, 'dataset', datasetId, [ns '.intensities']);

% Calculate speed for the mean positions
all_positions = horzcat(data.position);
negative_slopes = all_positions(end, :) <0;
all_positions(:, negative_slopes) = - all_positions(:, negative_slopes);


% Create filters
longtimes_filter = [data.t2]  < 1e10;
r2_filter_0 = [data.r2] >= 0;
r2_filter_90 = [data.r2] >= .9;
r2_filter_95 = [data.r2] >= .95;
filters(1).name = 'Unfiltered';
filters(1).values = ones(size([data.t2]));
filters(2).name = 'Long times';
filters(2).values = longtimes_filter & r2_filter_0;
filters(3).name = 'Long times & .90';
filters(3).values = longtimes_filter & r2_filter_90;
filters(4).name = 'Long times & .95 ';
filters(4).values = longtimes_filter & r2_filter_95;
filters(5).name = '.90 ';
filters(5).values = r2_filter_90;
filters(6).name = '.95 ';
filters(6).values = r2_filter_95;

% Average and fit the normalized intensities
dataTypes(1).name = 'uncorrected';
dataTypes(1).field = 'dInorm';
if isfield(data, 'dInorm_corr')
    dataTypes(2).name = 'corrected';
    dataTypes(2).field = 'dInorm_corr';
end
maxTimes = [200, 300];
logmsg = 'Writing results for %s data of dataset %g with a cutoff of %gs\n';
for iType = 1 : numel(dataTypes)
    for maxTime = maxTimes

        % Create results table
        fprintf(1, logmsg, dataTypes(iType).name, datasetId, maxTime);
        resultsPath = fullfile(outputDir, ['PA_' dataTypes(iType).name...
             '_' num2str(maxTime) 's_' num2str(datasetId) '.txt']);
        fid = fopen(resultsPath ,'w');
        fprintf(fid, ['Name\tId\tFull-range speed (microns/min)\t'...
            'Half-range speed (microns/min)\t'...
            'Fast turnover time (s)\tSlow turnover time (s)\tr2\n']);
        for i = find(validData)',
            [bFit, r2] = fitPADecay(data(i).times,...
                data(i).(dataTypes(iType).field), maxTime);
            fprintf(fid, '%s\t%g\t%g\t%g\t%g\t%g\t%g\n', data(i).name, data(i).id,...
                data(i).speed_full, data(i).speed_half,...
                -1/bFit(2), -1/bFit(4), r2);
        end
        
        allInorm = horzcat(data.(dataTypes(iType).field));
        speed_full =  zeros(numel(filters), 1);
        speed_half =  zeros(numel(filters), 1);
        t1_filt = zeros(numel(filters), 1);
        t2_filt = zeros(numel(filters), 1);
        r2_filt = zeros(numel(filters), 1);
        for iFilter = 1: numel(filters)
            % Calculate speed for the mean positions
            all_positions_filt = mean(all_positions(:, filters(iFilter).values), 2);
            coeff = polyfit(data(1).times, all_positions_filt, 1);
            speed_full(iFilter) = abs(coeff(1)) *60;
            coeff = polyfit(data(1).times(1:end/2), all_positions_filt(1:end/2), 1);
            speed_half(iFilter) = abs(coeff(1)) *60;
            
            dInorm_filt = mean(allInorm(:, filters(iFilter).values), 2);
            [bFit, r2] = fitPADecay(data(1).times, dInorm_filt, maxTime);
            t1_filt(iFilter) = -1/bFit(2);
            t2_filt(iFilter) = -1/bFit(4);
            r2_filt(iFilter) = r2;
        end
        
        
        fprintf(fid, '\n');
        fprintf(fid, ['Filter\tnCells\tFull-range speed (microns/min)\t'...
            'Half-range speed (microns/min)\t'...
            'Fast turnover time (s)\tSlow turnover time (s)\tr2\n']);
        for i = 1 : numel(filters)
            fprintf(fid, '%s\t%s\t%g\t%g\t%g\t%g\t%g\n', filters(i).name,...
                sum(filters(i).values), speed_full(i), speed_half(i),...
                t1_filt(i), t2_filt(i), r2_filt(i));
        end
        fclose(fid);
        
        % Upload results file to OMERO
        uploadFileResults(session, resultsPath, 'dataset', datasetId,...
            [ns '.result.' dataTypes(iType).name '.' num2str(maxTime)]);
    end
end

