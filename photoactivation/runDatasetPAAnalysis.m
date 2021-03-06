function data = runDatasetPAAnalysis(session, datasetId, varargin)


ip = inputParser;
ip.addRequired('datasetId', @(x) isscalar(x) && isnumeric(x));
ip.addOptional('corrData', [], @(x) isnumeric(x) || isempty(x));
ip.addParamValue('invalidIds', [], @(x) isnumeric(x) || isempty(x));
ip.parse(datasetId, varargin{:});
corrData = ip.Results.corrData;
invalidIds = ip.Results.invalidIds;

%%
% Define default output directory
mainOutputDir = fullfile(getenv('HOME'), 'omero');

% Define main namespace
ns = 'photoactivation';

% define small and large fonts for graphical output
sfont = {'FontName', 'Helvetica', 'FontSize', 18};
lfont = {'FontName', 'Helvetica', 'FontSize', 22};

%% Read images

% Load dataset
fprintf(1,'Loading dataset %g\n', datasetId);
dataset = getDatasets(session, datasetId);
datasetname = char(dataset.getName().getValue());

% Retrieve list of images
images = toMatlabList(dataset.linkedImageList);
imageIds = arrayfun(@(x) x.getId().getValue(), images);
fprintf(1,'Found %g images in dataset %s\n', numel(imageIds), datasetname);

%% Read rois
roiService = session.getRoiService();
nImages = numel(imageIds);

clear data
data(nImages, 1) = struct('id',[],'boxes', [], 'points', []);
for iImage = 1 : nImages
    % Read image properties
    data(iImage).id = images(iImage).getId().getValue();
    data(iImage).name = char(images(iImage).getName().getValue());
    fprintf(1,'Reading image %g\n', data(iImage).id);
    
    % Read ROIs attached to the image
    roiResult = roiService.findByImage(data(iImage).id, []);
    if isempty(roiResult.rois), continue; end
    
    allRois = toMatlabList(roiResult.rois);
    nRois = numel(allRois);
    
    if nRois < 2, continue; end
    
    for iRoi = 1 : nRois
        % Only store ROIs containing rectangles
        shapes = allRois(iRoi).copyShapes;
        if isa(shapes.get(0), 'omero.model.RectI')
            iBox = numel(data(iImage).boxes)+1;
            shapes = toMatlabList(shapes);
            data(iImage).boxes(iBox).times = arrayfun(@(x) x.getTheT.getValue, shapes);
            data(iImage).boxes(iBox).x = arrayfun(@(x) x.getX().getValue, shapes);
            data(iImage).boxes(iBox).y = arrayfun(@(x) x.getY().getValue, shapes);
            data(iImage).boxes(iBox).width = arrayfun(@(x) x.getWidth().getValue, shapes);
            data(iImage).boxes(iBox).height = arrayfun(@(x) x.getHeight().getValue, shapes);
        elseif isa(shapes.get(0), 'omero.model.PointI')
            iPoint = numel(data(iImage).points)+1;
            shapes = toMatlabList(shapes);
            data(iImage).points(iPoint).times = arrayfun(@(x) x.getTheT.getValue, shapes);
            data(iImage).points(iPoint).x = arrayfun(@(x) x.getCx().getValue, shapes);
            data(iImage).points(iPoint).y = arrayfun(@(x) x.getCy().getValue, shapes);
        end
    end
    
    % Filter valid rois with no shapes
    fprintf(1,'Found %g rectangle ROIs and %g point ROIs \n',...
        numel(data(iImage).boxes),numel(data(iImage).points));
    
end

% Filter data with no ROI
% hasBoxes = ~cellfun(@isempty, {data.boxes});
% hasPoints = ~cellfun(@isempty, {data.points});
has2Centrosomes = cellfun(@numel, {data.points}) == 2;
data = data(has2Centrosomes);
nMeasurements = arrayfun(@(x) numel(x.points(1).x), data);
data = data(nMeasurements > 5);

%% Read image conten

iChan = 0;
close all
vmax = 2; % Maximum velocity of the primiary signal
nSigma = 2;

for i = 1:numel(data)
    if ismember(data(i).id, invalidIds), continue; end
    fprintf(1, 'Reading intensity from image %g:%s\n', data(i).id,...
        data(i).name);
    
    cacheImage = fullfile(mainOutputDir, [num2str(data(i).id) '.ome.tiff']);
    hasCacheImage = exist(cacheImage, 'file') == 2;
    if hasCacheImage
        fprintf(1,'using cache image at %s\n', cacheImage);
        r = bfGetReader(cacheImage);
        sizeX = r.getSizeX();
        sizeY = r.getSizeY();
        
        % Read metadata
        store = r.getMetadataStore();
        pixelSize =  store.getPixelsPhysicalSizeX(0).getValue();
        dT = double(store.getPixelsTimeIncrement(0));
        
    else
        image = getImages(session, data(i).id);
        
        sizeX = image.getPrimaryPixels.getSizeX.getValue;
        sizeY = image.getPrimaryPixels.getSizeY.getValue;
        pixelSize =  image.getPrimaryPixels.getPhysicalSizeX().getValue();
        dT = image.getPrimaryPixels.getTimeIncrement.getValue();
        
        disp('using OMERO.matlab');
    end
    
    % Define output directory
    outputDir = fullfile(mainOutputDir, num2str(data(i).id));
    if isdir(outputDir), rmdir(outputDir, 's'); end
    mkdir(outputDir);
    
    % Calculating min point (Radon transform
    x0 =  floor((sizeX + 1) / 2);
    y0 =  floor((sizeY + 1) / 2);
    
    % Read rois and times bounds
    points = data(i).points;
    nPoints = numel(points);
    tmin = min(cellfun(@min, {points.times}));
    tmax = max(cellfun(@max, {points.times}));
    times = (tmin:tmax)';
    
    
    profilesFig = figure('Visible','off');
    hold on;
    
    % Initializing output
    nTimes = numel(times);
    Ibkg = zeros(nTimes, 1);
    Isignal = zeros(nTimes, 1);
    Asignal = zeros(nTimes, 1);
    dsignal = zeros(nTimes, 1);
    Isignal2 = zeros(nTimes, 1);
    Asignal2 = zeros(nTimes, 1);
    dsignal2 = zeros(nTimes, 1);
    
    colors = hsv(nTimes+1);
    
    % Calculate angles
    dx = diff(horzcat(points.x),[],2);
    dy = diff(horzcat(points.y),[],2);
    alpha = atan2(dy, dx) * 180/pi;
    
    dx0 = (horzcat(points.x)-x0);
    dy0 = (horzcat(points.y)-y0);
    alpha0 = atan2(dy0, dx0) * 180/pi;
    dL = (dx0.^2 + dy0.^2).^(.5);
    P = dL .* cos((repmat(alpha,1,nPoints) - alpha0)*pi/180);
    d_centrosomes = diff(P, [], 2);
    
    heat_map = NaN(numel(times), round(max(d_centrosomes)));
    for iT = 1 : numel(times)
        t = times(iT);
        color = colors(iT, :);
        
        % Retrieve tile and read value
        fprintf(1,'Reading plane %g...', t);
        if hasCacheImage
            iPlane = loci.formats.FormatTools.getIndex(r, 0, iChan, t);
            I = double(bfGetPlane(r, iPlane + 1));
        else
            I = double(getPlane(session, data(i).id, 0, iChan, t));
        end
        
        
        % Calculating Radon transform
        fprintf(1, 'alpha: %g\n', alpha(iT));
        [R,x] = radon(I, -alpha(iT));
        figure('Visible','off');
        hold on
        plot(x, R, 'Color', color);
        
        % Calculating background
        Ifilt= filterGauss2D(I,5);
        mask = bwareaopen(Ifilt>50, 20);
        threshold = thresholdRosin(Ifilt(mask));
        Ibkg(iT) = mean(I(mask & Ifilt < threshold));
        I0=I;
        I0(mask) = Ibkg(iT);
        [Rbg, xbg] = radon(I0,-alpha(iT));
        plot(xbg, Rbg, '--', 'Color', color);
        print(gcf, '-dpng', fullfile(outputDir, ['Time' num2str(t) '.png']));
        close(gcf)
        
        % Fitting residuals
        dR = R - Rbg;
        figure(profilesFig);
        set(gcf, 'Visible', 'off');
        hold on
        d = x-P(iT, 1);
        plot(d, dR, 'Color', color);
        xlim([0 max(d_centrosomes)])
        i0 = find(round(d)==0);
        dmax = round(max(d_centrosomes));
        heat_map(iT, (1:dmax)) = R(i0:(i0 + dmax -1));

        % Integrate intensity
        if iT == 1
            % Look for absolute maximum on the first frame
            [~, imax] = max(dR);
        else
            % Look for absolute maximum in the neighborood of previous
            % maximumu
            delta = vmax * dT / 60 / pixelSize;
            localrange = find(d > dsignal(iT-1)-delta &...
                d < dsignal(iT-1) + delta);
            [~, imax] = max(dR(localrange));
            imax = localrange(imax);
        end
        
        % Fit signal with a 1 D Gaussian
        signalrange = imax-20:imax+20;
        p = fitGaussian1D(d(signalrange), dR(signalrange),...
            [50 max(dR(signalrange)) 20 0],'xAs');
        % Isignal(iT) = sqrt(2*pi) * p(3) *p(2); % Integrated intensity
        Asignal(iT) = p(2); % Amplitude
        norm_int = @(x) erf(x/(p(3)*sqrt(2)));
        Isignal(iT) = p(2) * (norm_int(p(3)) - norm_int(-p(3)))/2;
        dsignal(iT) = p(1); %position
        signal_fit =  exp(-(d-p(1)).^2./(2*p(3)^2))*p(2) + p(4);
        fprintf(1, 'Primary signal detected at position %g with intensity %g\n',...
            dsignal(iT), Isignal(iT));
        
        % Calculate residuals to fit other spindle signal
        residuals = dR - signal_fit;
        signalrange = d > p(1) - nSigma * p(3) & d < p(1) + nSigma*p(3);
        if p(1) > d_centrosomes(iT)/2
            fullrange = d > 0 & d < d_centrosomes(iT)/2;
        else
            fullrange = d > d_centrosomes(iT)/2 & d < d_centrosomes(iT);
        end
        signal2range = fullrange & ~signalrange;
        
        % Fit residuals with a 1D Gaussian
        if ~any(signal2range) || max(residuals(signal2range)) < 0.05 * Isignal(iT),
            disp('No secondary signal  detected');
            signal2_fit = zeros(size(d));
            Asignal2(iT) = 0;
            Isignal2(iT) = 0; % amplitude
            dsignal2(iT) = 0;
        else
            try
                p2 = fitGaussian1D(d(signal2range), residuals(signal2range),...
                    [mean(d(signal2range)) max(residuals(signal2range)) 20 0],'xAs');
                signal2_fit =  exp(-(d-p2(1)).^2./(2*p2(3)^2))*p2(2) + p2(4);
                % Isignal2(iT) = sqrt(2*pi) * p2(3) *p2(2); % integrated
                Isignal2(iT) = p2(2) * (norm_int(p2(3)) - norm_int(-p2(3)))/2;
                Asignal2(iT) = p2(2); % amplitude
                dsignal2(iT) = p2(1);
                fprintf(1, 'Secondary signal maximum detected at position %g with intensity %g\n',...
                    dsignal2(iT), Isignal2(iT));
            catch ME
                disp(ME.message);
                signal2_fit = zeros(size(d));
                Asignal2(iT) = 0;
                Isignal2(iT) = 0; % amplitude
                dsignal2(iT) = 0;
            end
        end
        
        % Plot fitted signals
        figure('Visible','off');
        hold on
        plot(d, dR,'ok');
        plot(d, signal_fit,'-k');
        plot(d, signal2_fit,'--k');
        xlim([0 max(d_centrosomes)])
        print(gcf, '-dpng', fullfile(outputDir, ['Time' num2str(t) '-fitted.png']));
        close(gcf);
    end
    
    % Tranform times into real units
    times = times * dT;
    data(i).times = times;
    
    % Save profiles locally and on the server
    profilesPath = fullfile(outputDir, ['Profiles_' num2str(data(i).id) '.eps']);
    print(profilesFig, '-depsc', profilesPath);
    close(profilesFig)
    uploadFileResults(session, profilesPath, 'image',  data(i).id, [ns '.profiles']);
    
    % Plot heat map
    heatMapFig = figure('Visible','off');
    heatMapPath = fullfile(outputDir, ['HeatMap_' num2str(data(i).id) '.eps']);
    imagesc(heat_map);
    print(heatMapFig, '-depsc', heatMapPath);
    close(heatMapFig)
    uploadFileResults(session, heatMapPath, 'image', data(i).id, [ns '.heatmap'])
    
    % Plot band position
    fluxFig = figure('Visible','off');
    plot(times, (dsignal-dsignal(1)) * pixelSize, 'ok');
    coeff_full = polyfit(times, (dsignal - dsignal(1)) * pixelSize, 1);
    coeff_half = polyfit(times(1:end/2),...
        (dsignal(1:end/2)-dsignal(1)) * pixelSize, 1);
    hold on
    plot(times(1:end/2), coeff_half(1) * times(1:end/2) + coeff_half(2), '-.k');
    plot(times, coeff_full(1) * times + coeff_full(2), '--k');
    box on
    set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
    xlabel('Time (s)', lfont{:});
    ylabel('Relative position (microns)', lfont{:});
    
    % Save flux results locally and on the server
    fluxPath = fullfile(outputDir, ['Flux_' num2str(data(i).id) '.eps']);
    print(fluxFig, '-depsc', fluxPath);
    close(fluxFig)
    uploadFileResults(session, fluxPath, 'image', data(i).id, [ns '.flux']);
    
    % Print flux data
    data(i).position = (dsignal-dsignal(1)) * pixelSize;
    data(i).speed_half = abs(coeff_half(1)) * 60;
    data(i).speed_full = abs(coeff_full(1)) * 60;
    fprintf(1, 'Band speed (half range): %g microns/min\n', data(i).speed_half);
    fprintf(1, 'Band speed (full range): %g microns/min\n', data(i).speed_full);
    
    %Fit function to ratio timeseries
    dI = Asignal - Asignal2;
    dInorm = dI/dI(1);
    
    [bFit, r2, yfit] = fitPADecay(times, dInorm);
    data(i).r2 = r2;
    
    % Save turnover results locally
    data(i).dInorm = dInorm;
    data(i).t1 = -1/bFit(2);
    data(i).t2 = -1/bFit(4);
    fprintf(1, 'Fast turnover time: %g s\n', data(i).t1);
    fprintf(1, 'Slow turnover time: %g s\n', data(i).t2);
    
    
    turnoverFig = figure('Visible','off');
    plot(times, dInorm, 'ok');
    hold on
    plot(times, yfit, '--k');
    box on
    set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
    xlabel('Time (s)', lfont{:});
    ylabel('Integrated intensity', lfont{:});
    limits = ylim();
    ylim([0 limits(2)]);
    
    % Compute corr
    if ~isempty(corrData);
        dInorm_corr = dInorm./corrData;
        [bFit_corr, r2_corr, yfit_corr] = fitPADecay(times, dInorm_corr);
        data(i).dInorm_corr = dInorm_corr;
        data(i).r2_corr = r2_corr;
        data(i).t1_corr = -1/bFit_corr(2);
        data(i).t2_corr = -1/bFit_corr(4);
        plot(times, dInorm_corr, 'or');
        plot(times, yfit_corr, '--r');
    end
    
    % Save flux results locally and on the server
    turnoverPath = fullfile(outputDir, ['Turnover_' num2str(data(i).id) '.eps']);
    print(turnoverFig, '-depsc', turnoverPath);
    close(turnoverFig)
    uploadFileResults(session, turnoverPath, 'image', data(i).id, [ns '.turnover']);
    
    % Close reader if using local image
    if hasCacheImage, r.close(); end
end

%% Dataset results

outputDir = fullfile(mainOutputDir, num2str(datasetId));
if isdir(outputDir), rmdir(outputDir, 's'); end
mkdir(outputDir);

save(fullfile(outputDir, 'analysis.mat'), 'data');