function data = analyzeDatasetFlux(session, datasetId, varargin)

ip = inputParser;
ip.addRequired('datasetId', @(x) isscalar(x) && isnumeric(x));
ip.addOptional('invalidIds', [], @(x) isnumeric(x) || isempty(x));
ip.parse(datasetId, varargin{:});
invalidIds = ip.Results.invalidIds;

%%
% Define default output directory
mainOutputDir = fullfile(getenv('HOME'), 'omero');

% Define main namespace
ns = 'photoactivation';

% define small and large fonts for graphical output
tfont = {'FontName', 'Helvetica', 'FontSize', 14, 'FontAngle', 'italic'};
sfont = {'FontName', 'Helvetica', 'FontSize', 18};
lfont = {'FontName', 'Helvetica', 'FontSize', 22};

%Double-exponential function for fitting
fitFun = @(b,x)(b(1) .* exp(b(2) .* x))+(b(3) .* exp(b(4) .* x));     
bInit = [.8 -.1 .2 -0.01];

fitOptions = statset('Robust','on','MaxIter',500,'Display','off');

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
    if ~isdir(outputDir), mkdir(outputDir); end
    
    % Calculating min point (Radon transform
    x0 =  floor((sizeX + 1) / 2);
    y0 =  floor((sizeY + 1) / 2);
    
    % Read rois and times bounds
    points = data(i).points;
    boxes = data(i).boxes;
    nPoints = numel(points);
    nBoxes = numel(boxes);
    tmin = min(cellfun(@min, {points.times}));
    tmax = max(cellfun(@max, {points.times}));
    times = (tmin:tmax)';
    
    
    profilesFig = figure;
    hold on;
    
    % Initializing output
    nTimes = numel(times);
    Ibkg = zeros(nTimes, 1);
    Ibox = zeros(nTimes, nBoxes);
    dmax = zeros(nTimes, 1);
    It = zeros(nTimes, 1);
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
        figure;
        hold on
        plot(x, R, 'Color', color);
        
        % Calculating background
        Ifilt= filterGauss2D(I,5);
        mask = bwareaopen(Ifilt>50, 20);
        threshold = thresholdRosin(Ifilt(mask));
        for iBox = 1: numel(boxes)
            xrange = boxes(iBox).x:boxes(iBox).x + boxes(iBox).width;
            yrange = boxes(iBox).y:boxes(iBox).y + boxes(iBox).height;
            Icrop = I(yrange, xrange);
            Ibox(iT, iBox) = mean(Icrop(:));
        end
%         Ibkg(iT) = min(Ibox(iT, :));
%         mean(I(mask & I < threshold));
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
        hold on
        d = x-P(iT, 1);
        plot(d, dR, 'Color', color);
        xlim([0 max(d_centrosomes)])
        
        % Integrate intensity
        [~, imax] = max(dR);
        xrange = imax-20:imax+20;
        
        figure;
        [p, dp] = fitGaussian1D(d(xrange), dR(xrange), [50 max(dR(xrange)) 20 0],'xAs');
        plot(d, dR,'ok');
        yfit =  exp(-(d(xrange)-p(1)).^2./(2*p(3)^2))*p(2) + p(4);
        hold on
        plot(d(xrange), yfit,'-k');
        print(gcf, '-dpng', fullfile(outputDir, ['Time' num2str(t) '-fitted.png']));
        close(gcf);
        
        % Integrate intensity
        It(iT) = sqrt(2*pi) * p(3) *p(2);
        dmax(iT) =p(1);
        plot(dmax(iT), p(2) + p(4), 'o', 'Color', color, 'MarkerFaceColor', color);
    end
    
    saveAndUploadEPS(profilesFig, outputDir, 'Profiles', session, data(i).id, [ns '.profiles'])
    
    % Plot band position
    fluxFig = figure;
    plot(times * dT, (dmax-dmax(1)) * pixelSize, 'ok');
    coeff = polyfit((times + 1) * dT, (dmax-dmax(1)) * pixelSize, 1);
    hold on
    plot(times * dT, coeff(1) * (times + 1) * dT + coeff(2), '--k');
    box on
    set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
    xlabel('Time (s)', lfont{:});
    ylabel('Relative position (microns)', lfont{:});
    
    % Save flux results locally
    saveAndUploadEPS(fluxFig, outputDir, 'Flux', session, data(i).id, [ns '.flux'])
    
    % Print flux data
    fprintf(1, 'Band speed: %g microns/min\n', abs(coeff(1)) * 60);
    data(i).speed = abs(coeff(1)) *60;

    %Fit function to ratio timeseries
    [bFit,resFit,~,covFit,mseFit] = nlinfit(times*dT,It/It(1),fitFun,bInit,fitOptions);
    %Get confidence intervals of fit and fit values
    [fitValues,deltaFit] = nlpredci(fitFun,times*dT,bFit,resFit,'covar',covFit,'mse',mseFit);
    
    turnoverFig = figure;
    plot(times * dT, It/It(1), 'ok');
    hold on
    plot(times * dT, fitValues, '--k');
    box on
    set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
    xlabel('Time (s)', lfont{:});
    ylabel('Integrated intensity', lfont{:});
    limits = ylim();
    ylim([0 limits(2)]);
    
    % Save turnover results locally
    fprintf(1, 'Turnover time: %g s\n', -1/bFit(2));
    data(i).turnover_time = -1/bFit(2);
     
    % Save turnover results onto the server
    saveAndUploadEPS(turnoverFig, outputDir, 'Turnover', session, data(i).id, [ns '.turnover'])
    
    % Close reader if using local image
    if hasCacheImage, r.close(); end
end

%%

outputDir = fullfile(getenv('HOME'), 'omero', num2str(datasetId));
if ~isdir(outputDir), mkdir(outputDir); end

% Create results table
resultsPath = fullfile(outputDir, ['Photoactivation_results_'...
    num2str(datasetId) '.txt']);
fid = fopen(resultsPath ,'w+');
fprintf(fid, 'Name\tId\tSpeed (microns/min)\tTurnover time (s)\n');
for i = 1:numel(data),
    fprintf(fid, '%s\t%g\t%g\t%g\n', data(i).name, data(i).id, data(i).speed,...
        data(i).turnover_time);
end
fclose(fid);

% Upload results file to OMERO
results_ns = [ns '.results'];
fa = getDatasetFileAnnotations(session, datasetId, 'include', results_ns);
if isempty(fa),
    fa = writeFileAnnotation(session, resultsPath, 'namespace', results_ns);
    linkAnnotation(session, fa, 'dataset', datasetId);
else
    updateFileAnnotation(session, fa, resultsPath);
end