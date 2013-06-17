%%

cacheDir = '/Users/sebastien/Documents/Julie/Photoactivation';
load(fullfile(cacheDir, 'data.mat'));

hasPoints = ~cellfun(@isempty, {data.points});
data(~hasPoints) = [];


%%
% define small and large fonts for graphical output
tfont = {'FontName', 'Helvetica', 'FontSize', 14, 'FontAngle', 'italic'};
sfont = {'FontName', 'Helvetica', 'FontSize', 18};
lfont = {'FontName', 'Helvetica', 'FontSize', 22};


%% Read image conten

iChan = 0;
msg = 'Reading images %g/%g';
hWaitbar=waitbar(0, sprintf(msg, 1, numel(data)));
close all

for i = 1:numel(data)
    fprintf(1,'Reading intensity from image %g:%s\n', data(i).id,...
        data(i).name);
    
    cacheImage = fullfile(cacheDir, [num2str(data(i).id) '.ome.tiff']);
    hasCacheImage = exist(cacheImage, 'file') == 2;
    if hasCacheImage
        fprintf(1,'using cache image at %s\n', cacheImage);
        r = bfGetReader(cacheImage);
        
        % Calculating min point (Radon transform
        x0 =  floor((r.getSizeX() + 1) / 2);
        y0 =  floor((r.getSizeY() + 1) / 2);
        
        % Read metadata
        store = r.getMetadataStore();
        pixelSize =  store.getPixelsPhysicalSizeX(0).getValue();
        dT = double(store.getPixelsTimeIncrement(0));
        
        % Create output directory
        outputDir = fullfile(cacheDir, num2str(data(i).id));
        if ~isdir(outputDir), mkdir(outputDir); end
    else
        disp('using OMERO.matlab');
    end
    
    % Read rois and times bounds
    points = data(i).points;
    boxes = data(i).boxes;
    nPoints = numel(points);
    nBoxes = numel(boxes);
    tmin = min(cellfun(@min, {points.times}));
    tmax = max(cellfun(@max, {points.times}));
    times = (tmin:tmax)';
    
    
    f1=figure;
    hold on;
    
    % Initializing output
    nTimes = numel(times);
    A = zeros(nTimes, 1);
    dA = zeros(nTimes, 1);
    Ibkg = zeros(nTimes, 1);
    Ibox = zeros(nTimes, nBoxes);
    dmax = zeros(nTimes, 1);
    dRmax = zeros(nTimes, 1);
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
            I = double(getPlane(s, data(i).id, 0, iChan, t-1));
        end
        
        
        % Calculating Radon transform
        fprintf(1, 'alpha: %g\n', alpha(iT));
        [R,x] = radon(I, -alpha(iT));
        figure;
        hold on
        plot(x, R, 'Color', color);
        
        % Calculating background
        mask = bwareaopen(filterGauss2D(I,1)>50,20);
        threshold = thresholdRosin(I(mask));
        for iBox = 1: numel(boxes)
            xrange = boxes(iBox).x:boxes(iBox).x + boxes(iBox).width;
            yrange = boxes(iBox).y:boxes(iBox).y + boxes(iBox).height;
            Icrop = I(yrange, xrange);
            Ibox(iT, iBox) = mean(Icrop(:));
        end
        Ibkg(iT) = min(Ibox(iT, :));
        mean(I(mask & I < threshold));
        Ibkg(iT) = mean(I(mask & I < threshold));
        I0=I;
        I0(mask) = Ibkg(iT);
        [Rbg, xbg] = radon(I0,-alpha(iT));
        plot(xbg, Rbg, '--', 'Color', color);
        print(gcf, '-dpng', fullfile(outputDir, ['Time' num2str(t) '.png']));
        close(gcf)
        
        % Fitting residuals
        dR = R - Rbg;
        figure(f1);
        hold on
        d = x-P(iT, 1);
        plot(d, dR, 'Color', color);
        xlim([0 max(d_centrosomes)])
        
        % Interpolate maximum
        %         [~, imax] = max(dR);
        %         coeff=polyfit(d(imax-7:imax+7), dR(imax-7:imax+7),2);
        %         dmax(iT) = - coeff(2)/(2*coeff(1));
        %         dRmax(iT) =  coeff(1) * dmax(iT)^2  + coeff(2) *dmax(iT) + coeff(3);
        %         plot(dmax(iT), dRmax(iT), 'o', 'Color', color, 'MarkerFaceColor', color);
        %         [p, dp] = fitGaussian1D(x, dR, [-50 10e4 20 1]);
        %         A(iT)=p(2);
        %         dA(iT)=dp(2);
        
        % Integrate intensity
        figure;
        xrange = (x>=P(iT, 1) & x<=2*P(iT, 2)/3);
        [p, dp] = fitGaussian1D(d(xrange), dR(xrange), [50 max(dR(xrange)) 20 1],'xAsc');
        plot(d(xrange), dR(xrange),'ok');
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
    
    print(f1, '-dpng', fullfile(outputDir, 'Projections.png'));
    
    figure;
    plot(times * dT, dmax * pixelSize, 'o');
    coeff = polyfit((times + 1) * dT, dmax * pixelSize, 1);
    hold on
    plot(times * dT, coeff(1) * (times + 1) * dT + coeff(2), '--');
    box on
    set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
    xlabel('Time (s)', lfont{:});
    ylabel('Relative position (microns)', lfont{:});
    print(gcf, '-dpng', fullfile(outputDir, 'BandPosition.png'));
    fprintf(1, 'Band speed: %g microns/min\n', abs(coeff(1)) * 60);
    
    fitFun = @(b,x)(b(1) .* exp(b(2) .* x))+(b(3) .* exp(b(4) .* x));     %Double-exponential function for fitting
    bInit = [.8 -.1 .2 -0.01];
    %Fit function to ratio timeseries
    fitOptions = statset('Robust','on','MaxIter',500,'Display','off');
    [bFit,resFit,jacFit,covFit,mseFit] = nlinfit(times*dT,It/It(1),fitFun,bInit,fitOptions);
    %Get confidence intervals of fit and fit values
    [fitValues,deltaFit] = nlpredci(fitFun,times*dT,bFit,resFit,'covar',covFit,'mse',mseFit);
    
    figure;
    plot(times * dT, It/It(1), 'ok');
    hold on
    plot(times * dT, fitValues, '--k');
    box on
    set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
    xlabel('Time (s)', lfont{:});
    ylabel('Integrated intensity', lfont{:});
    print(gcf, '-dpng', fullfile(outputDir, 'Intensity.png'));
    fprintf(1, 'Turnover time: %g s\n', -1/bFit(2));
    
    if hasCacheImage, r.close(); end
    %     waitbar(i/numel(data), hWaitbar, sprintf(msg, i+1, numel(data)));
end
% close(hWaitbar)
%%

figure;
plot(A/A(1));
hold on
plot(fitValues,'--r');