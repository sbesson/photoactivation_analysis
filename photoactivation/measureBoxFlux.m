%% Read image conten

iChan = 0;
msg = 'Reading tiles from images %g/%g';
hWaitbar=waitbar(0, sprintf(msg, 1, numel(data)));

for i = 1 : numel(data)
    fprintf(1,'Reading intensity from image %g\n', data(i).id);

    % Read rois and times bounds
    rois = data(i).rois;
    nRois = numel(rois);
    tmin = min(cellfun(@min, {rois.times}));
    tmax = max(cellfun(@max, {rois.times}));
    
    
    for t = tmin : tmax
        for iRoi = 1 : nRois
            roi = rois(iRoi);
            ind = find(roi.times == t);
            if ~isempty(ind)
                % Retrieve rectangle coordinates
                shape = roi.shapes(ind);
                x = shape.getX().getValue();
                y = shape.getY().getValue();
                height = shape.getHeight().getValue();
                width = shape.getWidth().getValue();
                
                % Retrieve tile and read value
                I = getTile(s, data(i).id, 0, iChan, t, x-1, y-1, width+1, height+1);
                data(i).rois(iRoi).values(ind) = mean(I(:));
            end
        end
    end
    waitbar(i/numel(data), hWaitbar, sprintf(msg, i+1, numel(data)));
end
close(hWaitbar)

%%
% define small and large fonts for graphical output
tfont = {'FontName', 'Helvetica', 'FontSize', 14, 'FontAngle', 'italic'};
sfont = {'FontName', 'Helvetica', 'FontSize', 18};
lfont = {'FontName', 'Helvetica', 'FontSize', 22};

%% Read image content

f = figure('PaperPositionMode', 'auto',...
    'Name', data(i).name, 'Position',[50 50 500 500]); % enable resizing

for i = 1:numel(data)
    dI = data(i).rois(1).values - data(i).rois(2).values;
    if dI(1)<0, dI=-dI; end
    
    % plot
    clf
    hold on;
    data(i).dI = dI/dI(1);
    plot(data(i).dI, 'LineWidth',2);
    axis square
    
    % Set thickness of axes, ticks and assign tick labels
    box on
    set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
    xlabel('Time (frames)', lfont{:});
    ylabel('Relative intensity', lfont{:});
    set(gca,'LooseInset',get(gca,'TightInset'))
    
    %     print(gcf, '-dpng', pngfile);
    %     fa = writeFileAnnotation(s, pngfile, 'mimetype', 'image/png');
    %     linkAnnotation(s, fa, 'image', data(i).id);
end
close(f)

%% Group by type

for i = 1:numel(data)
    t = regexp(char(data(i).name),'^(.+).lif','tokens');
    data(i).type = t{1}{1};
end

nTimePoints = cellfun(@numel, {data.dI});
data(nTimePoints<30) = [];
nTimePoints = max(nTimePoints);
t = 1 :nTimePoints;
types = unique({data.type});
nTypes = numel(types);

f = figure('PaperPositionMode', 'auto',...
    'Position',[50 50 500 500]); % enable resizing

h = -ones(nTypes, 1);
colors = hsv(nTypes);
for  i = 1 : numel(types)
    ind = strcmp({data.type},types{i});
    dI = vertcat(data(ind).dI);
    
    hold on;
    if size(dI, 1) > 1
        h(i) = errorbar(t, mean(dI), std(dI, 0, 1), 'LineWidth',2,...
            'Color',colors(i,:));
    else
        h(i) = plot(t, dI, 'LineWidth',2);
    end
    axis square
end

% Set thickness of axes, ticks and assign tick labels
box on
set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
xlabel('Time (frames)', lfont{:});
ylabel('Relative intensity', lfont{:});
xlim([0 nTimePoints])
% legend({data.type})
set(gca,'LooseInset',get(gca,'TightInset'))

path = '/Users/sebastien/Documents/Julie/Photoactivation';
filePath = fullfile(path, 'Comparison.png');
print(gcf, '-dpng', filePath);
fa = writeFileAnnotation(s, filePath, 'mimetype', 'image/png');
linkAnnotation(s, fa, 'dataset', datasetId)

set(get(f, 'Children'), 'XLimMode', 'manual', 'YLimMode', 'manual');
for i = 1:numel(types)
    set(h, 'Visible','off');
    set(h(i), 'Visible','on','Color', 'k');
    filePath = fullfile(path, [types{i} '.png']);
    print(gcf, '-dpng',  fullfile(path, [types{i} '.png']));
    fa = writeFileAnnotation(s, filePath, 'mimetype', 'image/png');
    linkAnnotation(s, fa, 'dataset', datasetId)
end