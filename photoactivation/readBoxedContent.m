%% OMERO.matlab initialization

% Create client/ession
c = loadOmero('~/Documents/MATLAB/sce_ice.config');
c.enableKeepAlive(60);
fprintf(1, 'Connecting to %s\n', char(c.getProperty('omero.host')));
s = c.createSession();

% Change group
disp('Changing group')
group = s.getAdminService.lookupGroup('Rape project');
s.setSecurityContext(group);

%% Read images

% Load dataset
datasetId = 902;
fprintf(1,'Loading dataset %g\n', datasetId);
dataset = getDatasets(s, datasetId);

% Retrieve list of images
images = toMatlabList(dataset.linkedImageList);
imageIds = arrayfun(@(x) x.getId.getValue, images);
fprintf(1,'Found %g images in dataset %s\n', numel(imageIds),...
    char(dataset.getName.getValue));

%% Read rois
roiService = s.getRoiService();
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
hasBoxes = ~cellfun(@isempty, {data.boxes});
hasPoints = ~cellfun(@isempty, {data.points});
data = data(hasBoxes | hasPoints);

