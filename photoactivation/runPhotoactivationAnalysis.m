% Configuration

% OMERO.matlab configuration file
configFile = fullfile(fileparts(which('startup.m')), 'ice.config');
if ~exist(configFile, 'file') == 2,
    [filename, path] = uigetfile('Select an OMERO configuration file');
    assert(~isequal(filename, 0), 'No file selected');
    configFile = fullfile(path, filename);
end
fprintf(1, 'Using OMERO configuration file located in %s\n', configFile);

% Project to use
defaultgroup = 'Rape project';
groupQuery = 'Specify the name of the group [default: %s]: ';
groupName = input(sprintf(groupQuery, defaultgroup), 's');
if isempty(groupName), groupName = defaultgroup; end

% Define id of projects or datasets to be analyzed
ids = input('Specify the ids of the projects or datasets to analyze:');
assert(~isempty(ids), 'No ids selected');

% Define id of projects or correction datasets to be used
defaultCorrId = 1102;
corrQuery = 'Specify the ids of the photobleaching dataset [default: %g]: ';
corrId = input(sprintf(corrQuery, defaultCorrId));
if isempty(corrId), corrId = defaultCorrId; end
assert(isscalar(corrId), 'Only one photobleaching Id allowed');

if ~isempty(corrId),
    defaultCorrGroup = 'Rape project';
    groupQuery = 'Specify the name of the correction group [default: %s]: ';
    corrGroupName = input(sprintf(groupQuery, defaultCorrGroup), 's');
    if isempty(corrGroupName), corrGroupName = defaultCorrGroup; end
end
    
% Wrong images (to exclude from the analysis)
invalidIds = [3991, 3412, corrId];
%% OMERO.matlab initialization

% Create client/session
[client, session] = loadOmero(configFile);
fprintf(1, 'Created connection to %s\n', char(client.getProperty('omero.host')));
fprintf(1, 'Created session for user %s using group %s\n',...
    char(session.getAdminService().getEventContext().userName),...
    char(session.getAdminService().getEventContext().groupName));


%%
if ~isempty(corrId)
    % Change group
    disp('Changing group')
    corrGroup = session.getAdminService.lookupGroup(corrGroupName);
    session.setSecurityContext(corrGroup);
    
    % Read file annotation
    fa = getDatasetFileAnnotations(session, corrId,...
        'include', 'photoactivation.intensities');
    outputDir = fullfile(getenv('HOME'), 'omero', num2str(corrId));
    if ~isdir(outputDir), mkdir(outputDir); end
    outputFile = fullfile(outputDir, char(fa.getFile().getName().getValue()));
    getFileAnnotationContent(session, fa, outputFile);
    
    % Read correction data
    corrData = dlmread(outputFile, '\t', 1, 1);
    meanCorrDat = mean(corrData, 2);
end

%%

% Change group
disp('Changing group')
group = session.getAdminService.lookupGroup(groupName);
session.setSecurityContext(group);

% List project or dataset IDs
project = getProjects(session, ids, false);
if ~isempty(project)
    fprintf(1, 'Loading datasets from project %g\n', ids);
    datasets = toMatlabList(project.linkedDatasetList);
    datasetIds = arrayfun(@(x) x.getId().getValue(), datasets)';
else
    datasetIds = ids;
end

try
    for datasetId = datasetIds
        analyzeDatasetFlux(session, datasetId, meanCorrDat,...
            'invalidIds', invalidIds);
    end
catch ME
    fprintf(1, 'Error while analyzing dataset %g: %s\n', datasetId, ME.message);
end

%%
disp('Closing session');
client.closeSession
