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
inputString = ['Specify the name of the group  [default: ' defaultgroup ']: '];
groupName = input(inputString, 's');
if isempty(groupName), groupName = defaultgroup; end


% Define id of projects or datasets to be analyzed
corrId = input('Specify the ids of the photobleaching dataset:');
assert(isempty(corrId) || isscalar(corrId), 'Only one photobleaching Id allowed');

% Define id of projects or datasets to be analyzed
ids = input('Specify the ids of the projects or datasets to analyze:');
assert(~isempty(ids), 'No ids selected');

% Wrong images (to exclude from the analysis)
invalidIds = [3991, 3412, photoactivationId];
%% OMERO.matlab initialization

% Create client/session
[client, session] = loadOmero(configFile);
fprintf(1, 'Created connection to %s\n', char(client.getProperty('omero.host')));
fprintf(1, 'Created session for user %s using group %s\n',...
    char(session.getAdminService().getEventContext().userName),...
    char(session.getAdminService().getEventContext().groupName));

% Change group
disp('Changing group')
group = session.getAdminService.lookupGroup(groupName);
session.setSecurityContext(group);

%%

% List project or dataset IDs
project = getProjects(session, ids, false);
if ~isempty(project)
    fprintf(1, 'Loading datasets from project %g\n', ids);
    datasets = toMatlabList(project.linkedDatasetList);
    datasetIds = arrayfun(@(x) x.getId().getValue(), datasets)';
else
    datasetIds = ids;
end

% Run the analysis on 
% if ~isempty(corrId)
%     phoactivationData = analyzeDatasetFlux(session, photoactivationId);
%     photoactivationCorrection = mean(horzcat(data.dInorm), 2);
% end

try
    for datasetId = datasetIds
        analyzeDatasetFlux(session, datasetId, invalidIds);
    end
catch ME
    fprintf(1, 'Error while analyzing dataset %g: %s\n', datasetId, ME.message);
end

%%
disp('Closing session');
client.closeSession
