% Configuration

% OMERO.matlab configuration file
configFile = fullfile(fileparts(which('startup.m')), 'ice.config');

% Project to use
groupName = 'Rape project';

% Define projects or datasets to be analyzed
projectId = 651;
datasetIds = [];

% Wrong images (to exclude from the analysis)
invalidIds = [3421,3460];

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

if ~isempty(projectId),
    project = getProjects(session, projectId, false);
    datasets = toMatlabList(project.linkedDatasetList);
    datasetIds = arrayfun(@(x) x.getId().getValue(), datasets)';
end

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
