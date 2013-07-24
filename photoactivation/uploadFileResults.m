function fa = uploadFileResults(session, path, type, id, namespace)

% Input check
ip = inputParser;
ip.addRequired('filepath', @(x) ischar(x) && exist(x, 'file') == 2);
ip.addRequired('type', @(x) ismember(x, {'image', 'dataset'}));
ip.addRequired('id', @isscalar);
ip.addRequired('namespace', @ischar);
ip.parse(path, type, id, namespace);

% Save projections onto the server
if strcmp(type, 'image')
    fa = getImageFileAnnotations(session, id, 'include', namespace);
elseif strcmp(type, 'dataset')
    fa = getDatasetFileAnnotations(session, id, 'include', namespace);
end

if isempty(fa),
    fa = writeFileAnnotation(session, filepath, 'namespace', namespace);
    fprintf(1, 'Created file annotation %g ', fa.getId().getValue());
    linkAnnotation(session, fa, type, id);
    fprintf(1, 'and linked it to %s %g\n', type, id);
else
    fprintf(1, 'Updating file annotation: %g\n', fa.getId().getValue());
    fprintf(1, 'Old SHA1: %s\n', char(fa.getFile().getSha1().getValue()));
    fa = updateFileAnnotation(session, fa, filepath);
    fprintf(1, 'New SHA1: %s\n', char(fa.getFile().getSha1().getValue()));
end