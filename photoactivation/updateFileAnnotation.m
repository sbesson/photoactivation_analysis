function updateFileAnnotation(session, path, type, id, ns)

% Input check
ip = inputParser;
ip.addRequired('path', @(x) ischar(x) && exist(x, 'file') == 2);
ip.addRequired('type', @(x) ismember(x, {'image', 'dataset'}));
ip.addRequired('id', @isscalar);
ip.addRequired('ns', @ischar);
ip.parse(path, type, id, ns);

% Save projections onto the server
if strcmp(type, 'image')
    fa = getImageFileAnnotations(session, id, 'include', ns);
elseif strcmp(type, 'dataset')
    fa = getDatasetFileAnnotations(session, id, 'include', ns);
end

if isempty(fa),
    fa = writeFileAnnotation(session, path, 'namespace', ns);
    fprintf(1, 'Created file annotation %g ', fa.getId().getValue());
    linkAnnotation(session, fa, type, id);
    fprintf(1, 'and linked it to %s %g\n', type, id);
else
    originalFile =  fa.getFile();
    fprintf(1, 'Updating file: %g\n', originalFile.getId().getValue());
    fprintf(1, 'Old SHA1: %s\n', char(originalFile.getSha1().getValue()));
    originalFile = updateOriginalFile(session, originalFile, path);
    fprintf(1, 'New SHA1: %s\n', char(originalFile.getSha1().getValue()));
end