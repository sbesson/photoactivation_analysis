function saveAndUploadEPS(fig, path, filename, session, id, ns)

fullPath = fullfile(path, [filename '_' num2str(id) '.eps']);
print(fig, '-depsc', fullPath);
close(fig);

nTries = 5;

% Save projections onto the server
fa = getImageFileAnnotations(session, id, 'include', ns);
if isempty(fa),
    fa = writeFileAnnotation(session, fullPath, 'namespace', ns);
    fprintf(1, 'Created file annotation %g ', fa.getId().getValue());
    linkAnnotation(session, fa, 'image', id);
    fprintf(1, 'and linked it to image %g\n', id);
else
    originalFile =  fa.getFile();
    fprintf(1, 'Updating file: %g\n', originalFile.getId().getValue());
    fprintf(1, 'Old SHA1: %s\n', char(originalFile.getSha1().getValue()));
    originalFile = updateOriginalFile(session, originalFile, fullPath);
    fprintf(1, 'New SHA1: %s\n', char(originalFile.getSha1().getValue()));
end