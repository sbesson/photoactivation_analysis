function saveAndUploadEPS(fig, path, filename, session, id, ns)

fullPath = fullfile(path, [filename '_' num2str(id) '.eps']);
print(fig, '-depsc', fullPath);
close(fig);

% Save projections onto the server
fa = getImageFileAnnotations(session, id, 'include', ns);
if isempty(fa),
    fa = writeFileAnnotation(session, fullPath, 'namespace', ns);
    linkAnnotation(session, fa, 'image', id);
else
    updateFileAnnotation(session, fa, fullPath);
end