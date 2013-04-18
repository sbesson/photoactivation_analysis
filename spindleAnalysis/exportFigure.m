function exportFigure(fig_handle, path, basename)

% Define formats for figure export
f(1).driver = '-dtiff';
f(1).ext = 'tif';
f(2).driver = '-depsc';
f(2).ext = 'eps';

for i = 1 : numel(f)
    print(fig_handle, f(i).driver, fullfile(path, [basename '.' f(i).ext]));
end