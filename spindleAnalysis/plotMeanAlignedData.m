function plotMeanAlignedData(s)

% define small and large fonts for graphical output
%tfont = {'FontName', 'Helvetica', 'FontSize', 14, 'FontAngle', 'italic'};
sfont = {'FontName', 'Helvetica', 'FontSize', 18};
lfont = {'FontName', 'Helvetica', 'FontSize', 22};


nElements = numel(s);
colors = hsv(nElements);
tmin = Inf;
tmax = -Inf;
time_max = max([s.alignedTimes]);

for i = 1 : nElements
    
    % Calculate mean and std values of aligned data
    alignedmean = nanmean(s(i).alignedData,2);
    alignedstd = nanstd(s(i).alignedData, [], 2);
    nPoints  =sum(~isnan(s(i).alignedData),2);
    alignedste = alignedstd./sqrt(nPoints);
    
    % Calculate mininimum and maximum times
    imin = max(1,find(~isnan(alignedmean),1,'first'));
    tmin = min(tmin, s(i).alignedTimes(imin));
    imax = min(find(~isnan(alignedmean(imin:end)),1,'last'), time_max);
    tmax = max(tmax, s(i).alignedTimes(imin+imax));

    hold on;
    errorbar(s(i).alignedTimes', alignedmean, alignedste,...
        '-','Color', colors(i,:), 'LineWidth', 1);
end

set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
xlabel('Time (min)', lfont{:});
ylabel('Spindle length (\mum)', lfont{:});
xlim([tmin tmax]);
legend({s.name},'Location', 'NorthWest','Interpreter','None')
set(gca,'LooseInset',get(gca,'TightInset'));
box on