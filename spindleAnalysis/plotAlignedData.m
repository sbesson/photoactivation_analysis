function plotAlignedData(s)

% define small and large fonts for graphical output
%tfont = {'FontName', 'Helvetica', 'FontSize', 14, 'FontAngle', 'italic'};
sfont = {'FontName', 'Helvetica', 'FontSize', 18};
lfont = {'FontName', 'Helvetica', 'FontSize', 22};

plot(s.alignedTimes, s.alignedData)

% Calculate mean and std values of aligned data
alignedmean = nanmean(s.alignedData,2);
%alignedstd = nanstd(s.alignedData, [], 2);
%nPoints  =sum(~isnan(s.alignedData),2);
%alignedste = alignedstd./sqrt(nPoints);

% Calculate mininimum and maximum times
imin = max(1,find(~isnan(alignedmean),1,'first'));
imax = min(find(~isnan(alignedmean(imin:end)),1,'last')+imin, ...
    max(s.alignedTimes));
set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
xlabel('Time (min)', lfont{:});
ylabel('Spindle length (\mum)', lfont{:});
xlim([s.alignedTimes(imin) s.alignedTimes(imax)]);
