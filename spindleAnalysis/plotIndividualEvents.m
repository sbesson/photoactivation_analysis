function plotIndividualEvents(s)

% define small and large fonts for graphical output
%tfont = {'FontName', 'Helvetica', 'FontSize', 14, 'FontAngle', 'italic'};
sfont = {'FontName', 'Helvetica', 'FontSize', 18};
lfont = {'FontName', 'Helvetica', 'FontSize', 22};

% Plot indidividual events
plot(s.times, s.data)
hold on;
for i = 1 : s.nSamples
    plot(s.eventTimes(i), s.eventValues(i), 'ok', 'MarkerFaceColor', 'k');
end
set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
xlabel('Time (min)', lfont{:});
ylabel('Spindle length (\mum)', lfont{:});
