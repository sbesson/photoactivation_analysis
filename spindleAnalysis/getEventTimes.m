function [events, log] = getEventTimes(events, times)

% Read event times
nEvents = numel(events);
log = cell(nEvents, 1);
for i = 1 : nEvents
    for j = 1 : numel(events(i).index)
        events(i).times(j) = times(events(i).index(j));
    end
    % Read event times
    log{i} = sprintf([events(i).name ' mean time: %g +/- %g\n'], ...
        mean(events(i).times), std(events(i).times));
end
log = [log{:}];    