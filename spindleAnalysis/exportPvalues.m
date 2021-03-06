% Load saved analysis data
s = load('~/Desktop/Ammended data - accounts for EMCCD swap/analysis.mat');
conditions = s.conditions;

for i =1: numel(conditions)
    index = conditions(i).alignedTimes == 70;
    conditions(i).data70 = conditions(i).alignedData(index, :);
    conditions(i).data70(isnan(conditions(i).data70)) = [];
end

%% T-Test
p = zeros(numel(conditions));
fid = fopen('~/Desktop/Ammended data - accounts for EMCCD swap/ptable_ttest.txt','w');
for i = 1: numel(conditions)
    fprintf(fid, '\t%s', conditions(i).name);
end
fprintf(fid, '\n');

for i =1: numel(conditions)
    fprintf(fid, '%s', conditions(i).name);
    for j = 1: i -1
        [~, p(i,j)] = ttest2(conditions(i).data70, conditions(j).data70);
        fprintf(fid, '\t%g', p(i,j));
    end
    fprintf(fid, '\n');
end
fclose(fid);

%% ANOVA

% Reorder conditions to list control first
index = 1 : numel(conditions);
iControl = find(strcmp({conditions.name},'control'));
index(iControl) = [];
index = [iControl index];
conditions = conditions(index);

% Combine data into groups and perform ANOVA test
raw_data=[conditions(:).data70];
groupNames=arrayfun(@(i)repmat({conditions(i).name},1,  numel(conditions(i).data70)), 1:18,'unif',0);
groupNames=horzcat(groupNames{:});
[p,t,stats] = anova1(raw_data,groupNames);
[c,m,h,nms] = multcompare(stats);

%%
p_anova = zeros(numel(conditions));
fid = fopen('~/Desktop/Ammended data - accounts for EMCCD swap/ptable_anova.txt','w');
for i = 1: numel(conditions)
    fprintf(fid, '\t%s', conditions(i).name);
end
fprintf(fid, '\n');

for i =1: numel(conditions)
    fprintf(fid, '%s', conditions(i).name);
    for j = 1: i -1
        fprintf(fid, '\t%g', c(find(c(:,1)==j & c(:,2)==i, 1), 6));
    end
    fprintf(fid, '\n');
end
fclose(fid);
