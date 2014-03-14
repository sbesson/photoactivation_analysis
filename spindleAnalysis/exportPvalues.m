load('~/Desktop/Ammended data - accounts for EMCCD swap/analysis.mat');

for i =1: numel(conditions)
%     for j = 1: numel(conditions(i).series)
    index = conditions(i).alignedTimes == 70;
    conditions(i).data70 = conditions(i).alignedData(index, :);
    conditions(i).data70(isnan(conditions(i).data70)) = [];
%     end
end

%%
p = zeros(numel(conditions));
fid = fopen('~/Desktop/Ammended data - accounts for EMCCD swap/ptable.txt','w');
for i =1: numel(conditions)
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

