Ipath = '~/Downloads';
close all

corrI = dlmread(fullfile(Ipath, 'Intensities_1102.txt'),'\t',1,0);
corr_mean = mean(corrI(:,2:end),2);

outpath = fullfile(Ipath, 'fit.txt');
fid = fopen(outpath,'w');
    
f=figure; 
hold on
colors =hsv(4);
dids = [1101, 1107,1108, 1110];
    fn =  @(b,x) b(1) .* exp(-x/b(2)) + b(3) .* exp(-x/b(4)) +b(5);

for iDataset = 1:4
    did = dids(iDataset);
    dataPath = fullfile(Ipath, ['Intensities_' num2str(did) '.txt']);
    fid2 =fopen(dataPath,'r');
    header_line = fgetl(fid2);
    headers =  strsplit(header_line,'\t');
    headers = headers(2:end);
    fclose(fid2);
    I = dlmread(dataPath,'\t',1,0);
    times = I(:,1);
    I = I(:,2:end);
    Icorr = I./repmat(corr_mean, 1, size(I,2));
    
    fprintf(1, 'Fitting');
    fprintf(fid, '\tpf\ttf (s)\tps\tts (s)\n');
%     clear r2
%     for i = 1 : size(I,2)
%         fprintf(1,'.');
%         fprintf(fid, headers{i});
%         %         [b,r2(i)] = fitPADecay(times,Icorr(:,i));
%         try 
%         nlm = fitnlm(times, Icorr(:,i), fn ,[.8 16 .1 100 .1]);
%         r2(i)=nlm.Rsquared.Adjusted;
%         catch
%             r2(i)=0;
%         end
% %       if r2(fprintf(fid, '\t%g\t%g\t%g\t%g\n', b(1), b(2), b(3), b(4));
% 
%     end
    
if iDataset ==1
    index = Icorr(end,:)>.05 & Icorr(end,:)<.4;
elseif iDataset ==2
    index = Icorr(end,:)>.02 & Icorr(end,:)<.41;
elseif iDataset ==3
    index = Icorr(end,:)>.0 & Icorr(end,:)<.4;

else
    index = Icorr(end,:)>0 & Icorr(end,:)<.4;
end
%      index = r2 > .95;
    Icorr = Icorr(:,index);
    Imean = mean(Icorr,2);
    Iste = std(Icorr,[],2)./sqrt(size(Icorr,1)-1);
    
    fprintf(1,'\n');
    fprintf(fid, '\n\n');
    nlm = fitnlm(times, Imean, fn ,[.8 16 .1 100 .1]);
    %     [b, r2, ce] = fitPADecay(times, Imean);
    b=nlm.Coefficients.Estimate;
    fprintf(fid, 'Average fit\t%g\t%g\t%g\t%g\n',...
        b(1)/(b(1)+b(3)), b(2), b(3)/(b(1)+b(3)), b(4));
    fprintf(1, 'Average fit\t%g\t%g\t%g\t%g\n',...
        b(1)/(b(1)+b(3)), b(2), b(3)/(b(1)+b(3)), b(4));
    se=nlm.Coefficients.SE;

    fprintf(fid, 'Average fit SE\t%g\t%g\t%g\t%g\n', se(1), se(2), se(3), se(4));
    fprintf(fid, '\n\n');
    
%     fprintf(1,'\n');
%     fprintf(fid, '\n\n');
%     [b,~, ce] = fitPADecay(times, Imean);
%     fprintf(fid, 'Average fit\t%g\t%g\t%g\n',...
%         b(1)/(b(1)+b(3)), b(2), b(3));
%     fprintf(fid, 'Average fit CI\t%g\t%g\n', ce(1), ce(2), ce(3));
%     fprintf(fid, '\n\n');

    figure;
    plot(times,Icorr);
    figure(f); 
    errorbar(times,Imean, Iste, '-','Color', colors(iDataset,:));

end
legend({'Control', 'CLASP2', 'Kid', 'CLASP2+Kid'})
fclose(fid);