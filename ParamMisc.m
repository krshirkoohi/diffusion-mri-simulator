%% plot graphs and etc.
%% sim2
clear; loadFolder = dir("Results/Sim2");
fileNames = extractfield(loadFolder(3:length(loadFolder)),'name')'; E = length(fileNames);
SD = zeros(5,10);
%
for param=1:5
    figure, hold on,
    for trial=1:10
        load(char(fileNames((param-1)*10+trial)));
        semilogy(bVal,real(S)), axis([0.05, 2.5e11, 1e-2,0.3]);
        SD(param,trial) = std(real(S'));
        %disp(fileNames(param+trial));
    end
    %hold on,legend("#1","#2","#3","#4","#5","#6","#7","#8","#9","#10");
    mean_SD = mean(SD(param,:));
    subtitle(sprintf('N_{RW}=1e{%d}, N_{ii}=%d, N_t=%d, N_{G_s}=%d, \x3C3\x302=%.2d',log10(simu.N),N_ii,length(t),length(seq.G_s),mean_SD),'FontSize',9), ylabel('Signal attenuation'),xlabel('B-Value')
    %title(sprintf('N_{RW}=1e{%d}, N_{ii}=%d, N_t=1e{%d}, N_{G_s}=%d, \x3C3\x302=%.2d',log10(simu.N),N_ii,log10(length(t)),length(seq.G_s),mean_SD));
%     annotation('textbox',[.0 1.0 .0 .0], ...
%         'String',sprintf('N_{RW}=1e{%d}, N_{ii}=%d, N_t=1e{%d}, N_{G_s}=%d, \x3C3\x302=%.2d',log10(simu.N),N_ii,log10(length(t)),length(seq.G_s),mean_SD),'EdgeColor','none')
    %saveas(gcf,sprintf('Results/Sim2 Nii=%d Nrw=1e%d Nt=1e%d.png',log10(simu.N),N_ii,log10(length(t))))
end
%% sim1
clear; loadFolder = (dir("Results/Sim1"));
fileNames = extractfield(loadFolder(3:length(loadFolder)),'name')'; E = length(fileNames);
MRAE = zeros(E,1);
compTime = zeros(E,1);
for p=1:E
    load(char(fileNames(p)));
    MRAE(p) = mean(abs((Sa-real(S'))./Sa));
    compTime(p) = elapsedTime;
end
[errVals,idx] = sort(MRAE); compTime = compTime(idx);
% plot best
for q=1:5
    load(char(fileNames(idx(q))))
    figure('Position', [200 200 1000 400]), 
    subplot(1,3,1),semilogy(bVal,real(S)),hold on,semilogy(bVal,Sa,'r-'), %axis([0, 2.5e11, 1e-2,1])
    legend('Monte carlo','Analytical'),title('Results'),subtitle(sprintf('N_{RW}=1e{%d}, N_{ii}=%d, N_t=%d, N_{G_s}=%d',log10(simu.N),N_ii,length(t),length(seq.G_s)),'FontSize',7), ylabel('Signal attenuation'),xlabel('B-Value')
    % error
    subplot(1,3,2),plot(bVal,abs((Sa-real(S'))./Sa)),title('Relative error'), %axis([0, 2.5e11, 1e-2,1])
    %saveas(gcf, sprintf('Results/Sim1 Nii=%d Nrw=1e%d Nt=1e%d',log10(simu.N),N_ii,log10(length(t))),'png')
end




