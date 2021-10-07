close all;
addpath(genpath(pwd));
disp('first generate a set of datasets (by sampling different number of samples from 5 different sampling functions)...')
data_names ={'TB', 'SF', 'CC', 'CG', 'Flower'};
for name_id=1:length(data_names)
    exponents = 14:20;
    for exponent_id = 1:length(exponents)
        dataSize = 2.^exponents(exponent_id);
        dataName = data_names{name_id};
        disp(['generate dataset: ',dataName,'_',num2str(dataSize),':'])
        fileName = ['data_',data_names{name_id},'_',num2str(dataSize),'.mat'];
        if ~exist(fileName,'file')
            figure;
            synthesizeLargescaleDatasets_withArbitrarySizes(dataName,dataSize);
            close(gcf)
        end
    end
end

% start test
disp('start test...')

 addpath(genpath(pwd));
N_array = 2.^(14:20);
method_array = {'FastLDPMST'};
dataNames = {'TB','SF','CC','CG','Flower'};

for data_id = 1:length(dataNames)
    dataName =dataNames{data_id};
    disp(dataName);
    result = []; time = 0; M_array = [];
    for m_id = 1:length(method_array)
        method = method_array{m_id};
        time = 0;
        for i = 1:length(N_array)
            disp(N_array(i))
            
            fileName = ['data_',dataName,'_',num2str(N_array(i)),'.mat'];
            if exist(fileName,'file')
                load(fileName)
            else
                figure;
                synthesizeLargescaleDatasets_withArbitrarySizes(dataName,N_array(i));
                close(gcf)
				load(fileName)
            end
            
            annotation_data = gt;
            data = fea;
            
            ClustN = length(unique(annotation_data));
            switch method
                case 'FastLDPMST'
                    [Label,time,Ini_clusterNum,supk,t_clusterMergeDistance,t_treeConstruction,t_edgeCutting] = FastLDPMST(data, ClustN); %%
                    M_array = [M_array Ini_clusterNum];
            end
            
            [NMI,ARI]= NMI_ARI(Label,annotation_data);
            result(m_id,i,1) = supk;
            result(m_id,i,2) = Ini_clusterNum;
            result(m_id,i,3) = t_clusterMergeDistance;
            result(m_id,i,4) = t_treeConstruction;
            result(m_id,i,5) = t_edgeCutting;
            result(m_id,i,6) = ARI;
            result(m_id,i,7) = time; % total time
            result(m_id,i,8) = N_array(i);
        end
    end
    
    %% 
    h = figure('position',[282.6 240.2 765.4 521.8]); fonS =12;
    markers = {'o','+','x','square','v','^','diamond','pentagram'};
    
    subplot(2,3,1);
    p = loglog(N_array,result(:,:,2),'o-');
    ylabel('M','FontSize',fonS);
    for i = 1:length(method_array)
        p(i).Marker = markers{i};
    end
    set(gca,'FontSize',fonS);
    % xlim([min(N_array)-1,max(N_array) + 10])
    xlabel({'N';'(a)'})
    % set(gca,'xtick',10.^(2:2:ceil(max(log10(N_array)))))
    set(gca,'xtick',10.^(round(min(log10(N_array))):ceil(max(log10(N_array)))))
    set(gca,'ytick',10.^(floor(min(log10(M_array))):ceil(max(log10(M_array)))))
    
    subplot(2,3,2);
    p = semilogx(M_array,result(:,:,3),'o-');
    ylabel('Runtime of step 5','FontSize',fonS);
    set(gca,'yscale','log');
    hold on;
    plot(M_array,M_array.*log(M_array),'--k');
    plot(M_array,M_array.*M_array,'-.k');
    text(M_array(4),8*M_array(4).*log(M_array(4)),'M\cdotlog(M)','rotation',15)
    text(M_array(4),10*M_array(4).*M_array(4),'M^2','rotation',30)
    xlim([min(M_array),max(M_array)])
    set(gca,'xtick',10.^(floor(min(log10(M_array))):ceil(max(log10(M_array)))))
    set(gca,'ytick',10.^(0:4:ceil(max(log10(M_array.*M_array)))))
    %     ylim([min(result(:,3))-0.1,max(result(:,3))])
    set(gca,'FontSize',fonS);
    for i = 1:length(method_array)
        p(i).Marker = markers{i};
    end
   
    xlabel({'M';'(b)'}) 
    
    subplot(2,3,3);
    p = semilogx(M_array,result(:,:,4),'o-');
    ylabel('Runtime of step 6','FontSize',fonS);
    set(gca,'yscale','log');
    hold on;
    plot(M_array,M_array.*log(M_array),'--k');
    plot(M_array,M_array.*M_array,'-.k');
    text(M_array(4),8*M_array(4).*log(M_array(4)),'M\cdotlog(M)','rotation',15)
    text(M_array(4),10*M_array(4).*M_array(4),'M^2','rotation',30)
    xlim([min(M_array),max(M_array)])
    set(gca,'xtick',10.^(round(min(log10(M_array))):ceil(max(log10(M_array)))))
    set(gca,'ytick',10.^(0:4:ceil(max(log10(M_array.*M_array)))))
    %     ylim([min(result(:,4))-0.1,max(result(:,4))])
    set(gca,'FontSize',fonS);
    for i = 1:length(method_array)
        p(i).Marker = markers{i};
    end
    title('  ','FontSize',20)
    % legend(method_array,'Position', [0.31,0.92,0.43,0.03],'Units', 'normalized','Orientation','horizontal')
    xlabel({'M';'(c)'})
    
    subplot(2,3,4);
    p = semilogx(M_array,result(:,:,5),'o-');
    ylabel('Runtime of step 7','FontSize',fonS);
    set(gca,'yscale','log');
    hold on;
    plot(M_array,M_array.*log(M_array),'--k');
    plot(M_array,M_array.*M_array,'-.k');
    text(M_array(4),8*M_array(4).*log(M_array(4)),'M\cdotlog(M)','rotation',15)
    text(M_array(4),10*M_array(4).*M_array(4),'M^2','rotation',30)
    xlim([min(M_array),max(M_array)])
    
    %     ylim([min(result(:,5))-0.1,max(result(:,5))])
    set(gca,'xtick',10.^(round(min(log10(M_array))):ceil(max(log10(M_array)))))
    set(gca,'ytick',10.^(0:4:ceil(max(log10(M_array.*M_array)))))
    
    set(gca,'FontSize',fonS);
    for i = 1:length(method_array)
        p(i).Marker = markers{i};
    end
    title('  ','FontSize',20)
    % legend(method_array,'Position', [0.31,0.92,0.43,0.03],'Units', 'normalized','Orientation','horizontal')
    xlabel({'M';'(d)'})
    
    subplot(2,3,5);
    p = semilogx(N_array,result(:,:,7),'o-');
    set(gca,'yscale','log');
    hold on; plot(N_array,N_array.*log(N_array),'--k');
    plot(N_array,N_array.*N_array,'-.k');
    text(N_array(3),8*N_array(3).*log(N_array(3)),'N\cdotlog(N)','rotation',15)
    text(N_array(3),10*N_array(3).*N_array(3),'N^2','rotation',30)
    xlabel('N','FontSize',fonS);ylabel('Total runtime','FontSize',fonS);
    for i = 1:length(method_array)
        p(i).Marker = markers{i};
    end
    set(gca,'FontSize',fonS);
    % xlim([min(N_array)-1,max(N_array) + 10])
    set(gca,'xtick',10.^(round(min(log10(N_array))):ceil(max(log10(N_array)))))
    set(gca,'ytick',10.^(0:4:ceil(max(log10(N_array.*N_array)))))
    %  ylim([min(result(:,7))-0.1,max(result(:,7))])
    xlabel({'N';'(e)'})
    
    subplot(2,3,6);
    p = semilogx(N_array,result(:,:,6),'o-');
    ylabel('ARI','FontSize',fonS);
    % xlim([min(N_array)-1,max(N_array) + 10])
    ylim([max(0,min(min(result(:,:,6)))-0.1) 1.01])
    
    set(gca,'FontSize',fonS);
    for i = 1:length(method_array)
        p(i).Marker = markers{i};
    end
    title('  ','FontSize',20)
    legend(method_array,'Position', [0.31,0.92,0.43,0.03],'Units', 'normalized','Orientation','horizontal')
    % legend(method_array,'Position', [0.31,0.92,0.43,0.03],'Units', 'normalized','Orientation','vertical')
    xlabel({'N';'(f)'})
    set(gca,'xtick',10.^(round(min(log10(N_array))):ceil(max(log10(N_array)))))
    
    result_file = ['demo2_result_',dataName,'.png'];
    disp(['Save plot result to ',result_file])
    saveas(gcf,result_file)
    disp(' '); disp(' ');
end
