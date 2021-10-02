%%
clear;close all;clc;  
addpath(genpath(pwd));
%% Datasets 
data_names ={'TB', 'SF', 'CC', 'CG', 'Flower'}; 
% data_names ={'TB'}; 
%% Methods 
method_names = {'FastLDPMST'};  
%% start testing
plot_flag = 1;
for name_id=1:length(data_names)
    Result_all = [];test_num = 1; 
    
    for method_id = 1:length(method_names)
        method = method_names{method_id};
        
        exponents = 14:20;
        for exponent_id = 1:length(exponents)
           %% load dataset
            clear data annotation_data
            dataSize = 2.^exponents(exponent_id);          
            dataName = data_names{name_id};
            disp([dataName,'_',num2str(dataSize),':'])            
            fileName = ['data_',data_names{name_id},'_',num2str(dataSize),'.mat'];
            if ~exist(fileName,'file') 
                figure;
                synthesizeLargescaleDatasets_withArbitrarySizes(dataName,dataSize);
                close(gcf)
            end            
            load(fileName); pause(2);
            
            annotation_data = gt;
            data = fea;
            clear fea; clear gt;  pause(2)
            if (min(annotation_data) == 0)
                annotation_data = annotation_data+1;
            end 
            
            data = double(data);
            [N,dim]=size(data);
            ClustN = length(unique(annotation_data));
            disp(['#objects: ',num2str(N),'; #features: ',num2str(dim)])
            
            %% parameter setting
            ratio = 0.018; %  [0.01,0.02] is recommended; not needed for manual cutting;
            mS= ratio*N;  
            %% compared methods
            switch method 
                case 'FastLDPMST' 
                    [Label,time] = FastLDPMST(data, ClustN, mS); %%
                otherwise
                    error('method is not included...please name the method appropriately.')
            end
            %% Evaluate results
            clear data
            Result_all(test_num).data_name = [dataName,'_',num2str(dataSize)];
            Result_all(test_num).N = N;
         
            if exist('annotation_data','var')
                disp('evaluate the clustering result...')
                [NMI,ARI]= NMI_ARI(Label,annotation_data);
                Result_all(test_num).NMI = roundn(NMI,-3);
                Result_all(test_num).ARI = roundn(ARI,-3);
            end
            Result_all(test_num).time = roundn(time,-3);
            Result_all(test_num).method = method;
            test_num = test_num + 1;
            if time > 1800 % half an hour
                break;
            end
        end
    end 
    %% save result 
    save(['compare_varying_size_',dataName,'.mat'],'Result_all','method_names','data_names','exponents')
end
 
%% plot comparison on scores
close all
figure('position',[52.2000 61.8000 1415 712.8000]);   
fonS = 17;
ha = tight_subplot(2,length(data_names),[.08, .05],[0.15, .04],[0.06, 0.01]); %tight_subplot(Nh, Nw, gap, marg_h, marg_w)

subplot_id = 1;
for name_id = 1:length(data_names) 
    dataName = data_names{name_id};  
    load(['compare_varying_size_',dataName,'.mat']) 
    method_all = [];
    data_name_all = [];
    for j = 1:length(Result_all)
        method_all{j} = Result_all(1,j).method;
        data_name_all{j} = Result_all(1,j).data_name;
    end
    
    for j = 1:length(method_all)
        switch method_all{j} 
            case 'FastLDPMST'
                method_all{j} = 'FastLDPMST'; 
        end
    end
    
    NMI_all = [Result_all.NMI];
    ARI_all = [Result_all.ARI];
    time_all = [Result_all.time];
 
    markers = {'+','o','*','x','s','d','^','v','>','<','p','h'};
    axes(ha(subplot_id));    
    
    %% Runtime
    method_names = unique(method_all);
    for method_id = 1:length(method_names)
        method = method_names{method_id};
        id = find(strcmp(method_all, method));
        x = 2.^exponents; 
        y = time_all(id);     
        loglog(x,y,'-','Marker',markers{mod(method_id,numel(markers))},'MarkerSize',7);
      
        hold on
    end
    
   
    set(gca,'FontSize',fonS);
    %     xlabel('N','FontSize',fonS)
    if name_id == 1
        ylabel('Runtime (seconds)','FontSize',fonS)
    end
    title(dataName,'FontSize',fonS)
    set(gca,'xtick',10.^(round(min(log10(x))):ceil(max(log10(x))))) 
    set(gca,'ytick',10.^(round(min(log10(time_all))):ceil(max(log10(time_all))))) 
    
    hold off
    
    axes(ha(subplot_id + length(data_names)));     subplot_id = subplot_id + 1;
    
    %% ARI
    method_names = unique(method_all);
    for method_id = 1:length(method_names)
        method = method_names{method_id};
        id = find(strcmp(method_all, method));
        x = 2.^exponents;
        y = ARI_all(id); 
        semilogx(x,y,'-','Marker',markers{mod(method_id,numel(markers))},'MarkerSize',7);
       
        hold on
    end

    set(gca,'FontSize',fonS)
    xlabel('N','FontSize',fonS)
    if name_id == 1
        ylabel('ARI','FontSize',fonS)
    end
 
    ylim([min(ARI_all)-0.1,1])
     set(gca,'xtick',10.^(round(min(log10(x))):ceil(max(log10(x))))) 
 
    hold off
end
legend(strrep(method_names,'_','\_'),'FontSize',fonS+2,'Position', [0.31,0.01,0.43,0.03],'Units', 'normalized','Orientation','horizontal')
 
