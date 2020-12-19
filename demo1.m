%%
clear;close all;clc
 
%% datasets
% data_names = {'AGG','Flame','Spiral','Jain','2G','2G_unbalance','S1','UB','R15','D1024','3Circles','compound','pathData','S1_001S1','dig1-10_uni','PalmData25_uni','Coil20Data_25_uni','Cancer','Iris','Seeds','Banknote','AML28','PMF','Mfeat','USPS'};
% data_names = {'AGG','Flame','Spiral','Jain','2G','2G_unbalance','S1','R15','D1024','3Circles','compound','S1_001S1'};
% data_names = {'AGG','Flame','Spiral','Jain','2G','2G_unbalance','S1','UB','R15','3Circles','compound','pathData','S1_001S1'};
data_names = {'AGG','Flame','Spiral','Jain','2G','2G_unbalance','S1','R15','3Circles'};

data_names={'AGG','Flame','Spiral','Jain','2G','2G_unbalance','S1','UB','R15','D1024','3Circles','compound','S1_001S1'}
% data_names={'AGG','Flame','Spiral','Jain','2G','2G_unbalance','S1','UB','R15','D1024','3Circles','compound','S1_001S1','Iris','Seeds','Banknote','AML28','PMF','TSNE20_single_cell_mrna_pollen','TSNE20_USPS','TSNE20_MNIST_smallScale','TSNE20_dig1-10_uni','TSNE20_PalmData25_uni','TSNE20_Coil20Data_25_uni','TSNE20_orlFace','TSNE20_Mfeat','TSNE20_COIL100','TSNE20_Isolet_all','TSNE20_zip','TSNE20_semeion'};
 
% data_names = {'UB'}; % Note: 'UB' can be handled by setting minsize=0.0018*N for 'auto' cut,whereas can be hanlded by setting minsize = 0.018*N for manual cut
% data_names = {'S1_001S1'};
% data_names = {'TB1M', 'SF2M'};
% data_names = {'MNIST'};
% data_names = {'COIL100'};
% data_names = {'1.3m_mouse_brain_cell_from_10x_after_PCA'};
% data_names ={'cytof_h1','cytof_h2','cytof_one','Samusik_01','bipolar','retina'};
% data_names ={'TB1M', 'SF2M', 'CC5M', 'CG10M', 'Flower20M'};
% data_names = {'compound','pathData'};
% data_names = {'data_Flower_100000000'};  % out-of-memory in the knn graph construction step.

% data_names = {'PalmData25_uni','Coil20Data_25_uni','COIL100','Cancer','Iris','Seeds','Banknote'};
% data_names = {'Samusik_01','Samusik_all','bipolar','retina','mouse_scRNA'};
% data_names = {'PenDigits'}; % 100 million
% data_names =[{'PenDigits','USPS', 'Letters', 'MNIST','Covertype','PalmData25_uni','Coil20Data_25_uni','COIL100','Cancer','Iris','Seeds','Banknote','AML28','PMF','Mfeat'},strcat('RCCdata_',{'MNIST','Coil100','YTF','YaleB','reuters','pendigits','shuttle','miceprotein'})];
%  data_names ={'PenDigits','MNIST','Covertype','cytof_h1','cytof_h2','cytof_one','Samusik_01','Samusik_all'};
% data_names = {'PenDigits','USPS', 'Letters', 'MNIST','Covertype','cytof_h1','cytof_h2','cytof_one','Samusik_01','Samusik_all'};
% data_names = {'PenDigits','MNIST','cytof_h1','cytof_h2','cytof_one'};
% data_names = {'Samusik_all'}; % not appropriate beacuase the dataset is sparse data.
% data_names = {'PalmData25_uni','Coil20Data_25_uni','COIL100','connect-4','covtype.scale.first10features','Sensorless'};

% % data_names = {'cytof_h1','cytof_h2','cytof_one','Samusik_01','Samusik_all','Levine_32dim','Levine_13dim','CellCycle','colon','muscle'};
data_names = {'gauss_spiral_circle_dataWithLabel','gauss_spiral_circle_data_in_noiseWithLabel'};
% data_names = {'Chameleon1WithLabel','Chameleon2WithLabel','Chameleon3WithLabel','Chameleon4WithLabel','gauss_spiral_circle_data_in_noiseWithLabel'};
 
%% method and setting
method_names = {'LDP-MST','FastLDPMST'}; % 'LDP-MST', 'FastLDPMST', 'QT_app_v4_usingModifiedDensity','QT_HighDim','NN_MST','NN_MST2',QT_app_v5_usingModifiedDensity
densityMethod = 'NaturalNeighborWithThreshold'; %   'NaturalNeighborBased','NaturalNeighborWithThreshold',or 'Simple-knn-based'
cutMethod = 'a1'; % 'a1':  nodeWeight-based autoCut; 'a2' % unionFind-based autoCut (better than 'a1'); 'm': manual cut; 'a2G','mG'
ratio = 0.018; %  [0.01,0.02] is recommended; not needed for manual cutting;
distance_function = 'euclidean'; % Note: % for hnsw, its distance function only supports: 'euclidean','l2','cosine','ip';
LDP_method = 'GGA'; % NNA or GGA.
logplot = 1; % Not Needed for autoCut: 1, plot decision graph in logscale; 0, not in logscale

plot_flag = 1;
%% plot setting

figure('visible','on');subplot_id = 0;
ha = tight_subplot(length(method_names)+1,length(data_names)/1,[.02, .02],[0.1, .05],[0.1, 0.1]); %tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%% start
record_num = 0;
for name_id=1:length(data_names)
    clear data annotation_data
    dataName = data_names{name_id};
    disp([num2str(name_id),', ',dataName,':'])
    [data,annotation_data,ClustN,dataName] = load_dataset(dataName);
    [N,dim]=size(data);
    
    if plot_flag
        if dim == 2
            disp('plot input dataset...')
            subplot_id = subplot_id + 1;
            if N>100000
                idx = randperm(N,10000);
                axes(ha(subplot_id));scatter(data(idx,1),data(idx,2),10,'k','filled')
                axis tight; set(gca,'xtick',[],'ytick',[],'FontSize',10,'Linewidth',.01);box on;
            else
                axes(ha(subplot_id));
                scatter(data(:,1),data(:,2),3,'k','filled');
                axis tight; set(gca,'xtick',[],'ytick',[],'FontSize',10,'Linewidth',.01);box on;
                title(strrep(dataName,'_','\_'))
            end
            if subplot_id == 1
                ylabel('Dataset','FontSize',12,'FontWeight','bold')
            end
        end
    end   
    
    minsize=ratio*N;
    initial_max_k = ceil(log2(N));
    %% compare with our faster version
    if dim< 20
        knnMethod = 'kd_tree';         disp('kd-tree exact fast knn searching technique is used for the dataset with Dimension lower than 20')
    else
        knnMethod = 'hnsw';            disp('hnsw (L2 distance) approximate fast knn searching technique is used for the dataset with Dimension larger than 20')
    end
    
    % test
    for method_id = 1:length(method_names)
        method = method_names{method_id};
        switch method
            case 'LDP-MST' % Cheng's method
                [Label,time] = LDPMST_cheng(data, ClustN, minsize);   %% code by Cheng
            case 'FastLDPMST'  
                [Label,time] = LDPMST_QT_exact_fast_less_memory_v2(data, ClustN, minsize,knnMethod,initial_max_k); %%
            otherwise
                error('method is not included...please name the method appropriately.')
        end
    
        %% evaluate and plot 
        diff_colors = linspecer(length(unique(Label)));
        if plot_flag
            if dim == 2
                disp('plot clustering result...') 
                if N>100000 
                    axes(ha(subplot_id + length(data_names)*method_id));  
                    scatter(data(idx,1),data(idx,2),3,diff_colors(Label(idx),:),'filled')
                    axis tight; set(gca,'xtick',[],'ytick',[],'FontSize',10,'Linewidth',.01);box on;
                else 
                    axes(ha(subplot_id + length(data_names)*method_id));
                    scatter(data(:,1),data(:,2),3,diff_colors(Label,:),'filled')
                    axis tight; set(gca,'xtick',[],'ytick',[],'FontSize',10,'Linewidth',.01);box on;
                end               
            end
            if subplot_id == 1
                ylabel(strrep(method_names{method_id},'_','\_'),'FontSize',12,'FontWeight','bold')
            end
        end
        record_num = record_num + 1;
        Result_all(record_num).data_name = data_names{name_id};
        Result_all(record_num).N = N;
        Result_all(record_num).DIM = dim;
        %     Result_all(record_num).memUsed = memUsed; % memUsed is a self-defined function: Problem: when opening matlab, even we do nothing, this is a big value (strange)
        if exist('annotation_data','var')
            disp('evaluate the clustering result...')
            if any(isnan(annotation_data))
                id_Nan = isnan(annotation_data);
                annotation_data = annotation_data(~id_Nan);
                Label = Label(~id_Nan);
            end
            [NMI_max,ARI,NMI_sqrt,AMI]= NMI_ARI(Label,annotation_data);
            Result_all(record_num).NMI_max = round(NMI_max*1000)/1000;
            Result_all(record_num).NMI_sqrt = round(NMI_sqrt*1000)/1000;
            Result_all(record_num).ARI = round(ARI*1000)/1000;
            Result_all(record_num).AMI = round(AMI*1000)/1000;
            Result_all(record_num).PC = length(unique(Label));
            Result_all(record_num).TC = length(unique(annotation_data));
            disp([Result_all(record_num).data_name, ': ARI = ', sprintf('%0.4f',Result_all(record_num).ARI),',  Runtime = ',sprintf('%.3f',time),'sec'])
        end
        
        Result_all(record_num).time = [sprintf('%.3f',time),'sec'];
        Result_all(record_num).cutMethod = cutMethod;
        if strcmp(cutMethod, 'm1')
            Result_all(record_num).ratio = 'NotNeeded';
        else
            Result_all(record_num).ratio = ratio;
        end
 
        Result_all(record_num).LDP_method = LDP_method;
 
        Result_all(record_num).densityMethod = densityMethod;
        Result_all(record_num).method = method;
        disp([Result_all(record_num).data_name, ',  Runtime = ',sprintf('%.3f',time),'sec'])
        
        LabelArray = Label;
        
    end
end
%%
disp(' ******************** All Results ************************ ')
if exist('Result_all','var')
    disp(struct2table(Result_all, 'AsArray', true)) %struct2table function may not exist in low matlab version; if so, then use the following commented codes
    %     disp('NMI   Time(s)')
    %     for i = 1:length(Result_all)
    %         disp([Result_all(i).data_name, ':', sprintf('%8.2f',Result_all(i).NMI),[sprintf('%8.2f',Result_all(i).time)],'sec'])
    %     end
end
 
