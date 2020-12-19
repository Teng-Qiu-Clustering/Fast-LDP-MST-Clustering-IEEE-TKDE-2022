%%
clear;close all;clc

%% datasets
data_names={'AGG','Flame','Spiral','Jain','2G','2G_unbalance','S1','R15','3Circles','S1_001S1'};

% data_names = {'cytof_h1','cytof_h2','cytof_one','Samusik_01','Samusik_all','Levine_32dim','Levine_13dim','CellCycle','colon','muscle'};
% data_names = {'gauss_spiral_circle_dataWithLabel','gauss_spiral_circle_data_in_noiseWithLabel'};
data_names = {'PenDigits','MNIST'};
% data_names = {'data_TB1M'};
%% method and setting
method_names = {'FastLDPMST'};  
% method_names = {'LDP-MST','FastLDPMST'};
 
ratio = 0.018; %  [0.01,0.02] is recommended; not needed for manual cutting; 
%% start
record_num = 0;
for name_id=1:length(data_names)
    clear data annotation_data
    dataName = data_names{name_id};
    disp([num2str(name_id),', ',dataName,':'])
    [data,annotation_data,ClustN,dataName] = load_dataset(dataName);
    [N,dim]=size(data);
    
    minsize=ratio*N;
    initial_max_k = ceil(log2(N));
    %% compare with our faster version
    if dim< 20
        knnMethod = 'kd_tree';         disp('kd-tree exact fast knn searching technique is used for the dataset with Dimension lower than 20')
    else
        knnMethod = 'hnsw';            disp('hnsw (L2 distance) approximate fast knn searching technique is used for the dataset with Dimension larger than 20')
    end
     
    for method_id = 1:length(method_names)
        method = method_names{method_id};
        switch method
            case 'LDP-MST' % Cheng's method
                [Label,time] = LDPMST_cheng(data, ClustN, minsize);   %% code provided by Cheng
            case 'FastLDPMST'  
                [Label,time] = FastLDPMST(data, ClustN, minsize,knnMethod,initial_max_k); %%
            otherwise
                error('method is not included...please name the method appropriately.')
        end
        %% evaluate and plot
        diff_colors = linspecer(length(unique(Label)));
        
        if dim == 2
            figure;
            if N>100000
                idx = randperm(N,10000);
                subplot(1,2,1);scatter(data(idx,1),data(idx,2),10,'k','filled')
                subplot(1,2,2);scatter(data(idx,1),data(idx,2),3,diff_colors(Label(idx),:),'filled')
            else
                subplot(1,2,1);scatter(data(:,1),data(:,2),3,'k','filled');
                subplot(1,2,2);scatter(data(:,1),data(:,2),3,diff_colors(Label,:),'filled')
            end
            axis tight; set(gca,'xtick',[],'ytick',[],'FontSize',10,'Linewidth',.01);box on;
        end
        
        record_num = record_num + 1;
        Result_all(record_num).data_name = data_names{name_id};
        Result_all(record_num).N = N;
        Result_all(record_num).DIM = dim;
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
        Result_all(record_num).ratio = ratio; 
        Result_all(record_num).method = method;
        disp([Result_all(record_num).data_name, ',  Runtime = ',sprintf('%.3f',time),'sec'])
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

