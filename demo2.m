%% Generate test datasets
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

disp('Start testing ...')
%% start testing
method_names = {'FastLDPMST'};
%% start testing
plot_flag = 1;

Result_all = [];test_num = 1;

for method_id = 1:length(method_names)
    method = method_names{method_id};
    for name_id=1:length(data_names)
        
        
        exponents = 14:20;
        for exponent_id = 1:length(exponents)
            %% load dataset
            clear data annotation_data
            dataSize = 2.^exponents(exponent_id);
            dataName = data_names{name_id};
            disp(['data ', num2str(test_num),': ',dataName,'_',num2str(dataSize),':'])
            fileName = ['data_',data_names{name_id},'_',num2str(dataSize),'.mat'];
            if ~exist(fileName,'file')
                figure;
                synthesizeLargescaleDatasets_withArbitrarySizes(dataName,dataSize);
                close(gcf)
            end
            load(fileName);
            
            annotation_data = gt;
            data = fea;
            clear fea; clear gt;
            if (min(annotation_data) == 0)
                annotation_data = annotation_data+1;
            end
            
            data = double(data);
            [N,dim]=size(data);
            nC = length(unique(annotation_data)); % number of clusters
            disp(['#objects: ',num2str(N),'; #features: ',num2str(dim)])
            
            %% parameter setting
            ratio = 0.018; %  [0.01,0.02] is recommended; not needed for manual cutting;
            MinSize= ratio*N;
            %% compared methods
            switch method
                case 'FastLDPMST'
                    [Label,time] = FastLDPMST(data, nC, MinSize); %%
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
            disp(" "); disp(" ");
            if time > 1800 % half an hour
                break;
            end
        end
    end
    %% save result
    disp(struct2table(Result_all, 'AsArray', true))
   
end
 