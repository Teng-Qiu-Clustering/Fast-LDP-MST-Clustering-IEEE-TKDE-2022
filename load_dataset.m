function [data,annotation_data,ClustN,dataName] = load_dataset(dataName)
if strcmp(dataName,'PenDigits') || strcmp(dataName,'USPS') || strcmp(dataName,'Letters') || strcmp(dataName,'MNIST') || strcmp(dataName,'Covertype') || strcmp(dataName,'TB1M') || strcmp(dataName,'SF2M') || strcmp(dataName,'CC5M') || strcmp(dataName,'CG10M') || strcmp(dataName,'Flower20M')
    load(['D:\360¼«ËÙä¯ÀÀÆ÷ÏÂÔØ\matlab_and_dataset_for_MATLAB_source_code_for_U-SPEC_and_U-SENC\','data_',dataName,'.mat'],'gt','fea')
    annotation_data = gt;
    data = fea;
    clear fea; clear gt;
elseif any(strcmp(dataName,{'Test_1_mECS','Test_2_Kolod','Test_3_Pollen','Test_4_Usoskin','Test_5_Zeisel'}))
    load(['F:\Download_code\TPE\tpe\R\SIMLR-SIMLR\MATLAB\data\',dataName,'.mat'],'in_X','true_labs')
    data = in_X; annotation_data = true_labs;
elseif any(strcmp(dataName,{'cytof_h1','cytof_h2','cytof_one','Samusik_01','Samusik_all','bipolar','retina','mouse_scRNA','pbmc_68k_pca50vars','FlowCAP_ND','FlowCAP_WNV','Levine_32dim','Levine_13dim','CellCycle','colon','muscle'}))
    load(['F:\mc\mc\object_detection\c\GAP2\data_after_Processing\',dataName,'.mat'],'data','annotation_data')
elseif any(strcmp(dataName,strcat('RCCdata_',{'Coil100','miceprotein','MNIST','pendigits','reuters','shuttle','YaleB','YTF'})))
    load(['F:\Download_code\TPE\tpe\R\Data_from_RCC_PNAS\',dataName(9:end),'.mat'],'X','gtlabels')
    data = X; annotation_data = double(gtlabels');
elseif any(strcmp(dataName,{'TB', 'SF', 'CC', 'CG', 'Flower'}))
    load(['D:\large-scale-datasets\','data_',dataName,'_',num2str(dataSize),'.mat'],'gt','fea')
    annotation_data = gt;
    data = fea;
    clear fea; clear gt;
elseif strcmp(dataName,'data_Flower_100000000')
    load(['D:\large-scale-datasets\',dataName,'.mat'],'gt','fea')
    annotation_data = gt;
    data = fea;
    clear fea; clear gt;
else
    load(['F:\NND\8) in-tree in the kernel space\Datasets_random_permutated\','Random_permuted_', dataName,'.mat'],'data','annotation_data')
end

if exist('annotation_data','var')&&(min(annotation_data) == 0)
    annotation_data = annotation_data+1;
end

switch dataName
    case 'gauss_spiral_circle_data'
        ClustN = 9;  
    case 'gauss_spiral_circle_data_in_noise'
        ClustN = 9; 
    case 'Chameleon1'
        ClustN = 6;
    case 'Chameleon2'
        ClustN = 6;
    case  'Chameleon3'
        ClustN = 9;
    case  'Chameleon4'
        ClustN = 8; 
    otherwise
        ClustN = length(unique(annotation_data(~isnan(annotation_data))));
end

switch dataName % change name to be a more concise form
    case 'gauss_spiral_circle_data'
        dataName = 'GaSpCi';
    case 'gauss_spiral_circle_dataWithLabel'
        dataName = 'GaSpCi';
    case 'gauss_spiral_circle_data_in_noise'
        dataName = 'GaSpCiNo';
    case 'gauss_spiral_circle_data_in_noiseWithLabel'
        dataName = 'GaSpCiNo';
    case 'Chameleon1WithLabel'
        dataName = 'Chameleon1';
    case 'Chameleon2WithLabel'
        dataName = 'Chameleon2';
    case 'Chameleon3WithLabel'
        dataName = 'Chameleon3';
    case 'Chameleon4WithLabel'
        dataName = 'Chameleon4'; 
end

[N,dim]=size(data);
disp(['dataName: ',dataName, ';#objects: ',num2str(N),'; #features: ',num2str(dim),';#Clusters: ',num2str(ClustN)])
