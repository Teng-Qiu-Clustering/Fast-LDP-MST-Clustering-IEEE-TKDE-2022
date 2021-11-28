function [data,annotation_data,ClustN,dataName] = load_data(dataName)
if any(strcmp(dataName,{'PenDigits','MNIST'})) 
    load(['data_',dataName,'.mat'],'gt','fea')
    annotation_data = gt;
    data = fea;
    clear fea; clear gt;
elseif any(strcmp(dataName,{'cytof_h1','cytof_h2','cytof_one','Samusik_01','Samusik_all'}))
    load([dataName,'.mat'],'data','annotation_data') 
elseif any(strcmp(dataName,{'gauss_spiral_circle_dataWithLabel','gauss_spiral_circle_data_in_noiseWithLabel','AGG','Flame','Spiral','Jain','2G','2G_unbalance','S1','R15','3Circles','S1_001S1'}))
    load(['Random_permuted_', dataName,'.mat'],'data','annotation_data')
elseif any(strcmp(dataName,{'data_TB1M','data_SF2M','data_TB_16384','data_TB_32768','data_TB_65536','data_TB_131072','data_TB_262144','data_TB_524288','data_TB_1048576','data_SF_16384','data_SF_32768','data_SF_65536','data_SF_131072','data_SF_262144','data_SF_524288','data_SF_1048576','data_CC_16384','data_CC_32768','data_CC_65536','data_CC_131072','data_CC_262144','data_CC_524288','data_CC_1048576','data_CG_16384','data_CG_32768','data_CG_65536','data_CG_131072','data_CG_262144','data_CG_524288','data_CG_1048576','data_Flower_16384','data_Flower_32768','data_Flower_65536','data_Flower_131072','data_Flower_262144','data_Flower_524288','data_Flower_1048576','data_TB_1000000', 'data_SF_1000000','data_CC_1000000', 'data_CG_1000000', 'data_Flower_1000000','data_TB_100000', 'data_SF_100000','data_CC_100000', 'data_CG_100000', 'data_Flower_100000'}))
    load([dataName,'.mat']) % 
    annotation_data = gt;
    data = fea;
    clear fea; clear gt;
else
    load([dataName,'.mat'],'data','annotation_data')
end

if exist('annotation_data','var')&&(min(annotation_data) == 0)
    annotation_data = annotation_data+1;
end

switch dataName
    case 'gauss_spiral_circle_data'
        ClustN = 9;
    case 'gauss_spiral_circle_data_in_noise'
        ClustN = 9;
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
end

[N,dim]=size(data);
disp(['dataName: ',dataName, ';#objects: ',num2str(N),'; #features: ',num2str(dim),';#Clusters: ',num2str(ClustN)])
