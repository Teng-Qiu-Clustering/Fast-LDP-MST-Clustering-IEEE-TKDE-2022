function [data,annotation_data,ClustN,dataName] = load_dataset(dataName)
if any(strcmp(dataName,{'PenDigits','MNIST'})) 
    load(['data_',dataName,'.mat'],'gt','fea')
    annotation_data = gt;
    data = fea;
    clear fea; clear gt;
elseif any(strcmp(dataName,{'cytof_h1','cytof_h2','cytof_one','Samusik_01','Samusik_all'}))
    load([dataName,'.mat'],'data','annotation_data') 
elseif any(strcmp(dataName,{'gauss_spiral_circle_dataWithLabel','gauss_spiral_circle_data_in_noiseWithLabel','AGG','Flame','Spiral','Jain','2G','2G_unbalance','S1','R15','3Circles','S1_001S1'}))
    load(['Random_permuted_', dataName,'.mat'],'data','annotation_data')
elseif any(strcmp(dataName,{'data_TB1M'}))
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
