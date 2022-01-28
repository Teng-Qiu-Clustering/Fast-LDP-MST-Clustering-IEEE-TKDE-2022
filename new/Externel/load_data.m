function [data,annotation_data,ClustN,dataName] = load_data(dataName)
if any(strcmp(dataName,{'data_TB1M','data_SF2M','data_TB_16384','data_TB_32768','data_TB_65536','data_TB_131072','data_TB_262144','data_TB_524288','data_TB_1048576','data_SF_16384','data_SF_32768','data_SF_65536','data_SF_131072','data_SF_262144','data_SF_524288','data_SF_1048576','data_CC_16384','data_CC_32768','data_CC_65536','data_CC_131072','data_CC_262144','data_CC_524288','data_CC_1048576','data_CG_16384','data_CG_32768','data_CG_65536','data_CG_131072','data_CG_262144','data_CG_524288','data_CG_1048576','data_Flower_16384','data_Flower_32768','data_Flower_65536','data_Flower_131072','data_Flower_262144','data_Flower_524288','data_Flower_1048576','data_TB_1000000', 'data_SF_1000000','data_CC_1000000', 'data_CG_1000000', 'data_Flower_1000000','data_TB_100000', 'data_SF_100000','data_CC_100000', 'data_CG_100000', 'data_Flower_100000','data_TB_10000000','PenDigits','MNIST'}))
    load([dataName,'.mat'],'fea','gt') % 
    annotation_data = gt;
    data = fea;
    clear fea; clear gt;
else
    load([dataName,'.mat'],'data','annotation_data')
end

if exist('annotation_data','var')&&(min(annotation_data) == 0)
    annotation_data = annotation_data+1;
end
 
ClustN = length(unique(annotation_data(~isnan(annotation_data)))); 

[N,dim]=size(data);
disp(['dataName: ',dataName, ';#objects: ',num2str(N),'; #features: ',num2str(dim),';#Clusters: ',num2str(ClustN)])
