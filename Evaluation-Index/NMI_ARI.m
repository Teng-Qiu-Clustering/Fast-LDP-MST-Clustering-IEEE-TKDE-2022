function [NMI_max,ARI,NMI_sqrt,AMI,RI,HI,AVI,EMI]= NMI_ARI(Label,annotation_data)

% Proprcessing: make the labels started from 1 (This is not necessary For
% NMI and ARI; Note: this could be needed for computing other evaluation
% metric (e.g., the Accuracy metric used in some literature)

% temp=zeros(size(Label));u=unique(Label);
% for tt=1:length(u),temp(Label==u(tt))=tt; end
% Label = temp;

[NMI_max,AMI,AVI,EMI,NMI_sqrt] = ANMI_analytical_11(annotation_data,Label); 
% Note that AMI could be NaN and when other indexes are far from 1, AMI could be larger than 1; 

% external clustering index (%ARI: adjusted Rand index, RI: the unadjusted Rand index, HI: "Hubert's" index.)
[ARI,RI,HI]=valid_RandIndex(annotation_data,Label);
%internal index (DB:Davies-Bouldin, CH:Calinski-Harabasz
%     [DB, CH, KL, Ha, ST]=valid_internal_deviation(data,Label,1);
%     DI=dunns(length(unique(Label)),D,Label);


