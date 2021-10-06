%% generate data
addpath(genpath(pwd));
disp('generate the test dataset (if this demo is first ran)');
dim = 1;  % dim can be either 1 or 2;
N_array = 2.^(5:22);  
ClustN = 2; 
for i = 1:length(N_array)
    if dim == 1
        fileName = ['One_Dim_uniform_data_',num2str(N_array(i)),'.mat'];
        if ~exist(fileName,'file')
            len = N_array(i)/ClustN;
            data1 = (0:len-1)'+7;  data2 = -(0:len-1)'-7;
            data = [data1; data2]; disp(['size of data:',num2str(length(data))]);
            annotation_data1 = 1*ones(size(data1,1),1);
            annotation_data2 = 2*ones(size(data2,1),1);
            annotation_data = [annotation_data1;annotation_data2];
            %         if N_array(i) < 1000
            %           figure;plot(data,zeros(length(data)),'xb');
            %           box off; hold on; plot([min(data)-5,max(data)+5],[0,0],'k-')
            %         end 
            save(['One_Dim_uniform_data_',num2str(N_array(i)),'.mat'],'data','annotation_data')
        end
    elseif dim == 2
        fileName = ['Two_Dim_uniform_data_',num2str(N_array(i)),'.mat'];
        if ~exist(fileName,'file')
            len = sqrt(N_array(i)/ClustN);
            x = 1:len; y = 1:len;
            [xx yy]=meshgrid(x,y);
            %         figure; plot(xx,yy,'*')
            data1 = [xx(:) yy(:)];  data2 = [xx(:)+len+40 yy(:)];
            data = [data1;data2]; disp(['size of data:',num2str(size(data,1))]);
            annotation_data1 = 1*ones(size(data1,1),1);
            annotation_data2 = 2*ones(size(data2,1),1);
            annotation_data = [annotation_data1;annotation_data2];
            save(['Two_Dim_uniform_data_',num2str(N_array(i)),'.mat'],'data','annotation_data')
        end
    end
end
%% test 
disp('start test...')
method_array = {'FastLDPMST'};
result = []; time = 0; M_array = [];
for m_id = 1:length(method_array)
    method = method_array{m_id};
    time = 0;
    for i = 1:length(N_array)
        disp(N_array(i))         
        if dim == 1 
            load(['One_Dim_uniform_data_',num2str(N_array(i)),'.mat'],'data','annotation_data')
        elseif dim == 2
            load(['Two_Dim_uniform_data_',num2str(N_array(i)),'.mat'],'data','annotation_data')
        end
        
        %% randomly permute the orders of samples (this is used to validate whether the test data is sensitive to the sample order for this test dataset);
        N = size(data,1);
        idx_new = randperm(N);
        data = data(idx_new,:);
        annotation_data = annotation_data(idx_new);
        
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
    end
end
 
%% 
h = figure('position',[282.6 240.2 765.4 521.8]); fonS =12;
markers = {'o','+','x','square','v','^','diamond','pentagram'};
 
subplot(2,3,1);
p = semilogx(N_array,result(:,:,2)*100./N_array,'o-'); 
ylabel('Ratio of root nodes (%)','FontSize',fonS);
for i = 1:length(method_array)
    p(i).Marker = markers{i};
end 
set(gca,'FontSize',fonS);
xlim([min(N_array)-1,max(N_array) + 10])
xlabel('N') 
set(gca,'xtick',10.^(2:2:ceil(max(log10(N_array)))))  
 

subplot(2,3,2);
p = loglog(M_array,result(:,:,3),'o-');
ylabel('Runtime of step 5','FontSize',fonS); 
hold on; plot(M_array,M_array.*log(M_array),'--k');
plot(M_array,M_array.*M_array,'-.k');
text(M_array(10),7*M_array(10).*log(M_array(10)),'M\cdotlog(M)','rotation',15)
text(M_array(10),7*M_array(10).*M_array(10),'M^2','rotation',30)
xlim([min(M_array)-1,max(M_array) + 10])
set(gca,'xtick',10.^(2:2:ceil(max(log10(M_array)))))
set(gca,'ytick',10.^(0:4:ceil(max(log10(M_array.*M_array)))))

set(gca,'FontSize',fonS);
for i = 1:length(method_array)
    p(i).Marker = markers{i};
end
title('  ','FontSize',20)
xlabel('M')
 
subplot(2,3,3);
p = loglog(M_array,result(:,:,4),'o-');
ylabel('Runtime of step 6','FontSize',fonS); 
hold on; 
plot(M_array,M_array.*log(M_array),'--k');
plot(M_array,M_array.*M_array,'-.k');
text(M_array(10),7*M_array(10).*log(M_array(10)),'M\cdotlog(M)','rotation',15)
text(M_array(10),7*M_array(10).*M_array(10),'M^2','rotation',30)
xlim([min(M_array)-1,max(M_array) + 10])
set(gca,'xtick',10.^(2:2:ceil(max(log10(M_array)))))
set(gca,'ytick',10.^(0:4:ceil(max(log10(M_array.*M_array)))))

set(gca,'FontSize',fonS);
for i = 1:length(method_array)
    p(i).Marker = markers{i};
end
title('  ','FontSize',20)
xlabel('M')

subplot(2,3,4);
p = loglog(M_array,result(:,:,5),'o-');
ylabel('Runtime of step 7','FontSize',fonS); 
hold on; 
plot(M_array,M_array.*log(M_array),'--k');
plot(M_array,M_array.*M_array,'-.k');
text(M_array(10),7*M_array(10).*log(M_array(10)),'M\cdotlog(M)','rotation',15)
text(M_array(10),7*M_array(10).*M_array(10),'M^2','rotation',30)
xlim([min(M_array)-1,max(M_array) + 10])
 
set(gca,'xtick',10.^(2:2:ceil(max(log10(M_array)))))
set(gca,'ytick',10.^(0:4:ceil(max(log10(M_array.*M_array)))))

set(gca,'FontSize',fonS);
for i = 1:length(method_array)
    p(i).Marker = markers{i};
end  
xlabel('M')
 
subplot(2,3,5);
p = loglog(N_array,result(:,:,7),'o-'); 
hold on; plot(N_array,N_array.*log(N_array),'--k');
plot(N_array,N_array.*N_array,'-.k'); 
text(N_array(10),7*N_array(10).*log(N_array(10)),'N\cdotlog(N)','rotation',15)
text(N_array(10),7*N_array(10).*N_array(10),'N^2','rotation',30)
xlabel('N','FontSize',fonS);ylabel('Total runtime','FontSize',fonS);
for i = 1:length(method_array)
    p(i).Marker = markers{i};
end
set(gca,'FontSize',fonS);
xlim([min(N_array)-1,max(N_array) + 10])
    set(gca,'xtick',10.^(2:2:ceil(max(log10(N_array))))) 
    set(gca,'ytick',10.^(0:4:ceil(max(log10(N_array.*N_array))))) 
%  ylim([min(result(:,7))-0.1,max(result(:,7))])
xlabel('N') 

subplot(2,3,6);
p = loglog(N_array,result(:,:,6),'o-');
ylabel('ARI','FontSize',fonS);
xlim([min(N_array)-1,max(N_array) + 10]) 
ylim([0 1.1])

set(gca,'FontSize',fonS);
for i = 1:length(method_array)
    p(i).Marker = markers{i};
end
title('  ','FontSize',20)
legend(method_array,'Position', [0.31,0.92,0.43,0.03],'Units', 'normalized','Orientation','horizontal')
xlabel('N')
set(gca,'xtick',10.^(2:2:ceil(max(log10(N_array))))) 

clear data
% 
disp('save figure to "demo3_result.png"');
saveas(gcf,'demo3_result.png')