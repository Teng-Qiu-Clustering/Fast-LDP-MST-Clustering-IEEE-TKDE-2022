
% Written by Teng Qiu
% A: the data set

function [neighborIds,supk,nb,rho,local_core,cores,cl,ini_cluster_number] = LDP_Searching_by_QT(A,knnMethod,initial_max_k)
[N,dim]=size(A);
disp('determine the inital k nearest neighbors...');

% % Note: knnsearch uses the exhaustive search method by default to find the k-nearest neighbors: The number of columns of X is more than 10.
% % knnsearch uses the kd-tree search method by default to find the k-nearest neighbors: The number of columns of X is less than or equal to 10.
 
    
distance_function = 'euclidean';
 switch knnMethod
    case 'kd_tree'
        constructed_search_tree = createns(A,'NSMethod','kdtree','Distance',distance_function);
        [neighborIds, knnD] = knnsearch(constructed_search_tree,A,'k',initial_max_k);
    case 'hnsw' 
        file_name = 'HnswConstructionforCurrentData'; % file_name can be named in other ways that one like.
        MatlabHnswConstruct(single(A),file_name,distance_function); % for hnsw, its distance function only supports: 'euclidean','l2','cosine','ip';
        [neighborIds, knnD] = MatlabHnswSearch(single(A),initial_max_k,file_name,distance_function);
        neighborIds = double(neighborIds);
        knnD = double(knnD);
end

%初始化基本数据
supk=1;
flag=0;
nb=zeros(1,N);  %自然邻居个数
%NNN=zeros(N,N); %各点的自然邻居集
count=0;        %自然最近邻数为零的数据量连续相同的次数
count1=0;       %前一次自然最近邻数为零的数据量 

disp('Search natural neighbors...');
while flag==0
    for i=1:N
        q=neighborIds(i,supk+1);
        nb(q)=nb(q)+1;
    end
    supk=supk+1;
    count2=sum(nb==0);
    %计算nb(i)=0的点的数量连续不变化的次数
    if count1==count2
        count=count+1;
    else
        count=1;
    end
    if count2==0 || (supk>2 && count>=2) || supk == initial_max_k   %邻居搜索终止条件
        flag=1;
    end
    count1=count2;
 
end


%计算自然最近邻的各种特征量
supk=supk-1;               %最终K值，也是自然最近邻居的平均数
max_nb=max(nb);         %自然邻居的最大数目

% if initial_max_k < max_nb + 1
%     initial_max_k = max_nb+1;
%     switch knnMethod
%         case 'kd_tree'
%             [neighborIds, knnD] = knnsearch(constructed_search_tree,A,'k',initial_max_k);
%         case 'hnsw'
%             [neighborIds, knnD] = MatlabHnswSearch(single(A),initial_max_k,file_name,distance_function);
%             neighborIds = double(neighborIds);
%             knnD = double(knnD);
%     end
% end
 

disp(['max_nb = ',num2str(max_nb),' supk = ',num2str(supk)]);

disp('Search local density peaks...');
dist_sum = sum(knnD,2);
rho=nb./dist_sum';

[max_density_unused,max_ind] = max(rho(neighborIds(:,1:supk+1)),[],2);
local_core=zeros(N,1);
for i=1:N
    local_core(i) = neighborIds(i,max_ind(i));
end

% find root for each node
% disp('find root for each node...');
% local_core_update = local_core(local_core);
% while any(local_core~=local_core_update)
%     local_core=local_core_update;
%     local_core_update = local_core(local_core);
% end

disp('find root for each node...');
local_core_update = local_core; 
while 1
    for i = 1:N
        local_core(i) = local_core(local_core(i));
    end
    if all(local_core==local_core_update)
        break
    else 
        local_core_update = local_core; 
    end
end


disp('initial clustering labeling (cluster label: 1,2,...#roots) ...');
cores = find(local_core' == 1:N); % cores: i.e., root nodes
ini_cluster_number = length(cores);
cl=zeros(N,1);
cl(cores) = 1:ini_cluster_number;
cl=cl(local_core);

neighborIds = neighborIds(:,1:supk+1);
end



