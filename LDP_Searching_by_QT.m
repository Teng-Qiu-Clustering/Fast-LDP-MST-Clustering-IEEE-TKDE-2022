function [neighborIds,supk,rho,r,rs,c] = LDP_Searching_by_QT(data,knnMethod,initial_max_k)
% Written by Teng Qiu (UESTC)
% InPut:
%    data: data set;
%    initial_max_k: threshold of parameter k;
% OutPut:
%    neighborIds: neighborIds(i,:) is the neighbors of node i; note: neighborIds(i,1) = i;
%    supk: value of parameter k; 
%    rho: density vector; r: root label vector; c: cluster label vector; 

%%
disp('determine the inital k nearest neighbors...');
distance_function = 'euclidean';
switch knnMethod
    case 'kd_tree'
        % Note: knnsearch uses the exhaustive search method by default to find the k-nearest neighbors: The number of columns of X is more than 10.
        % knnsearch uses the kd-tree search method by default to find the k-nearest neighbors: The number of columns of X is less than or equal to 10.
        constructed_search_tree = createns(data,'NSMethod','kdtree','Distance',distance_function);
        [neighborIds, knnD] = knnsearch(constructed_search_tree,data,'k',initial_max_k);
    case 'hnsw'
        file_name = 'HnswConstructionforCurrentData'; % file_name can be named in other ways that one like.
        MatlabHnswConstruct(single(data),file_name,distance_function); % for hnsw, its distance function only supports: 'euclidean','l2','cosine','ip';
        [neighborIds, knnD] = MatlabHnswSearch(single(data),initial_max_k,file_name,distance_function);
        neighborIds = double(neighborIds);
        knnD = double(knnD);
end
%% Search natural neighbors
disp('Search natural neighbors...');  
[N,~]=size(data);
nb=zeros(1,N); 
supk=1; count1=0;   
while 1
    for i=1:N
        q=neighborIds(i,supk+1);
        nb(q)=nb(q)+1;
    end
 
    count2=sum(nb==0);
    if (count1==count2) || (supk == initial_max_k - 1)
        break
    else
        supk=supk+1;
        count1=count2;
    end
end
 
neighborIds = neighborIds(:,1:supk+1);
%% Search local density peaks
disp('Search local density peaks...');
dist_sum = sum(knnD,2);
rho=nb./dist_sum'; 
[~,max_ind] = max(rho(neighborIds),[],2);
pr=zeros(N,1);
for i=1:N
    pr(i) = neighborIds(i,max_ind(i));
end
%% find root for each node
disp('find root for each node...');  
rs = [];
r = pr;
for i = 1:N
    if r(i)~=i
        parent=i;
        passed_nodes = i;
        while r(parent)~=parent % search root 
            parent=r(parent);
            passed_nodes = [passed_nodes,parent];
        end
        r(passed_nodes)=parent; % update root label
    else
        rs = [rs,i]; % rs: i.e., root node vector
    end
end
%% initial clustering labeling
disp('initial clustering labeling (cluster label: 1,2,...#roots) ...');

c=zeros(N,1); % c: cluster label vector
c(rs) = 1:length(rs);
c=c(r);
end



