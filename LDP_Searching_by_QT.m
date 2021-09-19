function [neighborIds,supk,rho,r,rs,c,pr,W] = LDP_Searching_by_QT(data,knnMethod,initial_max_k,dataName)
% Written by Teng Qiu (UESTC)
% InPut:
%    data: data set;
%    initial_max_k: threshold of parameter k;
% OutPut:
%    neighborIds: neighborIds(i,:) is the neighbors of node i; note: neighborIds(i,1) = i;
%    supk: value of parameter k; 
%    rho: density vector; r: root label vector; c: cluster label vector; 

%%
[N,~]=size(data);
disp('determine the inital k nearest neighbors...');
distance_function = 'euclidean';
disp(['knnMethod: ', knnMethod])
time_on_fastKNN_start = tic;
switch knnMethod
    case 'kd_tree'
        % Note 1: knnsearch uses the exhaustive search method by default to find the k-nearest neighbors: The number of columns of X is more than 10.
        % knnsearch uses the kd-tree search method by default to find the k-nearest neighbors: The number of columns of X is less than or equal to 10.
        % Note 2: the initial neighbor is the searched node itself
        constructed_search_tree = createns(data,'NSMethod','kdtree','Distance',distance_function);
        [neighborIds, knnD] = knnsearch(constructed_search_tree,data,'k',initial_max_k);
    case 'hnsw'
        file_name = ['HnswConstructionforCurrentData',dataName]; % file_name can be named in other ways that one like.
        MatlabHnswConstruct(single(data),file_name,distance_function); % for hnsw, its distance function only supports: 'euclidean','l2','cosine','ip';
        [neighborIds, knnD] = MatlabHnswSearch(single(data),initial_max_k,file_name,distance_function);
        neighborIds = double(neighborIds);
        knnD = double(knnD);
        % Note: we find that hnsw cannot guarantee that the 1st nearest neighbor of a
        % node is itself.The following loop is used to remedy this bug in hnsw.  
        for i=1:N
            if neighborIds(i,1) ~= i
                neighborIds(i,2:end) = neighborIds(i,1:end-1);
                neighborIds(i,1) = i;
                knnD(i,2:end) = knnD(i,1:end-1);
                knnD(i,1) = 0;
            end
        end  
    case 'NNDescent' % initialized by RP trees
        [neighborIds, knnD] = KnnFind.Approximate(double(data),[],'K',initial_max_k);
    case 'RP_kdTree'
        Reduced_Dim = 10;
        T=rand(size(data,2),Reduced_Dim);
        %   Second type of RP
        T(T<(1/3))=-sqrt(3);
        T(T>=(2/3))= sqrt(3);
        T(T<(2/3)&T>=(1/3))=sqrt(3);
        data=data*T/sqrt(Reduced_Dim);
        
         constructed_search_tree = createns(data,'NSMethod','kdtree','Distance',distance_function);
        [neighborIds, knnD] = knnsearch(constructed_search_tree,data,'k',initial_max_k);
end
time_on_fastKNN = toc(time_on_fastKNN_start);
disp(['Time Cost on fast knn: ',num2str(time_on_fastKNN),'s']);
%% Search natural neighbors
disp('Search natural neighbors...');  

nb=zeros(1,N); 
supk=1; count1=0;   
while 1
    for i=1:N
        q=neighborIds(i,supk+1);
        nb(q)=nb(q)+1;
    end
 
    count2=sum(nb==0);
    if (count1==count2) || (supk+1 == initial_max_k)
        break
    else
        supk=supk+1;
        count1=count2;
    end
end
 
neighborIds = neighborIds(:,1:supk+1);
 
knnD = knnD(:,1:supk+1);
% knnD = knnD(:,1:min(max(nb),initial_max_k));
%% Search local density peaks
disp('Search local density peaks...');
dist_sum = sum(knnD,2);
rho=nb./(dist_sum'+eps); 
 
[~,max_ind] = max(rho(neighborIds),[],2);
pr =zeros(N,1);
W = zeros(N,1);
for i=1:N
    pr(i) = neighborIds(i,max_ind(i));
    W(i) = knnD(i,max_ind(i));
end
%% find root for each node
disp('find root for each node...');  
rs = find(pr == (1:N)'); % rs: i.e., root node vector
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
    end
end
%% initial clustering labeling
disp('initial clustering labeling (cluster label: 1,2,...#roots) ...');

c=zeros(N,1); % c: cluster label vector
c(rs) = 1:length(rs);
c=c(r);
end



