function [c,TotalTime,Ini_clusterNum,supk,time_clusterDistance,time_on_MSF,time_on_cutting] = FastLDPMST(data,clu_num,minsize,initial_max_k,knnMethod,dataName)
%Input:
%  data: data set (nrow: samples; col: features);
%  clu_num: number of clusters;
%Output:
%  c: cluster label vector;
%  TotalTime: running time (seconds);

% Written by Teng Qiu (UESTC)

%%
[N,dim]=size(data);
if nargin < 2
    error('at least two inputs are needed');
end
if nargin < 3
    minsize = 0.01*N;
end
if nargin < 4
    initial_max_k = ceil(log2(N));
end
if nargin < 5
    if dim <= 10
        knnMethod = 'kd_tree';       disp('kd-tree exact fast knn searching technique is used for the dataset with Dimension lower than 20')
    else
        knnMethod = 'NNDescent';       disp('hnsw (L2 distance) approximate fast knn searching technique is used for the dataset with Dimension larger than 20')
    end
end
if nargin < 6
    dataName = '__';
end

TotalTime_start = tic;
%% initial clustering by LDP
time_initial_clustering_start = tic;
[knnIndex,supk,rho,r,rs,c,pr,W] = LDP_Searching_by_QT(data,knnMethod,initial_max_k,dataName);
time_initial_clustering = toc(time_initial_clustering_start);
disp(['Time Cost on initial clustering (including fast knn): ',num2str(time_initial_clustering),'s']);

A_cores = data(rs,:); 
clear data
M = length(rs);
Ini_clusterNum = M;
if M > clu_num
    disp('determine the size of each component of the initial forest...');
    nc=zeros(1,M); 
    for i = 1:N
        nc(c(i)) = nc(c(i)) + 1;
    end   
    
    %% compute shared sparse distance matrix (containing the following two sub-steps)
    time_clusterDistance_start = tic;
    disp('first: determine multiple labels of each node...')
    K = supk + 1;
    % method 1 (the following two lines):
    % ReverseNeighbor_Labels = sparse(knnIndex(:),c(repmat((1:N)',K,1)),ones(N*K,1),N,M);
    % [I,J] = find(ReverseNeighbor_Labels); % J: cluster label vector, where J(t) is a label of node I(t)
    % method 2 (the following three lines); Methods 1 and 2 are equvalent in both effectiveness and efficiency
    Temp = [knnIndex(:),c(repmat((1:N)',K,1))];
    clear knnIndex;
    Temp2 = unique(Temp,'rows'); % C = unique(A,'rows') for the matrix A returns the unique rows of A.The rows of the matrix C will be in sorted order.
    clear Temp;
    I = Temp2(:,1); J = Temp2(:,2);
    clear Temp2;
    % Note: the values in I have been sorted in ascending order in unique function;
   
    FirstAppear_idx = zeros(N,1);
    FirstAppear_value = zeros(N,1);
    L_size = zeros(N,1);
    i =  1;  
    FirstAppear_idx(i) = 1;
    FirstAppear_value(i) = I(1);  %% note: in most cases (i.e., if knnIndex(j,1)= j,for all j = 1 to N), FirstAppear_value(j) = j, for all j = 1 to N;
    for t = 2:length(I)
        if I(t) ~= I(t-1)           
            i = i + 1;
            FirstAppear_value(i) = I(t);
            FirstAppear_idx(i) = t;
            L_size(FirstAppear_value(i-1)) = FirstAppear_idx(i) - FirstAppear_idx(i-1);
        end 
    end   
    L_size(FirstAppear_value(i)) = length(I) + 1 - FirstAppear_idx(i); % note: +1 is necessary here
 
    % note: in most cases (i.e., if knnIndex(j,1)= j,for all j = 1 to N), i = N
    
    disp(['Then: compute shared sparse distance matrix of size:',num2str(M),' x ',num2str(M),'...']);
    
    Total_Num = sum(L_size.*(L_size - 1))/2;
    A = zeros(Total_Num,1);
    B = zeros(Total_Num,1);
    C = zeros(Total_Num,1); 
%     D = ones(Total_Num,1);
   
    tt = 1;
    for m = 1:i
        j = FirstAppear_value(m);
        if L_size(j)~= 1
            L_j = J(FirstAppear_idx(m):(FirstAppear_idx(m)+L_size(j)-1)); 
%             L_j = sort(L_j,'descend'); % this has been fulfilled (sorted in ascending order) in by unique function (it makes no difference either in ascending or descending order);
            for s = 1:L_size(j)-1
                for t = s+1:L_size(j) 
                    A(tt) = L_j(s);
                    B(tt) = L_j(t);
                    C(tt) = rho(j); 
                    tt = tt + 1; 
                end
            end
        end
    end
    clear I J;
    % method 1
%     rho_matrix= sparse(A,B,C,M,M);
%     count_matrix = sparse(A,B,D,M,M);
%     clear A B C D;
%     [I,J,V] = find(rho_matrix);
%     [~,~,V2] = find(count_matrix);
%     for t = 1:length(I)
%         d=sqrt(sum((A_cores(I(t),:)-A_cores(J(t),:)).^2));
%         V(t) = d/V(t)/V2(t); % slightly faster than V(t) = d/V(t)/count_matrix(I(t),J(t));
%     end

% method 2
    Y = [A,B];
    [CC,~,IC] = unique(Y,'rows'); % Y = CC(IC,:);
    L = size(CC,1);
    V = zeros(L,1);
    V2 = zeros(L,1);
    for i = 1:length(IC)
        V(IC(i)) = V(IC(i)) + C(i);
        V2(IC(i)) = V2(IC(i)) + 1;
    end
    for t = 1:L
        d=sqrt(sum((A_cores(CC(t,1),:)-A_cores(CC(t,2),:)).^2));
        V(t) = d/V(t)/V2(t); % slightly faster than V(t) = d/V(t)/count_matrix(I(t),J(t));
    end
    
    
    G = graph(CC(:,1),CC(:,2),V,M);
    time_clusterDistance = toc(time_clusterDistance_start);
    disp(['Time Cost on compute cluster distance: ',num2str(time_clusterDistance),'s']);
    %% Construct MSF
    time_on_MSF_start = tic;
    disp('Construct MSF on root nodes (rs) only ...');   
    [T,pred] = minspantree(G,'Type','forest','Method','sparse');
    ST = sparse(T.Edges.EndNodes(:,1),T.Edges.EndNodes(:,2),T.Edges.Weight,M,M);
    if any(isnan(pred))
        % % ****** Note:  there is problem for the output 'pred' when UG is unconnected for matlab 2013; this problem has been modified by matlab 2016 and above
        error('a higher version of matlab is needed (e.g., matlab 2016 and later version than that)');
    end
    % transform MST to in-tree-based MSF
    disp('transform MSF to in-tree-based MSF...')
    idx_roots = find(pred == 0);
    pred(idx_roots) = idx_roots;
    I_base = pred;
    W = zeros(1,M);
    for i = 1:M
        W(i) = max(ST(i,I_base(i)),ST(I_base(i),i));
    end
    time_on_MSF = toc(time_on_MSF_start);
    disp(['Time Cost on building in-tree-based MSF: ',num2str(time_on_MSF),'s']);
    %% Cut the tree
    disp('Cut the tree...')
    time_on_cutting_start = tic;
    cores_cl = Cut_MST_QT_v2(I_base,W,nc,minsize,clu_num);
    time_on_cutting = toc(time_on_cutting_start);
    disp(['Time Cost on cutting: ',num2str(time_on_cutting),'s'])
    %% final cluster assignment
    disp('final cluster assignment...');
    c=zeros(N,1);
    c(rs)=cores_cl;
    c=c(r);
else
    disp('the initially generated number of clusters is too small, and thus there is no need to merge the clusters')
    c = Cut_MST_QT_v2(pr',W,ones(N,1),minsize,clu_num);
    time_clusterDistance = nan;
    time_on_MSF = nan;
    time_on_cutting = nan;
    disp('there is no need to merge !!!!!!');
end


TotalTime=toc(TotalTime_start);
disp(['Time Cost on whole method: ',num2str(TotalTime),'s']);
disp('End');
end


