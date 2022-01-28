function [c,TotalTime] = FastLDPMST(data,nC,MinSize,K,knnMethod)
% Input:
%   data: test dataset (row: samples; column: features);
%   nC: number of clusters;
%   MinSize: cluster size threshold (specifying the minimal cluster size);
%   K: a threshold for parameter k;
%   knnMethod: name of the fast knn method (e.g., 'kd_tree' or 'RP-forest')
% Output:
%   c: cluster label vector; c(i): cluster label of node i; c(i) ranges from 1 to nC;
%   TotalTime: running time (in seconds);

% Written by Teng Qiu (UESTC, Chengdu, China, 2022) for the following paper:
% "Fast LDP-MST: an efficient density-peak-based clustering method for large-size datasets"

% Note 1: for the pseudocodes of all the algorithms in our paper, we use
% "for" loops or nested "for" loops to describe the whole procedure for the
% clarity purpose. However, for Matlab, it is known that vector- and matrix-based code are
% more efficient than loop-based code. According to the above feature of Matlab,
% we have tried our best to avoid loop from occurring during writing the following Matlab codes.
% As a result, many "for" loops in the pseudocodes have been replaced by other equivalent vector- or
% matrix-based codes.

% Note 2: in Matlab, the efficient dataset structures that can be used is
% limited: mainly vectors and matrices. Thus, in many places (e.g.,
% Alg. 2, Alg. 3 and Alg. 4), we have to consider the code in a vector's or matrix's paradigm and make
% some transformation of the speudocodes in our paper to the vector or matrix's form (if the transformation is possible).
% For instance, in the pseudocodes of Alg. 3, we use concatenation operation (e.g., Pass <- [Pass,PN]) to successively store
% all the passed nodes, which is, however, an inefficient operation for Matlab (since Matlab does not
% contains an efficient dataset to handle such operation).
% To solve the above problem, we first use a vector to pre-allocate the maximum space for the passed nodes,
% and then use the assignment operation rather than the concatenation
% operation to fulfill the same goal. In Alg. 2, we also use concatenation
% operation (e.g., rs <- [rs, i]) to find the root nodes one by one. In our
% Matlab, we use an equivalent vector-based code (i.e., rs = find(pr == (1:N)')) to avoid the concatenation
% operation as well.

% Note 3: to get high speed, we have also tried out best to avoid calling Matlab's built-in functions
% (such as union) multiple times in the loop statement (e.g., line 13 in Alg. 4). Alternatively, we
% first store all the possible results (i.e., pairs of p,q) in vectors (at the cost of
% sacrificing some space), and then use the "unique" function to extract
% all the unique pairs. And the results are not stored in a binary set Xi
% with each element being a pair of element: (p,q), representing the row and column indexes of a non-zero
% element in a matrix; instead, they are stored in a matrix with two columns
% (each row of the matrix represents a pair of elements).

%% error control of the input
[N,dim]=size(data);
if nargin < 2
    error('At Least Two Inputs Are Needed');
end
if nargin < 3
    MinSize = 0.01*N; % empirical parameter setting.
end
if nargin < 4
    K = ceil(log2(N)); % empirical parameter setting.
end
if nargin < 5
    if dim <= 10
        knnMethod = 'kd_tree'; % for low-dimensional datasets
    else
        knnMethod = 'RP-forest'; % for high-dimensional datasets.
    end
end

%% Steps 1 to 4
TotalTime_start = tic;
% the following sub-function (involving steps 1 to 4) obtains the initial clustering result
[neighborIds,f,rs,c] = LDP(data,knnMethod,K); % steps 1 to 4
%  neighborIds: a matrix; neighborIds(i,:) stores the k nearest neighbors of node i; neighborIds(i,1) = i, for all i;
%  f: density vector; f(i) stores the density of node i;
%  rs: a vector storing the indexes (in ascending order) of all the root nodes;
%  c: cluster label vector; c(i) stores the cluster label of node i;

%% Steps 5 to 8
Rsamples = data(rs,:); % Rsamples stores the samples corresponding to the root nodes;
clear data; % save some space when processing very large datasets
M = length(rs); % number of root nodes = the number of initially obtained clusters (i.e.,
% the number of supernodes)
if M > nC
    %% Step 5: compute adjacency (or weight) matrix of a weighted graph
    % step 5 contains the following two sub-steps (step 5.1 and step 5.2).
    step5_time_clusterDistance_start = tic;
    % step 5.1
    disp('step 5.1: determine multiple cluster labels of each node...')
    % The lines 4 and 5 of the pseudocodes of Alg. 4 show a nested
    % "for" loop, which should be avoided in Matlab (since Matlab is more efficient on
    % processing vectors than loops). In the following, we first determine the reverse neighbor i of each
    % node j (not that in lines 4 and 5 of Alg. 4, j is a neighbor of node
    % i; conversely, i is a reverse neighbor of node j). The goal of lines
    % 3 to 5 is to union the cluster label of all the reverse neighbors of
    % nodes j (note that the cluster labels of some reverse neighbors could
    % be the same).
    k = size(neighborIds,2);
    Temp = [neighborIds(:),c(repmat((1:N)',k,1))]; % Each row contains two elements, storing a node j (1st element) and the cluster label of its reverse neighbor (2nd element).
    clear neighborIds; % free up some space.
    
    % Note that in Alg. 4, ML is a data structure containing N sets (i.e., ML(j),j = 1 to N). 
    % In the following, we will not try to store the information by sets.
    % Insteads, the vectors will be used for storing the information in ML.
    % Specifically, we will use "unique" to remove the repeated row;
    % consequently, the "union" operation (in line 5 of Alg. 4) is not needed.
    % And the output matrix Temp2 (= [I J]) stores all the information associated
    % with ML, where J(t) stores the reverse cluster label of node I(t), t>=1.
    
    Temp2 = unique(Temp,'rows'); % C = unique(A,'rows') for the matrix A returns the unique rows of A. The rows of the matrix C will be in sorted order.
    clear Temp;
    I = Temp2(:,1); % Note: the values in I have been sorted in ascending order in "unique" function;
    J = Temp2(:,2); % Note: the values in J have been sorted in ascending order in "unique" function;
    clear Temp2;
    
    % Then, we determine the size of ML(j) indirectly based on the repeatedness of the elements j in vector I.
    % Note that the values in I have been sorted in ascending order
    % by "unique" function. In order words, the same values are 
    % adjacent in vector I (e.g., [1,1,1,2,2,3,3]). 
    % Thus, to determine the size of ML(2), we need to determine (from left to right) the 
    % first index (i.e., 4) and last index (i.e., 5) of 2 in vector I, then the size of ML(2) is
    % 5 - 4 + 1 = 2. 
    
    start_idx = zeros(N,1);
    end_idx = zeros(N,1);
    ML_size = zeros(N,1);
    j =  1;
    start_idx(j) = 1;
    for t = 2:length(I)
        if I(t) ~= I(t-1)
            end_idx(j) = t - 1;
            ML_size(j) = end_idx(j) - start_idx(j) + 1;
            j = j + 1;
            start_idx(j) = t;
        end
    end
    end_idx(j) = length(I); % in this line, j should be N.  
    ML_size(j) = end_idx(j) - start_idx(j) + 1;
    
    % step 5.2
    disp(['step 5.2: compute a sparse weight matrix of size:',num2str(M),' x ',num2str(M),'...']);
    Total_Num = sum(ML_size.*(ML_size - 1))/2; % Total_Num: the total number of pairs of (p,q) that requires to be considered in line 9 of Alg. 4
    Pairs = zeros(Total_Num,2); % Pre-allocate the space. Each row of Pairs will be used to store a pair of elements (p,q) in any set ML(j).
    Pairs_f = zeros(Total_Num,1); % Pre-allocate the space. Pairs_f(t) will be used to store the density information f(j) corresponding to a pair of (p,q)
    
    % in the following, we first store all the pairs of elements (p,q) in
    % each set ML(j) in matrix Pairs (based on the cluster label vector J) and store the corresponding densities in
    % vector Pairs_f.
    
    tt = 1;
    for j = 1:N
        if ML_size(j)~= 1
            ML_j = J(start_idx(j):end_idx(j)); % note: the values in J has been sorted in ascending order (according to the unique function);
            for s = 1:ML_size(j)-1
                for t = s+1:ML_size(j)
                    % first store in Pairs(tt,1) and Pari(tt,2) all
                    % different elements p,q in ML(j),respectively. Note:
                    % there could exist repeated rows in Pairs, which will be removed later.
                    Pairs(tt,1) = ML_j(s);
                    Pairs(tt,2) = ML_j(t);
                    Pairs_f(tt) = f(j); % store the corresponding density f(j);
                    tt = tt + 1;
                end
            end
        end
    end
    clear I J;
        
    % Then, we compute the non-zero elements of matrix F, Count, and D, as well as
    % the set Xi (in line 13) of Alg. 4 in the following way.
    % First, we get all the different pairs of (p,q) from matrix Pairs (using "unique"
    % function in Matlab); the output matrix DiffPairs stores the
    % coordinates of the non-zero elements in the weight matrix D (note that
    % matrix DiffPairs can be viewed as the matrix form of the set Xi
    % in line 13 of Alg. 4, i.e, each row of DiffPairs can be viewed as an element
    % in set Xi); Then, we accumulate the densities of those
    % repeated pairs stored in Pairs so as to get non-zero elements of
    % matrix F in Alg. 4 and we also count the number of those
    % repeated pairs so as to get the non-zero elements of matrix Count in
    % Alg.4. Note that the above goals are fulfilled based on another output of "unique" function,
    % i.e., vector IC in the following. And note that we store the non-zero elements
    % of matrix F and matrix Count in vectors F_v and Count_v, respectively.
    
    [DiffPairs,~,IC] = unique(Pairs,'rows'); % get all the different rows of matrix Pairs, each representing the row and column index of a non-zero element in weight matrix D;
    % note: Pairs = DiffPairs(IC,:);
    P = DiffPairs(:,1); % row index of the non-zero elements of D;
    Q = DiffPairs(:,2); % column index of the non-zero elements of D;
    % P stores all the values of p, and Q stores all the corresponding
    % values of q in Alg. 4.
    
    Num_nonZeroE = length(P); % Num_nonZeroE: the number of all the non-zero elements in weight matrix D.
    F_v = zeros(Num_nonZeroE,1); % Pre-allocate the space.
    Count_v = zeros(Num_nonZeroE,1); % Pre-allocate the space.
    for t = 1:length(IC)
        % note: Pairs = DiffPairs(IC,:);
        % IC(t) specifies which row of DiffParis is the same as the t-th
        % row of Pairs
        % Accordingly, the the following two lines of code can be viewed as the 1-dimensional
        % indexing version of the corresponding 2-dimensional indexing
        % version in lines 10 and 11 of Alg. 4. 
        F_v(IC(t)) = F_v(IC(t)) + Pairs_f(t); % vector F_v stores all the non-zero elements of matrix F in Alg.4
        Count_v(IC(t)) = Count_v(IC(t)) + 1; % vector Count_v stores all the non-zero elements of matrix Count in Alg.4
    end
    
    D_v = zeros(Num_nonZeroE,1); % D_v will store the values of all the non-zero values of weight matrix D in Alg.4
    for t = 1:Num_nonZeroE
        d = sqrt(sum((Rsamples(P(t),:)-Rsamples(Q(t),:)).^2)); % numerator of Eq. 6
        D_v(t) = d/F_v(t)/Count_v(t);  % compute a non-zero element of weight matrix D according to Eq. 6 (i.e., the non-zero distances between clusters or supernodes);
        % Note: D_v, F_v, Count_v are the three vectors storing all the
        % non-zero elements of matrices D, F, and Count in Alg. 4, respectively.
    end
    step5_time_clusterDistance = toc(step5_time_clusterDistance_start);
    disp(['  Time Cost on Step 5: ',num2str(step5_time_clusterDistance),'s']);
    
    % Note: an easier and more straight-forward way of determining matrices F, Count, and D of Alg.4 is by virtue of sparse function in Matlab.
    % Specifically, define matrices F, Count, and D as sparse matrices and then do the assignment and query operation
    % based on the pseudocode of Alg.4. However, our experiments show that
    % such easier sparse-matrix-based programming
    % is slightly slower than the above non-sparse-function-based programming in Matlab; for instance,
    % on a dataset with millions of samples, the non-sparse-function-based code spent around 10 second, in contrast to
    % around 30 seconds for the sparse-function-based code (omitted here).
    % This is due to the following limitation for the sparse matrix in Matlab.
    % In Matlab, assigning a value to an element of a sparse matrix is not
    % as efficient as querying an element of a sparse matrix; in Alg. 4, there involves such code (line 10):
    % F(p,q) <- F(p,q) + f(j), which involves the assignment operations on F. Since line 10
    % of Alg. 4 exists in a nested loop, this will amplify the above limitation of sparse matrix.
    
    %% Step 6: Construct MSF
    step6_time_on_MSF_start = tic;
    disp('step 6: Construct MSF ...');
    G = graph(P,Q,D_v,M); % construct the sparse weighted graph based on the non-zero elements of the weight matrix D
    [T,pred] = minspantree(G,'Type','forest','Method','sparse'); % 'sparse' means that the Kruskal's (rather than Prim's) method is used
    ST = sparse(T.Edges.EndNodes(:,1),T.Edges.EndNodes(:,2),T.Edges.Weight,M,M);
    if any(isnan(pred))
        % % Note: there is problem for the output 'pred' when UG is
        % unconnected for Matlab 2013; this problem has been modified by Matlab 2016 and higher versions
        error('a higher version of Matlab is needed (e.g., Matlab 2016 and later version than that)');
    end
    % transform MST to in-tree-based MSF
    disp('  transform MSF to in-tree-based MSF...')
    idx_roots = find(pred == 0);
    pred(idx_roots) = idx_roots;
    pr1 = pred;
    EW = zeros(1,M); % Pre-allocate the space for edge weight vector of an in-tree-based MSF.
    for i = 1:M
        EW(i) = max(ST(i,pr1(i)),ST(pr1(i),i));  % ST is not a symmetric matrix, and thus max is used to get the non-zero element
    end
    step6_time_on_MSF = toc(step6_time_on_MSF_start);
    disp(['  Time Cost on Step 6: ',num2str(step6_time_on_MSF),'s']);
    
    %% Step 7: Cut the forest
    disp('Step 7: Cut the tree...')
    step7_time_on_cutting_start = tic;
    [pr2,rs2] = EdgeCutting(pr1,EW,c,MinSize,nC); % Supplementary Alg. S2; pr1: pr' in Alg. S1; pr2: pr'' in Alg. S2.
    step7_time_on_cutting = toc(step7_time_on_cutting_start);
    disp(['  Time Cost on Step 7: ',num2str(step7_time_on_cutting),'s'])
    
    %% Step 8: update cluster labels of all the nodes
    disp('Step 8: final cluster assignment...');
    step8 = tic;
    sp_c = ClusterLabeling(pr2,rs2); % first, get the cluster labels of all supernodes (line 1 of supplementary Alg. S3); sp_c(i): the cluster label of the supernode i;
    c = sp_c(c); % then, get the cluster label of all the (normal) nodes (line 2 of supplementary Alg. S3).
    step8_time = toc(step8);
    disp(['  Time Cost on Step 8: ',num2str(step8_time),'s'])
    
else
    % this case is explained in the 2nd footnote of the paper (Section 4.9).
    warning('  in this case, steps 5 to 8 are not needed')
end


TotalTime=toc(TotalTime_start);
disp(['Time Cost on whole method: ',num2str(TotalTime),'s']);
disp('End');


