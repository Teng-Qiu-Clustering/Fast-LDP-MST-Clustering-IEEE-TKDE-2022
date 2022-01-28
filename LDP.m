function [neighborIds,f,rs,c] = LDP(data,knnMethod,K)
% Written by Teng Qiu (UESTC), Chengdu, China, 2022
% Implementation of steps 1 to 4. 
% Input:
%    data: data set (each row denotes a sample; each column denotes a dimension);
%    K: threshold of parameter k;
% Output:
%    neighborIds: a matrix, where neighborIds(i,:) stores the indexes of all the neighbors (in ascending order of the distance) of node i; note: neighborIds(i,1) = i;
%    f: density vector;
%    r: root label vector;
%    rs: a vector storing the indexes (in ascending order) of all the root nodes
%    c: cluster label vector;

%% Step 1: detemine the natural neighbors
step1 = tic;
disp('step 1.1: determine the inital k nearest neighbors...');
[N,~]=size(data);
distance_function = 'euclidean';
disp(['  knnMethod: ', knnMethod])
switch knnMethod
    case 'kd_tree' 
        constructed_search_tree = createns(data,'NSMethod','kdtree','Distance',distance_function);
        [neighborIds, knnD] = knnsearch(constructed_search_tree,data,'k',K);
    case 'RP-forest' % initialized by RP-forest and then modifed by NNDescent
        if ~isa(data,'double'), data = double(data); end % for KnnFind.Approximate, the input data should be double-type
        [neighborIds, knnD] = KnnFind.Approximate(data,[],'K',K); 
end
% Note 1: the 1st nearest neighbor of each node i is node i itself,
% that is, neighborIds(i,1) = i and knnD(i,1) = 0. 
% The following if statement further guarantee the above premise of the first neighbor.
if ~all(neighborIds(:,1) == (1:N)')
    error('for each node i, its 1st nearest node should be itself');
end

disp('step 1.2: Search natural neighbors...');
nb=zeros(1,N);
k=1; count1=0;
while 1
    for i=1:N
        q=neighborIds(i,k+1);
        nb(q)=nb(q)+1;
    end
    count2=sum(nb==0);
    if (count1==count2) || (k+1 == K)
        k=k+1; 
        break
    else
        k=k+1; 
        count1=count2;
    end
end

neighborIds = neighborIds(:,1:k);  
knnD = knnD(:,1:k);

time_on_step1 = toc(step1);
disp(['  Time Cost on Step 1: ',num2str(time_on_step1),'s']);

% Step 2: compute density
step2 = tic;
disp('step 2: compute density...')
dist_sum = sum(knnD,2); % Accumulate each row of matrix knnD; 
f=nb./(dist_sum'+eps); % f(i): density of node i; eps is added to avoid the denominator being 0.
time_on_step2 = toc(step2);
disp(['  Time Cost on Step 2: ',num2str(time_on_step2),'s']);

% Step 3: determine the parent node
step3 = tic;
disp('step 3: determine the parent node...')
[~,max_ind] = max(f(neighborIds),[],2); 
% note: to determine the parent node (i.e. pr(i)) of each node i, 
% we need to find the node with the maximal density among 
% the k nearest nodes stored in the i-th row of neighborIds.
% However, the i-th element, max_ind(i), in the above output 
% does not directly store pr(i); instead,
% max_ind(i) stores the index of pr(i) in neighborIds(i,:). Thus,
% to get pr(i), we need to further query neighborIds(i,max_ind(i)) as follow.
pr =zeros(N,1); % Pre-allocate the space for the parent node vector. 
for i=1:N
    pr(i) = neighborIds(i,max_ind(i)); % pr(i): parent node of node i;
end
rs = find(pr == (1:N)'); % rs: a vector storing the indexes (in ascending order) of all the root nodes
% note: function "find" returns the indexes of the vector pr where pr(i) == i. 
% The above vector-based code for determining rs is equivalent 
% to the following lines of Matlab code: 
% rs = [];
% for i = 1:N
%     if pr(i) == i
%         rs = [rs,i];
%     end
% end

% Note: since Matlab is more efficient in processing vectors than the loop, in 
% our Matlab code, we chose the above vector-based code instead of the loop-based code. 

time_on_step3 = toc(step3);
disp(['  Time Cost on Step 3: ',num2str(time_on_step3),'s']);

% Step 4: get initial clustering assignments
step4 = tic;
disp('step 4: get initial clustering assignments...');
c = ClusterLabeling(pr,rs);
time_on_step4 = toc(step4);
disp(['  Time Cost on Step 4: ',num2str(time_on_step4),'s']);
