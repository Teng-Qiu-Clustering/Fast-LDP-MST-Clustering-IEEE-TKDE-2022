function c = ClusterLabeling(pr,rs)
% Written by Teng Qiu (UESTC), Chengdu, China, 2022
% Input:
%    pr: parent node vector;
%    rs: root node vector;
% Output:
%     c: cluster label vector;

% first, determine the root label of each node;
r = pr; % initialize the root label vector r (we will update it in the following loop)
N = length(pr);
passed_nodes = zeros(1,N); % Pre-allocate the maximum space for the passed nodes, since there are maximally N passed nodes.
% note: the above pre-allocation aims to avoid using concatenation operation 
% (shown in the pseudocode of Alg.3) to store all the passed nodes;
% since concatenation operation is inefficient for Matlab.
for i = 1:N
    if r(i)~=i
        parent=i;
        t = 1; % a variable to count the number of passed nodes
        passed_nodes(t) = i;
        while r(parent)~=parent % search root
            parent=r(parent);
            t = t + 1;
            passed_nodes(t) = parent;
        end
        r(passed_nodes(1:t))=parent; % update root label of all the passed nodes (note that here "parent" stores the index of the reached root).
    end
end

% then determine the clustering labeling (cluster label: 1,2,...#roots)
c=zeros(N,1); % Preallocate the space for the cluster label vector; c(i): cluster label of node i;
c(rs) = 1:length(rs); % first, assign cluster labels to the root nodes;
c=c(r); % then, assign cluster labels to non-root nodes based on the root nodes they have reached