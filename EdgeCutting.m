function [pr2,rs2] = EdgeCutting(pr1,EW,c,MinSize,nC)
% Written by Teng Qiu (UESTC), Chengdu, China, 2022
% % implementation of Supplementary Alg. S2 
% Input: 
%     pr1: i.e., pr' in Alg. S1, which is the parent supernode vector for the supernodes; pr1(i) denotes the parent supernode of supernode i;
%     EW: i.e., ew' in Alg. S2, denoting the edge weight vector; EW(i) denotes the edge between supernode i and its parent supernode pr1(i);
%     c: initial cluster label vector (obtained in step 4). 
%     MinSize: cluster size threshold for the cutting (specifying the minimal cluster size);
%     nC: expected number of clusters;
% Output: 
%     pr2: i.e., pr'' in Alg. S2, denoting the updated parent supernode vector after edge cutting. 
%     rs2: i.e., rs'' in Alg. S2, denoting the root supernode vector. 

% Note: the input NW is the result (i.e.,initial weight of each supernode) of step 1 of the "ComputeNW" algorithm
disp('  Cut (stage 1): determine node weight in the constructed forest...')
CS1_st = tic;
% determine the initial weight of each supernode (1st step of Alg. S1)
M = length(EW); N = length(c);
NW=zeros(1,M); % NW(i) will store the number of nodes in cluster i.
for i = 1:N
    NW(c(i)) = NW(c(i)) + 1;
end

% Compute inNeighbors (i.e., step 2 of Alg. S1)
inNei = cell(1,M); % note: each cell is not initialized by the array size; this is a problem for speed as the elements are gradually added; but this is not a problem for c++ using advanced data structure; For Matlab, one can also consider to use a vector to count the number of elements of each cell first, based on which each cell can be initialized; 
for i = 1:M
    if pr1(i) ~= i
        inNei{1,pr1(i)}=[inNei{1,pr1(i)},i];
    end
end
 
% linearlize the nodes in the in-tree-based forest such that all the in-neighbors of each
% node have larger indexes than this node (this is fulfilled based on the
% Breadth First Search algorithm. (Step 3 of Alg. S1)
i = 1;
if size(pr1,1)>size(pr1,2)
    pr1 = pr1';
end
roots = find(pr1 == (1:M));
 
Ind = zeros(M,1);
ct = 0;
for j = 1:length(roots)       
    ct = ct + 1;
    Ind(ct) = roots(j); 
    while i<=ct
        temp = inNei{1,Ind(i)};
        if ~isempty(temp)     
            for m = 1:length(temp)                
                ct = ct + 1;
                Ind(ct) = temp(m);
            end
        end
        i = i + 1;
    end
end
                
% update the weight of each node on the in-tree-based forest (Step 4 of Alg. S1)
for j = M:-1:1  
    i = Ind(j);
    if pr1(i)~=i
        NW(pr1(i))=NW(pr1(i))+NW(i);
    end
end
 
%% cut the edges according to cluster number and MinSize
[~,ind]=sort(EW,'descend'); % Step 2 of Alg. S2; 
% the following is the step 3 of Alg. S2
cut_edge_num = 0; 
t = 1;
num_roots_initial = length(roots);
num_of_edges_in_graph = M - num_roots_initial;  
if nC >= num_roots_initial
    num_of_edges_required_to_remove = nC - num_roots_initial; 
else
    num_of_edges_required_to_remove = 0;
    warning('there could exist over-partitioning problem; it is suggest to increase the value of parameter k or increase cluster number');
end
passed_node = zeros(M,1);
CS1 = toc(CS1_st); 
disp(['  cost time on 1st stage of Cut: ',num2str(CS1),'s'])

disp('  Cut (stage 2): check the edges one by one in decreasing order of edge weight...')
CS2_st = tic;
while cut_edge_num ~= num_of_edges_required_to_remove && t <= num_of_edges_in_graph
    start_node = ind(t);
    end_node = pr1(start_node);
    if NW(start_node) > MinSize
        % search the root node of end node
        ct = 1; 
        % lines 10 to 15
        passed_node(ct) = end_node;  
        while end_node ~= pr1(end_node)
            end_node = pr1(end_node);  
            ct = ct + 1;
            passed_node(ct) = end_node;
        end
        
        % lines 16 to 24
        root_node_reached = end_node;
        if NW(root_node_reached)-NW(start_node) > MinSize
            NW(passed_node(1:ct)) = NW(passed_node(1:ct)) - NW(start_node);
            pr1(start_node) = start_node;
            cut_edge_num = cut_edge_num + 1;  
        end
    end
    t = t + 1;
end
pr2 = pr1;
rs2 = find(pr2 == (1:M)); 
CS2 = toc(CS2_st);
disp(['  time cost on the 2nd stage of Cut: ',num2str(CS2),'s']) 
end
