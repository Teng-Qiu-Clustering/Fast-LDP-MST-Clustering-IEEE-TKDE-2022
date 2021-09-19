function cl = Cut_MST_QT_v2(I_base,EW,NW,minsize,cluterN)
 
% Note: the input NW stores the weights of the root nodes in the initially
% constructed in-tree-based forest (i.e., the result of step 1 of the "ComputeNW" algorithm)

% Compute inNeighbors (i.e., step 2 of the "ComputeNW" algorithm)
M = length(NW);
inNei = cell(1,M);
for i = 1:M
    if I_base(i) ~= i
        inNei{1,I_base(i)}=[inNei{1,I_base(i)},i];
    end
end

% linearlize the nodes in the in-tree-based forest such that all the in-neighbors of each
% node have larger indexes than this node (this is fulfilled based on the
% "breadth first search" technology. 
i = 1;
roots = find(I_base == (1:M));
Ind = [];
for j = 1:length(roots)
    Ind = [Ind,roots(j)];
    while i<=length(Ind)
        if length(inNei{1,Ind(i)}) ~= 0
            Ind = [Ind,inNei{1,Ind(i)}];
        end
        i = i + 1;
    end
end
                
% update the weight of each node on the in-tree-based forest
for j = M:-1:1  
    i = Ind(j);
    if I_base(i)~=i
        NW(I_base(i))=NW(I_base(i))+NW(i);
    end
end
 
% disp('NW = ')
% disp(NW);
%% cut the edges according to cluster number and minsize

[W_sort,ind]=sort(EW,'descend');
cut_edge_num = 0;
start_node_of_removed_edge = []; % with regard to the in-tree
t = 1;
roots_last = roots;
num_roots_initial = length(roots);
num_of_edges_in_graph = M - num_roots_initial;  
if cluterN > num_roots_initial
    num_of_edges_required_to_remove = cluterN - num_roots_initial; 
else
    num_of_edges_required_to_remove = 0;
     warning('there could exist over-partitioning problem; it is suggest to increase the value of parameter k or increase cluster number');
end
while cut_edge_num ~= num_of_edges_required_to_remove && t <= num_of_edges_in_graph
    start_node = ind(t);
    end_node = I_base(start_node);
    if NW(start_node) > minsize
        % search the root node of end node
        passed_node = end_node;
        flag = 1;
        while flag
            parent_end_node = I_base(end_node);
            if end_node == parent_end_node
                flag = 0;
                root_node_reached = parent_end_node;
            else
                end_node = parent_end_node;
            end
            passed_node = [passed_node,parent_end_node];
        end
        
        if NW(root_node_reached)-NW(start_node)>minsize
            NW(passed_node) = NW(passed_node) - NW(start_node);
            I_base(start_node) = start_node;
            cut_edge_num = cut_edge_num + 1;
            start_node_of_removed_edge = [start_node_of_removed_edge;start_node];
            roots_last = [roots_last,start_node];
        end
    end
    t = t + 1;
end
 
% disp(start_node_of_removed_edge)

I_root = I_base(I_base);
while norm(I_root-I_base) ~=0
    I_base = I_root;
    I_root = I_base(I_base);
end
cl=zeros(1,M);
cl(roots_last) = 1:length(roots_last);
cl=cl(I_root);
end
