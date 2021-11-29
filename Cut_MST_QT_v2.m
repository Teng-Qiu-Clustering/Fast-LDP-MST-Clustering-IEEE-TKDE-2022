function cl = Cut_MST_QT_v2(I_base,EW,NW,minsize,cluterN)
 
% Note: the input NW stores the weights of the root nodes in the initially
% constructed in-tree-based forest (i.e., the result of step 1 of the "ComputeNW" algorithm)
disp('Cut (stage1): determine node weight in the constructed forest...')
CS1_st = tic;
% Compute inNeighbors (i.e., step 2 of the "ComputeNW" algorithm)
M = length(NW);
inNei = cell(1,M); % note: each cell is not initilized by the array size; this is a problem for speed as the elements are gradually added; but this is not a problem for c++ using advanced data structure; For matlab, one can also consider to use a vector to count the number of elements of each cell first, based on which each cell can be initialized; 
for i = 1:M
    if I_base(i) ~= i
        inNei{1,I_base(i)}=[inNei{1,I_base(i)},i];
    end
end
 
% linearlize the nodes in the in-tree-based forest such that all the in-neighbors of each
% node have larger indexes than this node (this is fulfilled based on the
% "breadth first search" technology. 
i = 1;
if size(I_base,1)>size(I_base,2)
    I_base = I_base';
end
roots = find(I_base == (1:M));
 
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
                
% update the weight of each node on the in-tree-based forest
for j = M:-1:1  
    i = Ind(j);
    if I_base(i)~=i
        NW(I_base(i))=NW(I_base(i))+NW(i);
    end
end
 
%% cut the edges according to cluster number and minsize

[~,ind]=sort(EW,'descend');
cut_edge_num = 0;
% start_node_of_removed_edge = []; % with regard to the in-tree
t = 1;
roots_last = roots;
num_roots_initial = length(roots);
num_of_edges_in_graph = M - num_roots_initial;  
if cluterN >= num_roots_initial
    num_of_edges_required_to_remove = cluterN - num_roots_initial; 
else
    num_of_edges_required_to_remove = 0;
    warning('there could exist over-partitioning problem; it is suggest to increase the value of parameter k or increase cluster number');
end
passed_node = zeros(M,1);
CS1 = toc(CS1_st); disp(['cost time on 1st stage of Cut: ',num2str(CS1)])

disp('Cut (stage2): check the edges one by one in decreasing order of edge weight...')
CS2_st = tic;
while cut_edge_num ~= num_of_edges_required_to_remove && t <= num_of_edges_in_graph
    start_node = ind(t);
    end_node = I_base(start_node);
    if NW(start_node) > minsize
        % search the root node of end node
        ct = 1; 
        passed_node(ct) = end_node;  
        while end_node ~= I_base(end_node)
            end_node = I_base(end_node);  
            ct = ct + 1;
            passed_node(ct) = end_node;
        end
        root_node_reached = end_node;
        
        if NW(root_node_reached)-NW(start_node)>minsize
            NW(passed_node(1:ct)) = NW(passed_node(1:ct)) - NW(start_node);
            I_base(start_node) = start_node;
            cut_edge_num = cut_edge_num + 1;
%             start_node_of_removed_edge = [start_node_of_removed_edge;start_node];
            roots_last = [roots_last,start_node];
        end
    end
    t = t + 1;
end
CS2 = toc(CS2_st);
disp(['time cost on the 2nd stage of Cut: ',num2str(CS2)])
 
I_root = I_base(I_base);
while norm(I_root-I_base) ~=0
    I_base = I_root;
    I_root = I_base(I_base);
end
cl=zeros(1,M);
cl(roots_last) = 1:length(roots_last);
cl=cl(I_root);
end
