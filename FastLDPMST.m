function [c,time] = FastLDPMST(data,clu_num,minsize,knnMethod,initial_max_k)
%Input:
%  data: data set (nrow: samples; col: features);
%  clu_num: number of clusters;
%Output:
%  c: cluster label vector;
%  time: running time (seconds);

% Written by Teng Qiu (UESTC)
 
tic;
%% initial clustering by LDP
[knnIndex,supk,rho,r,rs,c] = LDP_Searching_by_QT(data,knnMethod,initial_max_k);
 
A_cores = data(rs,:);
M = length(rs);

disp('determine the size of each component of the initial forest...'); 
nc=zeros(1,M);
[N,~]=size(data);
for i = 1:N
    nc(c(i)) = nc(c(i)) + 1;
end
% clear data 
 
%% determine multiple label matrix ML
disp('determine multiple label matrix for each node...'); 
% method 3: vectorization of method 2
K = supk + 1; 
ML = sparse(knnIndex(:),c(repmat((1:N)',K,1)),ones(N*K,1),N,M); % Equivalent to: I = knnIndex(:);J = c(repmat((1:N)',K,1));V = ones(N*K,1); ML = sparse(I,J,V,N,M);

%% compute shared sparse distance matrix 1
disp(['compute distance matrix of size:',num2str(M),' x ',num2str(M),'...']);
% method 2:
disp('determine multiple labels of each node...')
[I,J] = find(ML); % J: cluster label vector, where J(t) is a label of node I(t)
cell_L = cell(1,N); %cell_L{1,i} is the set of neighboring set showing all the labels of node i
for t = 1:length(I)
    cell_L{1,I(t)} = [cell_L{1,I(t)},J(t)]; 
end
 
disp('sum over density...');
rho_matrix=sparse(M,M);
count_matrix = sparse(M,M);
for j = 1:N
    cell_L_j = cell_L{1,j};
    cell_L_size = length(cell_L_j);
    if cell_L_size~= 1 
        for s=1:cell_L_size-1            
            for t=s+1:cell_L_size
                a = max(cell_L_j(s),cell_L_j(t)); b = min(cell_L_j(s),cell_L_j(t));
                rho_matrix(a,b) = rho_matrix(a,b) + rho(j);
                count_matrix(a,b) = count_matrix(a,b) + 1;
            end
        end
    end
end
 
%% compute shared sparse distance matrix 2
disp(['compute shared sparse distance matrix of size:',num2str(M),' x ',num2str(M),'...']);
[I,J,V] = find(rho_matrix);
[~,~,V2] = find(count_matrix);
for t = 1:length(I)
    d=sqrt(sum((A_cores(I(t),:)-A_cores(J(t),:)).^2)); 
    V(t) = d/V(t)/V2(t);
end
 G = graph(I,J,V,M); 
  
%% Construct MSF
disp('Construct MSF on root nodes (rs) only ...');
[T,pred] = minspantree(G,'Type','forest','Method','sparse');

ST = sparse(T.Edges.EndNodes(:,1),T.Edges.EndNodes(:,2),T.Edges.Weight,M,M);
if any(isnan(pred))
    % % ****** Note:  there is problem for the output 'pred' when UG is unconnected for matlab 2013; this problem has been modified by matlab 2016 and above
    error('a higher version of matlab is needed (e.g., matlab 2016 and later version than that)');
end 

%% transform MST to in-tree-based MSF
disp('transform MSF to in-tree-based MSF...')
idx_roots = find(pred == 0);
pred(idx_roots) = idx_roots;
I_base = pred;
W = zeros(1,M); 
for i = 1:M
W(i) = max(ST(i,I_base(i)),ST(I_base(i),i));
end

%% Cut the tree
disp('Cut the tree...') 
cores_cl = Cut_MST_QT_v2(I_base,W,nc,minsize,clu_num);

%% final cluster assignment 
disp('final cluster assignment...');
c=zeros(N,1);
c(rs)=cores_cl;
c=c(r);
 
disp('End'); 
time=toc; 
end


