
%Input:
%A: the data set
%clu_num:the number of clusters
% -------------------------------------------------------------------------
%Output:
%cl: the clustering result
%time: the running time of LDP-MST

% Written by Teng Qiu

function [cl,time] = LDPMST_QT_exact_fast_less_memory_v2(A,clu_num,minsize,knnMethod,initial_max_k)

tic;
[N,dim]=size(A);

[knnIndex,supk,nb,rho,local_core,cores,cl,ini_cluster_number] = LDP_Searching_by_QT(A,knnMethod,initial_max_k);
 
A_cores = A(cores,:);
M = ini_cluster_number;

disp('determine the size of each component of the initial forest...'); 
nc=zeros(1,M);
for i = 1:N
    nc(cl(i)) = nc(cl(i)) + 1;
end
% clear A

 
disp('determine multiple label matrix for each node');
% method 1: very slow, due to the reason that we cannot specify the number
% of non-zeros elements in the sparse matrix before hand. On the other
% hand, if we define ML as a full matrix, this will cost a lot of space,
% since M is usually very large when N is large. One can use density
% smoothing to solve this problem. 

% % ML is a multiple label matrix. ML(i,j) denote
% % whether ML is assigned with label j (if ML(i,j) = 1) or not (if ML(i,j) =
% % 0).Thus, a node i could have multiple labels. 

% ML = sparse(N,M); 
% for i = 1:N
%     ML(i,cl(i)) = 1;
% end
% for i = 1:N
%     
%       if rem(i,10000) == 0 || i == N
%         disp([num2str(i),'/',num2str(N)]) 
%       end
%     
%     for m = 2:supk+1
%         j = knnIndex(i,m);
%         if cl(i) ~= cl(j)
%             ML(j,cl(i)) = 1; 
%         end 
%     end
% end

% method 2: 
% K = supk + 1;
% I = zeros(N*K,1);
% J = zeros(N*K,1);
% V = ones(N*K,1);
% t = 1;
% for i = 1:N
%     
%       if rem(i,10000) == 0 || i == N
%         disp([num2str(i),'/',num2str(N)]) 
%       end     
% 
%     for m = 1:K
%         j = knnIndex(i,m);
%         J(t) = cl(i);
%         I(t) = j;    
%         t = t + 1;
%     end
% end
% ML = sparse(I,J,V,N,M);

% method 3: vectorization of method 2
K = supk + 1;
% I = knnIndex(:);
% J = cl(repmat((1:N)',K,1));
% V = ones(N*K,1); 
% ML = sparse(I,J,V,N,M);
 ML = sparse(knnIndex(:),cl(repmat((1:N)',K,1)),ones(N*K,1),N,M); 
 
disp(['compute distance matrix of size:',num2str(ini_cluster_number),' x ',num2str(ini_cluster_number),'...']);
% Method 1: very slow for large dataset
% mu_distance_matrix=zeros(M,M);
% for s=1:M-1
%     for t=s+1:M
%         d=sqrt(sum((A_cores(s,:)-A_cores(t,:)).^2));
%         mu_distance_matrix(s,t) = d/sum(ML(:,s).*ML(:,t).*rho')/(ML(:,s)'*ML(:,t));
%     end
% end  
%  G = graph(mu_distance_matrix,'upper'); 

% method 2:
disp('determine multiple labels of each node...')
[I,J] = find(ML); % J: cluster label vector, where J(t) is a label of node I(t)
cell_L = cell(1,N); %cell_L{1,i} is the set of neighboring set showing all the labels of node i
for t = 1:length(I)
    cell_L{1,I(t)} = [cell_L{1,I(t)},J(t)]; 
end

disp('compute sum rho..');
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
 
disp(['compute shared sparse distance matrix of size:',num2str(ini_cluster_number),' x ',num2str(ini_cluster_number),'...']);
   
[I,J,V] = find(rho_matrix);
[I2,J2,V2] = find(count_matrix);
for t = 1:length(I)
    d=sqrt(sum((A_cores(I(t),:)-A_cores(J(t),:)).^2)); 
    V(t) = d/V(t)/V2(t);
end
 G = graph(I,J,V,ini_cluster_number); 
  
disp('Construct MST on cores only ...');
[T,pred] = minspantree(G,'Type','forest','Method','sparse');

ST = sparse(T.Edges.EndNodes(:,1),T.Edges.EndNodes(:,2),T.Edges.Weight,ini_cluster_number,ini_cluster_number);
if any(isnan(pred))
    % % ****** Note:  there is problem for the output 'pred' when UG is unconnected for matlab 2013; this problem has been modified by matlab 2016 and above
    error('a higher version of matlab is needed (e.g., matlab 2016 and later version than that)');
end 

disp('transform MST to in-tree...')
idx_roots = find(pred == 0);
pred(idx_roots) = idx_roots;
I_base = pred;
W = zeros(1,ini_cluster_number); 
for i = 1:ini_cluster_number
W(i) = max(ST(i,I_base(i)),ST(I_base(i),i));
end

disp('Cut the tree...')

cores_cl = Cut_MST_QT_v2(I_base,W,nc,minsize,clu_num);

disp('final cluster assignment...');
cl=zeros(N,1);
cl(cores)=cores_cl;
cl=cl(local_core);
 
disp('End'); 
time=toc 
end


