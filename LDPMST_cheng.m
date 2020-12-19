% -------------------------------------------------------------------------
%Aim:
%The matlab code of "Clustering with local density peaks-based minimum spanning tree"
% -------------------------------------------------------------------------
%Input:
%A: the data set
%clu_num:the number of clusters
% -------------------------------------------------------------------------
%Output:
%cl: the clustering result
%time: the running time of LDP-MST
% -------------------------------------------------------------------------
% Written by Dongdong Cheng
% Department of Computer Science, Chongqing University 
% December 2017

function [cl,time] = LDPMST_cheng(A,clu_num,minsize)
%优化求交集和并集的代码
%A 数据集
%stopK聚类数
tic;
[N,dim]=size(A);
disp('Search local density peaks');
%首先搜索核心点
% [index,supk,max_nb,rho,local_core,cores,cl,cluster_number ] = CoreSearch4(A);
[index,supk,nb,rho,local_core,cores2,cl,cluster_number ] = LDP_Searching(A);
% [dist,index,supk,max_nb,rho,local_core,cores,cl,cluster_number] = CoreSearch2(A);
% [cores,local_core,cl,supk,cluster_number,index] = KeyPoint(A);
%利用核心点之间的距离构造最小生成树
%（1）得到最小生成树
%第一种方法是利用每个簇与其他簇的k近邻的交集
% [rho_sorted,ordrho]=sort(rho,'ascend');
% alpha=0.01;
% rho_threshold=0;%rho_sorted(floor(N*alpha));
% for i=1:cluster_number
%     if cores(i)~=0
%     if rho(cores(i))<rho_threshold %排除密度较小的核心点
%         mind=inf;
%         p=0;
%         for j=1:cluster_number
%             if i~=j
%                 x=A(cores(i),:);
%                 y=A(cores(j),:);
%                 distance=sqrt(sum((x-y).^2));
% %                 distance=pdist2(A(cores(i),:),A(cores(j),:));
%             if mind>distance&&rho(cores(j))>rho_threshold
%                 mind=distance;
%                 p=j;
%             end
%             end
%         end
%         for j=1:N
%             if local_core(j)==cores(i)
%                 local_core(j)=cores(p);
%             end
%         end
%     end
%     end
% end
% cluster_number=0;
% cl=zeros(N,1);
% for i=1:N
%     if local_core(i)==i;
%        cluster_number=cluster_number+1;
%        cores2(cluster_number)=i;
%        cl(i)=cluster_number;
%     end
% end
disp('初始子簇个数为：');disp(cluster_number);

% 以下是得出准核心直接得到的子簇
% for i=1:N
%     cl(i)=cl(local_core(i));
% end
%画出核心点图以及相应的初始聚类结果
% figure;
% plot(A(:,1),A(:,2),'ko','MarkerSize',5,'MarkerFaceColor','k');
% hold on;
% for i=1:N
%     plot([A(i,1),A(local_core(i),1)],[A(i,2),A(local_core(i),2)],'color','g','LineStyle',':','linewidth',1.5);
%     hold on;
% end
% plot(A(local_core,1),A(local_core,2),'ro','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','r');
% drawcluster2(A,cl,cluster_number);
% hold on;
% plot(A(local_core,1),A(local_core,2),'ro','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','r');
%hold off;
disp('Compute shared neighbors distance');
cdata=cell(1,cluster_number);%保存每个簇中都有哪些点
cdataexp=cell(1,cluster_number);%保存每个簇中的点及每个点中的k近邻
nc=zeros(1,cluster_number);%保存属于某个核心点的点数
ncexp=zeros(1,cluster_number);
% core_dist=zeros(cluster_number,cluster_number);
% for i=1:cluster_number
%     for j=i+1:cluster_number
%         x=A(cores2(i),:);
%         y=A(cores2(j),:);
%         d=sqrt(sum((x-y).^2));
%         core_dist(i,j)=d;
%         core_dist(j,i)=d;
%     end
% end
core_dist=pdist2(A(cores2,:),A(cores2,:));
core_edist=core_dist;
maxd=max(max(core_dist));
sd=zeros(cluster_number,1);
for i=1:cluster_number
cdata{1,i}=find(local_core==cores2(i))';
[~,csize]=size(cdata{1,i});
nc(i)=csize;
temp=index(cdata{1,i},2:supk+1);
if csize~=1
% temp2=union(temp,temp)';
cdataexp{1,i}=union(cdata{1,i},temp)';
else
temp2=union(temp,temp);
cdataexp{1,i}=union(cdata{1,i},temp2);
end
% cdataexp{1,i}=union(cdata{1,i},temp2);
end
 
for i=1:cluster_number
    for j=i+1:cluster_number
        inset1=intersect(cdataexp{1,i},cdataexp{1,j});
%         inset2=intersect(cdata{1,i},cdataexp{1,j});
%         hold on;
%         plot(A(inset1,1),A(inset1,2),'o','MarkerSize',6,'MarkerFaceColor','y','MarkerEdgeColor','y');
        averho=sum(rho(inset1));
        [~,numinset1]=size(inset1);
%         [~,numinset2]=size(inset2);
%         fprintf('第%d个点和第%d个点的交集数为：%d,密度和为%f\n',i,j,numinset1,averho);
        if numinset1==0%&&numinset2==0
            core_dist(i,j)=maxd*(core_dist(i,j)+1);
            core_dist(j,i)=core_dist(i,j);
        else
            core_dist(i,j)=core_dist(i,j)/(averho*numinset1);
            core_dist(j,i)=core_dist(i,j);
        end
        
    end
end



disp('Construct MST on local density peaks');
[cores_cl]=LMSTCLU_OPT( A(cores2,:),core_dist,nc,minsize,clu_num);
cl=zeros(N,1);
for i=1:cluster_number
    cl(cores2(i))=cores_cl(i);
end
for i=1:N
    cl(i)=cl(local_core(i));
end
disp('Complete!');
time=toc
% ncluster=max(cl);
% figure;drawcluster2(A,cl,ncluster);
end


