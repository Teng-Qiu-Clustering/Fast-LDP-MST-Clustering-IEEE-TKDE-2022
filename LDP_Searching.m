
% -------------------------------------------------------------------------
%Aim:
%Search local density peaks
% -------------------------------------------------------------------------
%Input:
%A: the data set
% -------------------------------------------------------------------------
% Written by Dongdong Cheng
% Department of Computer Science, Chongqing University 
% December 2017

function [index,supk,nb,rho,local_core,cores,cl,cluster_number ] = LDP_Searching(A)
[N,dim]=size(A);
% dist=pdist2(A,A);
%  [~,index]=sort(dist,2);%对dist按行进行排序

kdtree=KDTreeSearcher(A(:,:),'bucketsize',1); % 1表示选用欧式距离
[index,~] = knnsearch(kdtree,A(:,:),'k',100);% 返回排好序的索引和距离
[~,ncolIndex] = size(index);
%初始化基本数据
r=1; 
flag=0;         
nb=zeros(1,N);  %自然邻居个数 
%NNN=zeros(N,N); %各点的自然邻居集
count=0;        %自然最近邻数为零的数据量连续相同的次数
count1=0;       %前一次自然最近邻数为零的数据量
count2=0;       %此次自然最近邻数为零的数据量

%搜索自然最近邻居
while flag==0
    for i=1:N
        k=index(i,r+1);
        nb(k)=nb(k)+1;
%         RNN(k,nb(k))=i;
    end
    r=r+1;
%     count2=0;
    [~,count2]=size(find(nb==0));
%     for i=1:N
%         if nb(i)==0
%             count2=count2+1;
%         end
%     end
    %计算nb(i)=0的点的数量连续不变化的次数
    if count1==count2
        count=count+1;
    else
        count=1;
    end
    if count2==0 || (r>2 && count>=2)   %邻居搜索终止条件
        flag=1;
    end
    

    if r == ncolIndex
        [index,~] = knnsearch(kdtree,A(:,:),'k',ncolIndex+50);%  there is a bug in the cheng's original code, modified by the code in this line:
    end

    count1=count2;
end

%计算自然最近邻的各种特征量
supk=r-1;               %最终K值，也是自然最近邻居的平均数
max_nb=max(nb);         %自然邻居的最大数目
min_nb=min(nb);         %自然邻居的最小数目
%NN=index(:,2:SUPk+1);   %各数据点的K近邻数据点集
%ratio_nb=nb./(N*SUPk);  %各数据点的自然最近邻居数目所占比例
%计算每个数据点的密度
%disp(SUPk);
%构造连接矩阵
disp('end of natural neighbor fiding')
disp(supk);
disp(max_nb);
% disp(min_nb);
rho=zeros(N,1);
Non=max_nb;

 
if Non>ncolIndex
    [index,~] = knnsearch(kdtree,A(:,:),'k',Non+1);% there is a bug in the cheng's original code, modified by the code in this line:
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %构造自然邻域图
% conn=zeros(N,N);
% for i=1:N
%     for j=2:supk+1
%         x=index(i,j);
%         dist=sqrt(sum((A(i,:)-A(x,:)).^2));
%         conn(i,x)=1/(1+dist);%距离的倒数作为两点的相似度
%         conn(x,i)=conn(i,x);
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N
    d=0;
    for j=1:Non+1
        x=index(i,j);
        dist=sqrt(sum((A(i,:)-A(x,:)).^2));
        d=d+dist;
    end
%     rho(i)=exp(-d^2/Non);
    rho(i)=nb(i)/d;
end
% [rho_sorted,ordrho]=sort(rho,'descend');%ordrho就是密度从大到小的顺序
local_core=zeros(N,1);%存放n个点的局部核心点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%找到每个点的代表点：搜索每个点的k邻域，将其中密度最大的点作为其代表点
% [rho_sorted,ordrho]=sort(rho,'descend');
% for i=1:N
%     rep=ordrho(i);
%     xrho=rho(rep);
%   for j=2:supk+1
%       if xrho<rho(index(ordrho(i),j))
%           xrho=rho(index(ordrho(i),j));
%           rep=index(ordrho(i),j);
%       end
%   end
%   local_core(ordrho(i))=rep;
% end
%  supk=5;
for i=1:N
    rep=i;
    xrho=rho(rep);
  for j=2:supk+1
      if xrho<rho(index(i,j))
          xrho=rho(index(i,j));
          rep=index(i,j);
      end
  end
  local_core(i)=rep;
end
% %画出每个点及其代表点
% figure;
% plot(A(:,1),A(:,2),'ko','MarkerSize',5,'MarkerFaceColor','k');
% hold on;
% for i=1:N
%     plot([A(i,1),A(local_core(i),1)],[A(i,2),A(local_core(i),2)],'linewidth',1.5,'color','k','LineStyle',':');
%     hold on;
% end
% plot(A(local_core,1),A(local_core,2),'rs','MarkerSize',6,'MarkerFaceColor','w');
%重新分配每个点的局部代表点
visited=zeros(N,1);
round=0;
for k=1:N
    if visited(k)==0
        parent=k;
        round=round+1;
        while local_core(parent)~=parent
            visited(parent)=round;
            parent=local_core(parent);
        end
        local_core(find(visited==round))=parent;
    end
end
% %画出每个点及其代表点
% hold on;
% plot(A(local_core,1),A(local_core,2),'ro','MarkerSize',8,'MarkerFaceColor','r');
% figure;
% plot(A(:,1),A(:,2),'ko','MarkerSize',5,'MarkerFaceColor','k');
% hold on;
% for i=1:N
%     plot([A(i,1),A(local_core(i),1)],[A(i,2),A(local_core(i),2)],'linewidth',1.5,'color','k','LineStyle',':');
%     hold on;
% end
% plot(A(local_core,1),A(local_core,2),'ro','MarkerSize',8,'MarkerFaceColor','r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%得到准核心点，准核心点以其自身为核心点
 cluster_number=0;
 cl=zeros(N,1);
for i=1:N
    if local_core(i)==i
       cluster_number=cluster_number+1;
       cores(cluster_number)=i;
       cl(i)=cluster_number;
    end
end
% disp('初始子簇个数为：');disp(cluster_number);
% 以下是得出准核心直接得到的子簇
for i=1:N
    cl(i)=cl(local_core(i));
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(6);
% % 画出核心点图以及相应的初始聚类结果
% plot(A(:,1),A(:,2),'.');
% hold on;
% for i=1:N
%     plot([A(i,1),A(local_core(i),1)],[A(i,2),A(local_core(i),2)]);
%     hold on;
% end
% drawcluster2(A,cl,cluster_number+1);
% hold on;
% plot(A(local_core,1),A(local_core,2),'ro','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','r');
% % title('搜索supk近邻的结果');
end



