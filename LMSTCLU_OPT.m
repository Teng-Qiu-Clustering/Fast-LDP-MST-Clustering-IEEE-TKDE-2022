% -------------------------------------------------------------------------
%Aim:
%Construct MST on local density peaks and cluster local densitypeaks
% -------------------------------------------------------------------------
%Input:
%A: local density peaks
%dist: shared neighbors-based distance matrix
%minsize:Minsize
%clu_num: the number of clusters
% -------------------------------------------------------------------------
%Output:
%cl: the clustering result of local density peaks
% -------------------------------------------------------------------------
% Written by Dongdong Cheng
% Department of Computer Science, Chongqing University 
% December 2017

function  [cl]=LMSTCLU_OPT( A,dist,nc,minsize,clu_num)
%A:原始数据集
%dist:构造最小生成树的距离矩阵
%minsize:最小簇的点数阈值
%clu_num:期望的聚类数
[N,dim]=size(A);
%首先构造最小生成树
%（2）得到最小生成树
SD=sparse(dist);
% UG=tril(SD);
% [ST,pred] = graphminspantree(UG,'METHOD','Prim');
ini_cluster_number = length(nc);
G = graph(SD);
[T,pred] = minspantree(G); 
ST = sparse(T.Edges.EndNodes(:,1),T.Edges.EndNodes(:,2),T.Edges.Weight,ini_cluster_number,ini_cluster_number);


% %画出最小生成树
% % figure(3);
% % plot(A(:,1),A(:,2),'.');
% % hold on;
% for i=1:N
%     j=pred(i);
%     if j~=0&&~isnan(j)
%         x=[A(i,1),A(j,1)];
%         y=[A(i,2),A(j,2)];
%         plot(x,y,'linewidth',3,'color','g');
%         hold on;
%     end
% end
% %hold off;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%确定最优聚类数
[xlabel,ylabel]=find(ST);
edges=zeros(N-1,1);
for i=1:N-1
    edges(i)=full(ST(xlabel(i),ylabel(i)));
end
[~,sortedind]=sort(edges,'descend');
% maxavecsp=0;
maxcv=-inf;

ST2=full(ST);
% pred2=pred;
% %画出最小生成树
%  figure;
%  plot(A(:,1),A(:,2),'.');
%  hold on;
% for i=1:N
%     j=pred(i);
%     if j~=0
%         x=[A(i,1),A(j,1)];
%         y=[A(i,2),A(j,2)];
%         plot(x,y);
%         hold on;
%     end
% end
i=0;
k=1;
while k<clu_num    
    %确定要去除的边
    s1=0;s2=0;
    p=xlabel(sortedind(i+1));
    q=ylabel(sortedind(i+1));
    while s1<minsize || s2<minsize
        s1=0;s2=0;
        visited=zeros(N,1);
        i=i+1;
        t=sortedind(i);
        xl(i)=xlabel(t);
        yl(i)=ylabel(t);
        %根据ST2矩阵和这两个点搜索连通子树
        p=xl(i);q=yl(i);
%         plot([A(p,1),A(q,1)],[A(p,2),A(q,2)],'r','linewidth',3);
%         hold on;
%         input('等一下');
        queue=zeros(N,1);
        front=1;%队列的头
        rear=1;%队列的尾
        queue(rear)=p;
        rear=rear+1;
%        fprintf('第一个子簇');
        while front~=rear
            temp=queue(front);
            s1=s1+nc(temp);
%             if s1>=minsize
%                 break;
%             end
            visited(temp)=1;
            front=front+1;
            for j=1:N
                if (ST2(temp,j)~=0||ST2(j,temp)~=0)&&j~=q&&visited(j)==0&&iscontain2(queue,front,rear,j)==0
                    queue(rear)=j;
                    rear=rear+1;
                end
            end
        end
%         fprintf('第二个子簇');
        queue=zeros(N,1);
        front=1;%队列的头
        rear=1;%队列的尾
        queue(rear)=q;
        rear=rear+1;
        while front~=rear
            temp=queue(front);
            s2=s2+nc(temp);
%             if s2>=minsize
%                 break;
%             end
            visited(temp)=1;
            front=front+1;
            for j=1:N
                if (ST2(temp,j)~=0||ST2(j,temp)~=0)&&j~=p&&visited(j)==0&&iscontain2(queue,front,rear,j)==0
                    queue(rear)=j;
                    rear=rear+1;
                end
            end
        end  
%         fprintf('s1的值为%d，s2的值为%d',s1,s2);
%         input('等一下');
    end
%     fprintf('\n');
    ST2(p,q)=0;
    ST2(q,p)=0;
    k=k+1; 
%     plot([A(p,1),A(q,1)],[A(p,2),A(q,2)],'y','linewidth',3);
%     hold on;
%     input('等一下');
end
%广度优先搜索得到最后结果
cl=zeros(N,1);
   ncl=0;
   sumedge=zeros(k,1);
   for i=1:N
       if cl(i)==0
           ncl=ncl+1;
           queue=zeros(N,1);
           front=1;%队列的头
           rear=1;%队列的尾
           queue(rear)=i;
           rear=rear+1;
           sumedge(ncl)=0;
           while front~=rear
               p=queue(front);
%                fprintf('p值有：%d',p);
               front=front+1;
               cl(p)=ncl;
               for j=1:N
                   if (ST2(p,j)~=0||ST2(j,p)~=0)&&cl(j)==0&&iscontain(queue,j)==0
                       queue(rear)=j;
                       rear=rear+1;
                   end
               end
           end
%            fprintf('\n');
       end
   end
end




