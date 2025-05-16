function f = Qg(U,adj,n,c)

W=adj; % 无权网络，权重矩阵W = 邻接矩阵adj
m=sum(W,2); % n个节点的度
m2 = sum(sum(W)); % m2=2M=156


Q = 0;

% 计算Qg = 每个社区内的计算结果累加
for k=1:c
    % 变量i,j遍历每个节点
    for i=1:n
        for j=1:n
            Q = Q + (W(i,j) - (m(i).*m(j))./m2).*U(k,i).*U(k,j);
        end
    end
    
end
Q = Q / m2;
    
f=Q;




    