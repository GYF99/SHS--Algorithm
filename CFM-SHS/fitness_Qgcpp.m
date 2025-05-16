%% the procedure to compute the modularity Qg values of individuals in pop

function fit = fitness_Qgcpp(pop,NP,adj,n,c)

fit = zeros(NP,1);

% 计算每个个体即模糊重叠社区划分对应的Qg值
for r=1:NP
    W=adj; % 无权网络，权重矩阵W = 邻接矩阵adj
    m=sum(W,2); % n个节点的度
    m2 = sum(sum(W)); % m2=2M=156   
    fit(r)=Qgcpp(pop(:,:,r),adj,m,m2,n,c);
%     fit(r,2)=Qg(pop(:,:,r),adj,n,c);
    
end 
