
%% the process of initalizing the population "pop"
% 第一维：1：c，代表c个社区
% 第二维：1：n，代表n个节点
% 第三维：1：NP，代表NP个个体，每个个体代表一个模糊社区划分，包含n个节点对c个社区的隶属度

function pop = inital_pop(n, c, NP)
    
 pop = rand(c,n,NP);  
% 对每个个体中的每一列分量进行归一化
% 限制每个节点对c个社区的隶属度值总和为1
for r=1:NP
    for i=1:n
        pop(:,i,r)=pop(:,i,r)./sum(pop(:,i,r));
    end
end
    
end

