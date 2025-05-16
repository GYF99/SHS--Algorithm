%consv4  不符合约束的节点隶属度调整为相同真实社区的邻居节点的平均值（整行赋值）
%最大社区修正后，节点的其他隶属度按照节点与社区的连接紧密程度重新排序
function y = bound_brain_SOSFCD_consv4(x,c,n,groundtruth,adj)

% 将每个元素约束到【0，1】
for k=1:c
    for i=1:n
        if x(k,i)>1
           x(k,i) = 0.9999;
        elseif x(k,i)<0
           x(k,i) = 0.0001;
        end    
    end
end
% 将隶属度矩阵的每一列归一化
% 每个节点对所有社区的隶属度总和 = 1
for i = 1 :n
    x(:,i) = x(:,i)./sum(x(:,i));
end

for i=1:n
    adj_nodes=[];
    adj_value=[];
    [max_value,max_index]=max(x(:,i));
    if  (max_index ~= groundtruth(i))
        % 取出第i个节点所有社区的隶属度值
        node_c_list = x(:,i);
        % 查询节点i的所有邻居节点在真实社区中的平均值,使用平均值代替
        adj_nodes = find(adj(i,:));
        len = 0;
        countlines=zeros(1,c);
        for v = adj_nodes
            [max_value_v,max_index_v]=max(x(:,v));
            if  (max_index_v == groundtruth(v)) && (groundtruth(i)==groundtruth(v))
                len = len + 1;
                adj_value(:,len)= x(:,v);                      
            end
            countlines(groundtruth(v))=countlines(groundtruth(v))+1;
        end
        if len ~= 0 
            average_value = mean(adj_value,2);  
            [average,indexa]=sort(average_value, 'descend');
            x(groundtruth(i),i)=average(1);
            %删除真实社区的隶属度值
            average(1)=[];
            indexa(1)=[];
            
            [counts,index]=sort(countlines, 'descend');   
            temp=find(index==groundtruth(i));
            index(temp)=[];
            for q=1:size(average,1)
                x(index(q),i)=average(q);
            end  
        else
            original_value = x(groundtruth(i),i); 
            x(groundtruth(i),i) = max_value;
            x(max_index,i) = original_value;
        end
   
    end
end  

% 节点社区固定，对应隶属度应始终最大

y=x;
