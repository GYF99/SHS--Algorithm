
%% the process of initalizing the population "pop"
% 第一维：1：c，代表c个社区
% 第二维：1：n，代表n个节点
% 第三维：1：NP，代表NP个个体，每个个体代表一个模糊社区划分，包含n个节点对c个社区的隶属度
%consv4  不符合约束的节点隶属度调整为相同真实社区的邻居节点的平均值（整行赋值）
%最大社区修正后，节点的其他隶属度按照节点与社区的连接紧密程度重新排序
function pop = inital_Brain_pop_consv4(n, c, NP,groundtruth,adj)
    
 pop = rand(c,n,NP);  
 %%%%替换%%%%%%%%%%%%%%%%%%%%%%
 
% for k =1:NP
%    % p=repmat(minVar,c,1)+rand(c,n).*repmat((maxVar-minVar),c,1);
%        % 加载数据
% %     data = load('cora_comm0.4645.mat');
%     data = load('LFR10000_d100_dmax500_1.mat');
%     % 获取字段的名称
%     field_name = fieldnames(data);
%     % 获取数据矩阵
%     matrix = data.(field_name{1});
% %     % 转置数据矩阵
%     transposed_matrix = matrix';
%     % 将转置后的矩阵赋值给变量 p
%     p = transposed_matrix;
%     pop(:,:,NP) = p;
% 
% end

 
% 每个个体c*n列，c个社区，n个节点
% 每一列=每个节点对c个社区的隶属度分布
% 限制每个节点对c个社区的隶属度值总和为1
for r=1:NP
   for i=1:n
        % 节点i对所有社区的u求和=1
        pop(:,i,r)=pop(:,i,r)./sum(pop(:,i,r));
   end
end
for r=1:NP
   for i=1:n
        adj_nodes=[];
        adj_value=[];
        % 节点社区固定，对应隶属度应始终最大
        [max_value,max_index]=max(pop(:,i,r));
        if  (max_index ~= groundtruth(i))
        % 取出第i个节点所有社区的隶属度值
        node_c_list = pop(:,i,r);
        % 查询节点i的所有邻居节点在真实社区中的平均值,使用平均值代替
        adj_nodes = find(adj(i,:));
        len = 0;
        countlines=zeros(1,c);
        
        for v = adj_nodes
            [max_value_v,max_index_v]=max(pop(:,v,r));
            if  (max_index_v == groundtruth(v)) && (groundtruth(i)==groundtruth(v))
                len = len + 1;
                adj_value(:,len)= pop(:,v,r);                      
            end
            countlines(groundtruth(v))=countlines(groundtruth(v))+1;
        end
        if len ~= 0 
            average_value = mean(adj_value,2);  
            [average,indexa]=sort(average_value, 'descend');
            pop(groundtruth(i),i,r)=average(1);
            %删除真实社区的隶属度值
            average(1)=[];
            indexa(1)=[];
            
            [counts,index]=sort(countlines, 'descend');   
            temp=find(index==groundtruth(i));
            index(temp)=[];
            for q=1:size(average,1)
                pop(index(q),i,r)=average(q);
            end  
        else
            original_value = pop(groundtruth(i),i,r); 
            pop(groundtruth(i),i,r) = max_value;
            pop(max_index,i,r) = original_value;
        end
   
        
        end                 
    end
end
    


