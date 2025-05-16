
%% the process of initalizing the population "pop"
% 第一维：1：c，代表c个社区
% 第二维：1：n，代表n个节点
% 第三维：1：NP，代表NP个个体，每个个体代表一个模糊社区划分，包含n个节点对c个社区的隶属度

function pop = inital_Brain_pop_1(n, c, NP,groundtruth)
    
 pop = rand(c,n,NP);  
 %%%%%%替换%%%%%%%%%%%%%%%%%%%%%%%%%
 for k =1:NP
   % p=repmat(minVar,c,1)+rand(c,n).*repmat((maxVar-minVar),c,1);
       % 加载数据
%     data = load('cora_comm0.4645.mat');
    data = load('LFR5000_d50_dmax250_1.mat');
    % 获取字段的名称
    field_name = fieldnames(data);
    % 获取数据矩阵
    matrix = data.(field_name{1});
%     % 转置数据矩阵
    transposed_matrix = matrix';
    % 将转置后的矩阵赋值给变量 p
    p = transposed_matrix;
    pop(:,:,NP) = p;

end

 
% 每个个体c*n列，c个社区，n个节点
% 每一列=每个节点对c个社区的隶属度分布
% 限制每个节点对c个社区的隶属度值总和为1
for r=1:NP
    for i=1:n
        % 节点i对所有社区的u求和=1
        pop(:,i,r)=pop(:,i,r)./sum(pop(:,i,r));
        % 节点社区固定，对应隶属度应始终最大
        [max_value,max_index]=max(pop(:,i,r));
        if  max_index ~= groundtruth(i)
            original_value = pop(groundtruth(i),i,r);  % 节点i真实社区中的u值
            pop(groundtruth(i),i,r) = max_value;  % 把最大u给节点i真实社区
            pop(max_index,i,r) = original_value;
        end                 
    end
end
    


