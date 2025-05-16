function y = bound_brain_SOSFCD(x,c,n,groundtruth)

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
% 节点社区固定，对应隶属度应始终最大
for i=1:n
     [max_value,max_index]=max(x(:,i));        
     if  max_index ~= groundtruth(i)
         original_value = x(groundtruth(i),i);  % 节点i真实社区中的u值
         x(groundtruth(i),i) = max_value;  % 把最大u给节点i真实社区
         x(max_index,i) = original_value;
     end                 
end

y=x;
