function y = bound_SOSFCD(x,c,n)

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
      
y=x;

    
    