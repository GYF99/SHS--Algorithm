
%% the biased initialization process of the population

function bias_pop = bias_init_pop( pop, c, n, NP, adj)

bias_pop = zeros(c,n,NP);

for r = 1:NP      % 对于种群中每个个体
    
    i_rand = floor(rand*n)+1;  % 随机选择一个节点
    
    % 将节点i_rand的社区隶属度，赋值给所有邻居节点
    
    for i = 1:n         % 对于个体中每个节点i
        
        if i~=i_rand         % 当前节点i不同于目标节点i_rand
            
            if adj(i_rand,i) == 1  % 若当前节点i和目标节点i_rand为邻居关系
               
               pop(:,i,r) = pop(:,i_rand,r);  
               
%                % 确定节点i的真实社区，有最大的隶属度
%                [max_value,max_index]=max(pop(:,i,r));
%               
%               if  max_index ~= groundtruth(i)
%                   original_value = pop(groundtruth(i),i,r);  % 节点i真实社区中的u值
%                   pop(groundtruth(i),i,r) = max_value;  % 把最大u给节点i真实社区
%                   pop(max_index,i,r) = original_value;
%               end                                                                 
                  
            end
            
        end
        
    end 
    bias_pop(:,:,r) = pop(:,:,r);
    
end
