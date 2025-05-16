
%% the biased initialization process of the population

function bias_pop = bias_init_brain_pop_consv4( pop, c, n, NP, adj, groundtruth)

bias_pop = zeros(c,n,NP);
for r=1:NP
   for i=1:n
        % 节点i对所有社区的u求和=1
        pop(:,i,r)=pop(:,i,r)./sum(pop(:,i,r));
   end
end
for r = 1:NP      % 对于种群中每个个体
    
    i_rand = floor(rand*n)+1;  % 随机选择一个节点
    
    % 将节点i_rand的社区隶属度，赋值给所有邻居节点
    
    for i = 1:n         % 对于个体中每个节点i
        
        if i~=i_rand         % 当前节点i不同于目标节点i_rand
            
            if adj(i_rand,i) == 1  % 若当前节点i和目标节点i_rand为邻居关系
               
               pop(:,i,r) = pop(:,i_rand,r);  
               adj_nodes=[];
               adj_value=[];
               % 确定节点i的真实社区，有最大的隶属度
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
        
    end 
    bias_pop(:,:,r) = pop(:,:,r);
    
end
