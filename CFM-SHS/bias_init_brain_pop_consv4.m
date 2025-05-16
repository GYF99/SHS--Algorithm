
%% the biased initialization process of the population

function bias_pop = bias_init_brain_pop_consv4( pop, c, n, NP, adj, groundtruth)

bias_pop = zeros(c,n,NP);
for r=1:NP
   for i=1:n
        % �ڵ�i������������u���=1
        pop(:,i,r)=pop(:,i,r)./sum(pop(:,i,r));
   end
end
for r = 1:NP      % ������Ⱥ��ÿ������
    
    i_rand = floor(rand*n)+1;  % ���ѡ��һ���ڵ�
    
    % ���ڵ�i_rand�����������ȣ���ֵ�������ھӽڵ�
    
    for i = 1:n         % ���ڸ�����ÿ���ڵ�i
        
        if i~=i_rand         % ��ǰ�ڵ�i��ͬ��Ŀ��ڵ�i_rand
            
            if adj(i_rand,i) == 1  % ����ǰ�ڵ�i��Ŀ��ڵ�i_randΪ�ھӹ�ϵ
               
               pop(:,i,r) = pop(:,i_rand,r);  
               adj_nodes=[];
               adj_value=[];
               % ȷ���ڵ�i����ʵ������������������
               [max_value,max_index]=max(pop(:,i,r));           
                if  (max_index ~= groundtruth(i))
                % ȡ����i���ڵ�����������������ֵ
                node_c_list = pop(:,i,r);
                % ��ѯ�ڵ�i�������ھӽڵ�����ʵ�����е�ƽ��ֵ,ʹ��ƽ��ֵ����
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
                    %ɾ����ʵ������������ֵ
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
