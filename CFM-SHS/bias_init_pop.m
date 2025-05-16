
%% the biased initialization process of the population

function bias_pop = bias_init_pop( pop, c, n, NP, adj)

bias_pop = zeros(c,n,NP);

for r = 1:NP      % ������Ⱥ��ÿ������
    
    i_rand = floor(rand*n)+1;  % ���ѡ��һ���ڵ�
    
    % ���ڵ�i_rand�����������ȣ���ֵ�������ھӽڵ�
    
    for i = 1:n         % ���ڸ�����ÿ���ڵ�i
        
        if i~=i_rand         % ��ǰ�ڵ�i��ͬ��Ŀ��ڵ�i_rand
            
            if adj(i_rand,i) == 1  % ����ǰ�ڵ�i��Ŀ��ڵ�i_randΪ�ھӹ�ϵ
               
               pop(:,i,r) = pop(:,i_rand,r);  
               
%                % ȷ���ڵ�i����ʵ������������������
%                [max_value,max_index]=max(pop(:,i,r));
%               
%               if  max_index ~= groundtruth(i)
%                   original_value = pop(groundtruth(i),i,r);  % �ڵ�i��ʵ�����е�uֵ
%                   pop(groundtruth(i),i,r) = max_value;  % �����u���ڵ�i��ʵ����
%                   pop(max_index,i,r) = original_value;
%               end                                                                 
                  
            end
            
        end
        
    end 
    bias_pop(:,:,r) = pop(:,:,r);
    
end
