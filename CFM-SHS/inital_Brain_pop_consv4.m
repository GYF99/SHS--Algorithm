
%% the process of initalizing the population "pop"
% ��һά��1��c������c������
% �ڶ�ά��1��n������n���ڵ�
% ����ά��1��NP������NP�����壬ÿ���������һ��ģ���������֣�����n���ڵ��c��������������
%consv4  ������Լ���Ľڵ������ȵ���Ϊ��ͬ��ʵ�������ھӽڵ��ƽ��ֵ�����и�ֵ��
%������������󣬽ڵ�����������Ȱ��սڵ������������ӽ��̶ܳ���������
function pop = inital_Brain_pop_consv4(n, c, NP,groundtruth,adj)
    
 pop = rand(c,n,NP);  
 %%%%�滻%%%%%%%%%%%%%%%%%%%%%%
 
% for k =1:NP
%    % p=repmat(minVar,c,1)+rand(c,n).*repmat((maxVar-minVar),c,1);
%        % ��������
% %     data = load('cora_comm0.4645.mat');
%     data = load('LFR10000_d100_dmax500_1.mat');
%     % ��ȡ�ֶε�����
%     field_name = fieldnames(data);
%     % ��ȡ���ݾ���
%     matrix = data.(field_name{1});
% %     % ת�����ݾ���
%     transposed_matrix = matrix';
%     % ��ת�ú�ľ���ֵ������ p
%     p = transposed_matrix;
%     pop(:,:,NP) = p;
% 
% end

 
% ÿ������c*n�У�c��������n���ڵ�
% ÿһ��=ÿ���ڵ��c�������������ȷֲ�
% ����ÿ���ڵ��c��������������ֵ�ܺ�Ϊ1
for r=1:NP
   for i=1:n
        % �ڵ�i������������u���=1
        pop(:,i,r)=pop(:,i,r)./sum(pop(:,i,r));
   end
end
for r=1:NP
   for i=1:n
        adj_nodes=[];
        adj_value=[];
        % �ڵ������̶�����Ӧ������Ӧʼ�����
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
    


