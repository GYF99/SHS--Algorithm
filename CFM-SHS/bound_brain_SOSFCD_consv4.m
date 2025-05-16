%consv4  ������Լ���Ľڵ������ȵ���Ϊ��ͬ��ʵ�������ھӽڵ��ƽ��ֵ�����и�ֵ��
%������������󣬽ڵ�����������Ȱ��սڵ������������ӽ��̶ܳ���������
function y = bound_brain_SOSFCD_consv4(x,c,n,groundtruth,adj)

% ��ÿ��Ԫ��Լ������0��1��
for k=1:c
    for i=1:n
        if x(k,i)>1
           x(k,i) = 0.9999;
        elseif x(k,i)<0
           x(k,i) = 0.0001;
        end    
    end
end
% �������Ⱦ����ÿһ�й�һ��
% ÿ���ڵ�������������������ܺ� = 1
for i = 1 :n
    x(:,i) = x(:,i)./sum(x(:,i));
end

for i=1:n
    adj_nodes=[];
    adj_value=[];
    [max_value,max_index]=max(x(:,i));
    if  (max_index ~= groundtruth(i))
        % ȡ����i���ڵ�����������������ֵ
        node_c_list = x(:,i);
        % ��ѯ�ڵ�i�������ھӽڵ�����ʵ�����е�ƽ��ֵ,ʹ��ƽ��ֵ����
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
            %ɾ����ʵ������������ֵ
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

% �ڵ������̶�����Ӧ������Ӧʼ�����

y=x;
