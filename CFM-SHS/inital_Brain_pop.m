
%% the process of initalizing the population "pop"
% ��һά��1��c������c������
% �ڶ�ά��1��n������n���ڵ�
% ����ά��1��NP������NP�����壬ÿ���������һ��ģ���������֣�����n���ڵ��c��������������

function pop = inital_Brain_pop(n, c, NP,groundtruth)
    
pop = rand(c,n,NP);  


 
% ÿ������c*n�У�c��������n���ڵ�
% ÿһ��=ÿ���ڵ��c�������������ȷֲ�
% ����ÿ���ڵ��c��������������ֵ�ܺ�Ϊ1
for r=1:NP
    for i=1:n
        % �ڵ�i������������u���=1
        pop(:,i,r)=pop(:,i,r)./sum(pop(:,i,r));
        % �ڵ������̶�����Ӧ������Ӧʼ�����
        [max_value,max_index]=max(pop(:,i,r));
        if  max_index ~= groundtruth(i)
            original_value = pop(groundtruth(i),i,r);  % �ڵ�i��ʵ�����е�uֵ
            pop(groundtruth(i),i,r) = max_value;  % �����u���ڵ�i��ʵ����
            pop(max_index,i,r) = original_value;
        end                 
    end
end
    


