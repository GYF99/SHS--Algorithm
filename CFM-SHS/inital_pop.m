
%% the process of initalizing the population "pop"
% ��һά��1��c������c������
% �ڶ�ά��1��n������n���ڵ�
% ����ά��1��NP������NP�����壬ÿ���������һ��ģ���������֣�����n���ڵ��c��������������

function pop = inital_pop(n, c, NP)
    
 pop = rand(c,n,NP);  
% ��ÿ�������е�ÿһ�з������й�һ��
% ����ÿ���ڵ��c��������������ֵ�ܺ�Ϊ1
for r=1:NP
    for i=1:n
        pop(:,i,r)=pop(:,i,r)./sum(pop(:,i,r));
    end
end
    
end

