
%% the process of initalizing the population "pop"
% ��һά��1��c������c������
% �ڶ�ά��1��n������n���ڵ�
% ����ά��1��NP������NP�����壬ÿ���������һ��ģ���������֣�����n���ڵ��c��������������

function pop = inital_Brain_pop_1(n, c, NP,groundtruth)
    
 pop = rand(c,n,NP);  
 %%%%%%�滻%%%%%%%%%%%%%%%%%%%%%%%%%
 for k =1:NP
   % p=repmat(minVar,c,1)+rand(c,n).*repmat((maxVar-minVar),c,1);
       % ��������
%     data = load('cora_comm0.4645.mat');
    data = load('LFR5000_d50_dmax250_1.mat');
    % ��ȡ�ֶε�����
    field_name = fieldnames(data);
    % ��ȡ���ݾ���
    matrix = data.(field_name{1});
%     % ת�����ݾ���
    transposed_matrix = matrix';
    % ��ת�ú�ľ���ֵ������ p
    p = transposed_matrix;
    pop(:,:,NP) = p;

end

 
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
    


