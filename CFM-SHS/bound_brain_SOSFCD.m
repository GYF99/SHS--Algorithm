function y = bound_brain_SOSFCD(x,c,n,groundtruth)

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
% �ڵ������̶�����Ӧ������Ӧʼ�����
for i=1:n
     [max_value,max_index]=max(x(:,i));        
     if  max_index ~= groundtruth(i)
         original_value = x(groundtruth(i),i);  % �ڵ�i��ʵ�����е�uֵ
         x(groundtruth(i),i) = max_value;  % �����u���ڵ�i��ʵ����
         x(max_index,i) = original_value;
     end                 
end

y=x;
