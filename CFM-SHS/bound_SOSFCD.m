function y = bound_SOSFCD(x,c,n)

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
      
y=x;

    
    