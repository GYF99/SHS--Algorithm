function f = Qg(U,adj,n,c)

W=adj; % ��Ȩ���磬Ȩ�ؾ���W = �ڽӾ���adj
m=sum(W,2); % n���ڵ�Ķ�
m2 = sum(sum(W)); % m2=2M=156


Q = 0;

% ����Qg = ÿ�������ڵļ������ۼ�
for k=1:c
    % ����i,j����ÿ���ڵ�
    for i=1:n
        for j=1:n
            Q = Q + (W(i,j) - (m(i).*m(j))./m2).*U(k,i).*U(k,j);
        end
    end
    
end
Q = Q / m2;
    
f=Q;




    