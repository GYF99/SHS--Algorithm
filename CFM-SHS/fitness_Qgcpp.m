%% the procedure to compute the modularity Qg values of individuals in pop

function fit = fitness_Qgcpp(pop,NP,adj,n,c)

fit = zeros(NP,1);

% ����ÿ�����弴ģ���ص��������ֶ�Ӧ��Qgֵ
for r=1:NP
    W=adj; % ��Ȩ���磬Ȩ�ؾ���W = �ڽӾ���adj
    m=sum(W,2); % n���ڵ�Ķ�
    m2 = sum(sum(W)); % m2=2M=156   
    fit(r)=Qgcpp(pop(:,:,r),adj,m,m2,n,c);
%     fit(r,2)=Qg(pop(:,:,r),adj,n,c);
    
end 
