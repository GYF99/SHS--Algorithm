%% the procedure to compute the modularity Qg values of individuals in pop

function fit = fitness_Qg(pop,NP,adj,n,c)

fit = zeros(NP,1);

% ����ÿ�����弴ģ���ص��������ֶ�Ӧ��Qgֵ
for r=1:NP
    
    fit(r)=Qg(pop(:,:,r),adj,n,c);
    
end 
