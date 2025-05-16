%% the procedure to compute the modularity Qg values of individuals in pop

function fit = fitness_Qg(pop,NP,adj,n,c)

fit = zeros(NP,1);

% 计算每个个体即模糊重叠社区划分对应的Qg值
for r=1:NP
    
    fit(r)=Qg(pop(:,:,r),adj,n,c);
    
end 
