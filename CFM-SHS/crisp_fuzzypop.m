function pop_crisp = crisp_fuzzypop( pop, n, NP )
% 根据每个节点的最大隶属度社区标号，构建NP个membership
pop_crisp = [];
for r = 1:NP
    membership = zeros(1,n);
    for i =1:n
        [temp,index]=max(pop(:,i,r));
        membership(i)=index;
    end
    pop_crisp = [pop_crisp; membership];       
end