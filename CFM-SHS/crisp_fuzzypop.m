function pop_crisp = crisp_fuzzypop( pop, n, NP )
% ����ÿ���ڵ�����������������ţ�����NP��membership
pop_crisp = [];
for r = 1:NP
    membership = zeros(1,n);
    for i =1:n
        [temp,index]=max(pop(:,i,r));
        membership(i)=index;
    end
    pop_crisp = [pop_crisp; membership];       
end