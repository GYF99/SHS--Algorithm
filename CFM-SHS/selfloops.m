
%% compute the number of selfloops in the network corresponding to "adj"

function sl=selfloops(adj)

sl = sum( diag(adj) );  % diag���Խ���Ԫ�����

