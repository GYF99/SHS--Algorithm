
%% compute the nuber of edges in the network corresponding to "adj"
function m = numedges(adj)

if issymmetric(adj)   % 若adj是对称矩阵，即adj==adj'，等同于自身的转置矩阵
    
    sl = selfloops(adj);  % adj的主对角线元素求和，0

    m = (sum(sum(adj))-sl)/2 + sl;  % 若邻接矩阵对对称的，则总边数为邻接矩阵元素之和除以2
    
    return

elseif ~issymmetric(adj)
    
    m=sum(sum(adj)); % 若邻接矩阵为非对称的，则总边数为邻接矩阵元素之和
    
    return
    
end