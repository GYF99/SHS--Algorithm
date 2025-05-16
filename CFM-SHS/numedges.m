
%% compute the nuber of edges in the network corresponding to "adj"
function m = numedges(adj)

if issymmetric(adj)   % ��adj�ǶԳƾ��󣬼�adj==adj'����ͬ�������ת�þ���
    
    sl = selfloops(adj);  % adj�����Խ���Ԫ����ͣ�0

    m = (sum(sum(adj))-sl)/2 + sl;  % ���ڽӾ���ԶԳƵģ����ܱ���Ϊ�ڽӾ���Ԫ��֮�ͳ���2
    
    return

elseif ~issymmetric(adj)
    
    m=sum(sum(adj)); % ���ڽӾ���Ϊ�ǶԳƵģ����ܱ���Ϊ�ڽӾ���Ԫ��֮��
    
    return
    
end