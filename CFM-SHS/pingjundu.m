% 读取文本文件
data = load('Wiki-Vote.txt'); % 替换 'your_file.txt' 为你的文件路径

% 构建邻接矩阵
numNodes = max(max(data)); % 获取节点数量
adjMatrix = zeros(numNodes); % 初始化邻接矩阵

% 根据连边信息更新邻接矩阵
for i = 1:size(data, 1)
    node1 = data(i, 1);
    node2 = data(i, 2);
    adjMatrix(node1, node2) = 1;
    adjMatrix(node2, node1) = 1; % 无向网络，所以要对称
end

% 计算节点度数
degrees = sum(adjMatrix, 2); % 对每一行求和，得到每个节点的度数

% 计算平均度
averageDegree = mean(degrees); % 平均度为所有节点度数的平均值

disp(['Average Degree: ', num2str(averageDegree)]);


% 计算中度数大于平均度的节点数
nodesAboveAverage = sum(degrees > averageDegree);

disp(['Number of nodes with degree greater than average: ', num2str(nodesAboveAverage)]);
