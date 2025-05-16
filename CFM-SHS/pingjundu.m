% ��ȡ�ı��ļ�
data = load('Wiki-Vote.txt'); % �滻 'your_file.txt' Ϊ����ļ�·��

% �����ڽӾ���
numNodes = max(max(data)); % ��ȡ�ڵ�����
adjMatrix = zeros(numNodes); % ��ʼ���ڽӾ���

% ����������Ϣ�����ڽӾ���
for i = 1:size(data, 1)
    node1 = data(i, 1);
    node2 = data(i, 2);
    adjMatrix(node1, node2) = 1;
    adjMatrix(node2, node1) = 1; % �������磬����Ҫ�Գ�
end

% ����ڵ����
degrees = sum(adjMatrix, 2); % ��ÿһ����ͣ��õ�ÿ���ڵ�Ķ���

% ����ƽ����
averageDegree = mean(degrees); % ƽ����Ϊ���нڵ������ƽ��ֵ

disp(['Average Degree: ', num2str(averageDegree)]);


% �����ж�������ƽ���ȵĽڵ���
nodesAboveAverage = sum(degrees > averageDegree);

disp(['Number of nodes with degree greater than average: ', num2str(nodesAboveAverage)]);
