% 读取数据
data = dlmread('p计算.txt');

% 获取两列数据
column1 = data(:, 1);
column2 = data(:, 2);

% 执行威尔科克森秩和检验
[p_value, ~, stats] = ranksum(column1, column2);

% 输出p值和统计信息
disp(['p值为: ', num2str(p_value)]);
% 判断是否显著
alpha = 0.05; % 显著性水平
if p_value < alpha
    disp('p值小于显著性水平，拒绝原假设，认为两个样本的中位数存在显著差异。');
else
    disp('p值大于显著性水平，接受原假设，认为两个样本的中位数没有显著差异。');
end
