% 文件路径
file_path = 'GRQCbestx.mat';

% 使用load函数加载.mat文件并将内容赋值给bestx变量
load(file_path, 'bestx');

% 现在你可以使用bestx变量访问矩阵数据
% disp(bestx);
 %结构洞节点排序
 % 计算每一列的标准差 std
selected_columns1 = [];
std_values = std(bestx);
% 按照标准差值从小到大进行排序
[sorted_std_values, idx] = sort(std_values);
% 定义前 n 个标准差的数量
top_k = 40;
 % 获取前 n 个标准差的列序号
selected_columns1 = [selected_columns1, idx(1:top_k)];