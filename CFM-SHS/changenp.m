% 导入文本文件
data = importdata('LFR10000_d100_dmax500_1.txt');

% 保存为.mat文件
save('LFR10000_d100_dmax500_1.mat', 'data');
