% �ļ�·��
file_path = 'GRQCbestx.mat';

% ʹ��load��������.mat�ļ��������ݸ�ֵ��bestx����
load(file_path, 'bestx');

% ���������ʹ��bestx�������ʾ�������
% disp(bestx);
 %�ṹ���ڵ�����
 % ����ÿһ�еı�׼�� std
selected_columns1 = [];
std_values = std(bestx);
% ���ձ�׼��ֵ��С�����������
[sorted_std_values, idx] = sort(std_values);
% ����ǰ n ����׼�������
top_k = 40;
 % ��ȡǰ n ����׼��������
selected_columns1 = [selected_columns1, idx(1:top_k)];