% ��ȡ����
data = dlmread('p����.txt');

% ��ȡ��������
column1 = data(:, 1);
column2 = data(:, 2);

% ִ�������ƿ�ɭ�Ⱥͼ���
[p_value, ~, stats] = ranksum(column1, column2);

% ���pֵ��ͳ����Ϣ
disp(['pֵΪ: ', num2str(p_value)]);
% �ж��Ƿ�����
alpha = 0.05; % ������ˮƽ
if p_value < alpha
    disp('pֵС��������ˮƽ���ܾ�ԭ���裬��Ϊ������������λ�������������졣');
else
    disp('pֵ����������ˮƽ������ԭ���裬��Ϊ������������λ��û���������졣');
end
