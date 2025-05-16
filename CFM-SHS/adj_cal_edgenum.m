   clear
%     load ENZYMES_g163_adj.txt
%     adj=ENZYMES_g163_adj;
%     load Ring_adj.txt
%     adj= Ring_adj;
%         load Y_shape_adj.txt
%     adj = Y_shape_adj; 
%     load Ring_K4_adj.txt
%     adj= Ring_K4_adj;
 load('cora-adj-2708.mat')  ;
 adj=A;
sl = selfloops(adj);  % adj的主对角线元素求和，0

    m = (sum(sum(adj))-sl)/2 + sl;  %