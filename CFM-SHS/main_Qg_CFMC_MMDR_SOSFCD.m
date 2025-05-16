%%%%CCFMM_Qg
clear all
% global net
looptime=1;
Qglist=zeros(1,looptime);
NMIlist=zeros(1,looptime);
mmdrnodes_num=zeros(1,looptime);
mmdrnodes_percent=zeros(1,looptime);
best_in_history_Qgloop=cell(looptime,1);
deltaTimeloop=cell(looptime,1);
selected_columns1 = [];
selected_columns2 = [];
selected_columns3 = [];
selected_columns4 = [];
% 初始化结果和运行时间的数组
results = zeros(looptime, 1);
runtimes = zeros(looptime, 1);
for iter=1:looptime 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% name
name=['CFMCMMDR_0414_Qg_NP10000_',num2str(iter)];
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % % 网络
% 需要操作的地方，把以下三行中的每个mu05改成mu06/07/08 对应μ=0.5/0.6/0.7/0.8  
%     load('LFR网络-操作有效性验证\mu08_N1000\LFRnetwork-adj-1000.mat')
%     community=load('LFR网络-操作有效性验证\mu08_N1000\community.dat');
%     groundtruth=community(:,2);    
%     name=['0813_CCFMM_Qg_mu08_',num2str(iter)];
%     adj=A;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%真实网络
%     load('Karate-adj-34.mat')
%     load Karate_groundtruth.txt  Karate_groundtruth -ASCII
%     groundtruth=Karate_groundtruth;   %记得修改
%     adj=A;
% %     
%     load('macaque.mat')
%     load macaque1_groundtruth.txt  macaque1_groundtruth -ASCII
%     groundtruth=macaque1_groundtruth;   %记得修改
%     adj=A;


%     load('brain47-adj-47.mat')
%     load brain47_groundtruth.txt  brain47_groundtruth -ASCII
%     groundtruth=brain47_groundtruth;   %记得修改
%     name=['0814_CCFMM_Qg_brain47_',num2str(iter)];
%     adj=A;
%     
%     load('Dolphins-adj-62.mat')
%     load dolphins_groundtruth.txt  dolphins_groundtruth -ASCII
%     groundtruth=dolphins_groundtruth;   %记得修改
%     name=['0813_CCFMM_Qg_Dolphins_',num2str(iter)];
%     adj=A;
    
%     load('PolBooks-adj-105.mat')
%     load polbooks_groundtruth.txt  polbooks_groundtruth -ASCII
%     groundtruth=polbooks_groundtruth;   %记得修改
%     name=['0814_CCFMM_Qg_Polbooks_',num2str(iter)];
%     adj=A;
% 
%     load('football-a-115.mat')
%     load football_groundtruth.txt  football_groundtruth -ASCII
%     groundtruth=football_groundtruth;   %记得修改
%     adj=A;
% %     
%     load('netscience-connect-adj-379.mat')
%     load('netscience_groundtruth-connect379-25.mat')
%     adj=A;
    
%     load('cornell-adj-195.mat')
%     load cornell_groundtruth.txt  cornell_groundtruth -ASCII
%     groundtruth=cornell_groundtruth;   %记得修改
%     name=['081802_CCFMM_Qg_cornell_',num2str(iter)];
%     adj=A;

%     load('email-Eu-core-adj-1005.mat')
%     load emailEucore_groundtruth.txt emailEucore_groundtruth -ASCII
%     groundtruth=emailEucore_groundtruth;   %记得修改
%     adj=A;
    
%     load('NetScience-adj-1589.mat');
%     load netscience_groundtruth.txt  netscience_groundtruth -ASCII
%     groundtruth=netscience_groundtruth;   %记得修改    
%     adj=A;
%     
%     load('cora-adj-2708.mat');
%     load cora_groundtruth.txt  cora_groundtruth -ASCII
%     groundtruth=cora_groundtruth;   %记得修改
%     
% %     load('cora_comm0.4645.mat')
%     adj=A;
    

%     load('Yeast-adj-2361.mat');
%     load Yeast_groundtruth.txt  Yeast_groundtruth -ASCII
%     groundtruth=Yeast_groundtruth;   %记得修改    
%     adj=A;
% 
%     load('Dblp_adj.mat')
%     load Dblp_groundtruth.txt  Dblp_groundtruth-ASCII
%     groundtruth=Dblp_groundtruth;   %记得修改
%     adj=full(A);

%     load wisconsin_adj.txt
%     adj=wisconsin_adj;
%     load wisconsin_groundtruth.txt  wisconsin_groundtruth -ASCII
%     groundtruth=wisconsin_groundtruth;   %记得修改
%     name=['wisconsin_CCFMM_0824_Runtime_Qg_NP10_',num2str(iter)];
%     
%     load ENZYMES_g163_adj.txt
%     adj=ENZYMES_g163_adj;
%     load ENZYMES_g163_groundtruth.txt  ENZYMES_g163_groundtruth -ASCII
%     groundtruth=ENZYMES_g163_groundtruth;   %记得修改
%     name=['0816_CCFMM_Qg_Enzymes_g163_',num2str(iter)];

%     load('GR-QC.mat')
%     load GRQC_groundtruth.txt  GRQC_groundtruth -ASCII
%     groundtruth=GRQC_groundtruth;   %记得修改
%     adj=A;

%     load('Wiki-Vote.mat')
%     load WikiVote_groundtruth.txt  WikiVote_groundtruth -ASCII
%     groundtruth=WikiVote_groundtruth;   %记得修改
%     adj=A;
% load('LFR10000_03-adj-10000.mat')
% load LFR10000_03_groundtruth.txt  LFR10000_03_groundtruth -ASCII
% groundtruth=LFR10000_03_groundtruth;   %记得修改
% adj=full(A);
% load('LesMis-adj-9.mat')
% load LesMis_groundtruth.txt  LesMis_groundtruth -ASCII
% groundtruth=LesMis_groundtruth;
% adj=A;

%     load('test1-adj-9.mat')
%     load test1_groundtruth.txt  test1_groundtruth -ASCII
%     groundtruth=test1_groundtruth;   %记得修改
%     adj=A;
load('LFR02-adj-1000.mat')
load LFR02_groundtruth.txt  LFR02_groundtruth -ASCII
groundtruth=LFR02_groundtruth;   %记得修改
adj=A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    g=1;
    p=0;
    edge=numedges(adj);
    n=numnodes(adj);
    m=numedges(adj);
    c=max(groundtruth(:));
    threshold_value=0.25;
    Gen = 500;
    mm=sum(adj,1);
    m2 = sum(sum(adj)); %||W||
    B=adj-(((mm').*(mm))./m2);
    NP = 100;
    deltaQtime=[];
    mmdrtime=[];
    unchangenode_Num=zeros(1,Gen);
    best_in_history_Qg=[];
    boundaryNode=culBoundaryNode(adj,groundtruth);
    Nodelist=(1:n);
    IntraNode=setdiff(Nodelist,boundaryNode);
    sigma=[1e-5,1e-3];
    canyujinhuaNode=ones(1,n);
    bubiancishu=zeros(1,n);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%初始化    
    pop = inital_Brain_pop_consv4(n, c, NP,groundtruth,adj);
    fit_Qg = fitness_Qgcpp(pop,NP,adj,n,c);

    bias_pop = bias_init_brain_pop_consv4( pop, c, n, NP, adj, groundtruth);
    bias_fit_Qg = fitness_Qgcpp(bias_pop,NP,adj,n,c);

    win_index = find(bias_fit_Qg > fit_Qg); 
    win_number = length(win_index);
    pop(:,:,win_index) = bias_pop(:,:,win_index);
    fit_Qg(win_index) = bias_fit_Qg(win_index);
    [bestfit,index]=max(fit_Qg); 
    bestx=pop(:,:,randsrc(1,1,find(fit_Qg==max(fit_Qg))'));
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % 3. Main loop
    pop_best_history=zeros(c,n,Gen); 

    while p < 1
    % 3.1 Mutualism 
    sigma_intra=ceil((1/4*Gen+(3/4*Gen-1/4*Gen)*(2-exp(g/Gen*log(2))))*0.1);
    sigma_bridge=ceil(1.3*sigma_intra);
    tstart1=tic;
    Mutu_pop=pop;
    Mutu_fit=fit_Qg;
    for i=1:NP
        % (1)Xbest
        [bestfit,index]=max(Mutu_fit); 
        bestx=Mutu_pop(:,:,randsrc(1,1,find(Mutu_fit==max(Mutu_fit))'));
        seed=setdiff(randperm(NP),i); 
        j=randsrc(1,1,seed);
        mutual_vector=0.5.*(Mutu_pop(:,:,i)+Mutu_pop(:,:,j));
        BF1=round(1+rand); 
        BF2=round(1+rand);
        Xinew=Mutu_pop(:,:,i)+rand(c,n).*(bestx-BF1.*mutual_vector); 
        Xjnew=Mutu_pop(:,:,j)+rand(c,n).*(bestx-BF2.*mutual_vector);
        bucanyu=find(canyujinhuaNode==0);
        if size(bucanyu,2)>0
            for v=1:size(bucanyu,2)
                Xinew(:,bucanyu(v))=Mutu_pop(:,bucanyu(v),i);
                Xjnew(:,bucanyu(v))=Mutu_pop(:,bucanyu(v),j);
            end
        end
        Xinew=bound_brain_SOSFCD_consv4(Xinew,c,n,groundtruth,adj);
        Xjnew=bound_brain_SOSFCD_consv4(Xjnew,c,n,groundtruth,adj);

        canyu=find(canyujinhuaNode==1);
        differindex=canyu;
        differindex2=canyu;              
        Xinew_Qg=deltaQgcpp(Mutu_pop(:,:,i),Xinew,adj,c,differindex,Mutu_fit(i),m2,B);
        Xjnew_Qg = deltaQgcpp(Mutu_pop(:,:,j),Xjnew,adj,c,differindex2,Mutu_fit(j),m2,B);  

        if Xinew_Qg > Mutu_fit(i)
           Mutu_pop(:,:,i)=Xinew;
           Mutu_fit(i)=Xinew_Qg;
        end
        if Xjnew_Qg > Mutu_fit(j)
           Mutu_pop(:,:,j)=Xjnew;
           Mutu_fit(j)=Xjnew_Qg;
        end  
    end

    % 3.2 Commensalism 
    Comm_pop=Mutu_pop;
    Comm_fit=Mutu_fit;
    for i=1:NP
        % (1)Xbest
        [bestfit,index]=max(Comm_fit); 
        bestx=Comm_pop(:,:,randsrc(1,1,find(Comm_fit==max(Comm_fit))'));
        seed=setdiff(randperm(NP),i); 
        j=randsrc(1,1,seed);
        Xinew=Comm_pop(:,:,i)+(rand(c,n)*2-1).*(bestx-Comm_pop(:,:,j));
        bucanyu=find(canyujinhuaNode==0);
        if size(bucanyu,2)>0
            for v=1:size(bucanyu,2)
                Xinew(:,bucanyu(v))=Comm_pop(:,bucanyu(v),i);
            end
        end
        Xinew=bound_brain_SOSFCD_consv4(Xinew,c,n,groundtruth,adj);


        canyu=find(canyujinhuaNode==1);
        differindex=canyu;
        Xinew_Qg = deltaQgcpp(Comm_pop(:,:,i),Xinew,adj,c,differindex,Comm_fit(i),m2,B);

        if Xinew_Qg > Comm_fit(i)
           Comm_pop(:,:,i)=Xinew;
           Comm_fit(i)=Xinew_Qg;
        end
    end 

    % 3.5 Parasitism  
    Para_pop=Comm_pop;
    Para_fit=Comm_fit;
    for i=1:NP
        seed=setdiff(randperm(NP),i); 
        j=randsrc(1,1,seed);
        Para_vector=Para_pop(:,:,i);
        canyuNodeNum=sum(canyujinhuaNode);
        canyuNode=find(canyujinhuaNode==1);
        seed=randperm(canyuNodeNum);
        pick=canyuNode(seed(1:ceil(rand*canyuNodeNum)));   
        Para_vector(:,pick) = inital_Brain_pop(length(pick), c, 1,groundtruth);
        bucanyu=find(canyujinhuaNode==0);
        if size(bucanyu,2)>0
            for v=1:size(bucanyu,2)
                Para_vector(:,bucanyu(v))=Para_pop(:,bucanyu(v),j);
            end
        end
        Para_vector = bound_brain_SOSFCD_consv4(Para_vector,c,n,groundtruth,adj);

        canyu=find(canyujinhuaNode==1);
        differindex=canyu;
        Para_vector_Qg = deltaQgcpp(Para_pop(:,:,j),Para_vector,adj,c,differindex,Para_fit(j),m2,B);

        
        if Para_vector_Qg > Para_fit(j)
           Para_pop(:,:,j)=Para_vector;
           Para_fit(j)=Para_vector_Qg;
        end
    end 

    % 3.6 Parasitism clean_uper  
    clean_pop = Para_pop;
    clean_fit = Para_fit;
    clean_pop_copy=clean_pop;
    clean_fit_copy=clean_fit;
    clean_pop_before = crisp_fuzzypop(clean_pop, n, NP );  
    clean_pop_after = clean_up_random( clean_pop_before, n, NP, adj, threshold_value );

    for r=1:NP
        for node = 1:n
            if clean_pop_after(r,node) ~= clean_pop_before(r,node)
                % node节点的社区划分调整了
                k = clean_pop_after(r,node);
                clean_pop(k,node,r) = clean_pop(k,node,r) + 0.5;
            end
        end
        bucanyu=find(canyujinhuaNode==0);
        if size(bucanyu,2)>0
            for v=1:size(bucanyu,2)
                clean_pop(:,bucanyu(v),r)=clean_pop_copy(:,bucanyu(v),r);
            end
        end        
        clean_pop(:,:,r) = bound_brain_SOSFCD_consv4(clean_pop(:,:,r),c,n,groundtruth,adj);


        canyu=find(canyujinhuaNode==1);
        differindex=canyu;
        clean_fit(r) = deltaQgcpp(clean_pop_copy(:,:,r),clean_pop(:,:,r),adj,c,differindex,clean_fit_copy(r),m2,B);

        if clean_fit(r) > Para_fit(r)
           Para_pop(:,:,r) = clean_pop(:,:,r);
           Para_fit(r) = clean_fit(r);
        end        
    end

    % 3.7 updata pop and fit
    pop = Para_pop;
    fit_Qg = Para_fit;
    [best_Qg,index]=max(fit_Qg); 
    bestx=pop(:,:,index);
    genbestx(:,:,g)=pop(:,:,index);
    genbestQ(g)=best_Qg; 
    best_in_history_Qg(g)=best_Qg; 
    unchangeindex=[];

    if g>1
    tstart2=tic;
    Ubest_before=genbestx(:,:,g-1);
    canyu=find(canyujinhuaNode==1);
    for q=1:size(canyu,2)
        difference=std(genbestx(:,canyu(q),g));
        delta=(sigma(2)-sigma(1))*difference+sigma(1);
        if abs(Ubest_before(:,canyu(q))-genbestx(:,canyu(q),g))<delta
            bubiancishu(canyu(q))=bubiancishu(canyu(q))+1;
        elseif abs(Ubest_before(:,canyu(q))-genbestx(:,canyu(q),g))>sigma(2)         
            bubiancishu(canyu(q))=0;
        end

    end

    if numel(intersect(find(bubiancishu>sigma_bridge),canyu))>0
        for q=1:size(canyu,2)
            if bubiancishu(canyu(q))>sigma_bridge && ismember(canyu(q),boundaryNode)==1
                unchangeindex=[unchangeindex canyu(q)];
            elseif bubiancishu(canyu(q))>sigma_intra && ismember(canyu(q),boundaryNode)==0
                unchangeindex=[unchangeindex canyu(q)];
            end
        end
        unchangenode_Num(g)=size(unchangeindex,2);
        for q=1:size(unchangeindex,2)
            i=unchangeindex(q);
            for r=1:NP
                pop(:,i,r)=bestx(:,i);
            end
            canyujinhuaNode(i)=0;
        end
        fit_Qg = fitness_Qgcpp(pop,NP,adj,n,c) ;     
    end
    mmdrtime=[mmdrtime toc(tstart2)];
    end

    if g>1
        deltaQtime(g)=deltaQtime(g-1)+toc(tstart1);
    else
        deltaQtime(g)=toc(tstart1);
    end
    % 3.11 whether the loop stop?
    if ~mod(g,1)
       clc
       fprintf('%d generations completed\n',g);
    end

    if (g>=Gen) ||(sum(canyujinhuaNode)==0) 
       p=1;
     elseif g>30&&(abs(best_in_history_Qg(g)-best_in_history_Qg(g-30))<1e-8)
        p=1;
    end
    if p==0
       g=g+1;  
    end
    
    end 
    fprintf('\n');
   
    Qg = max(best_in_history_Qg(g));
    Qglist(iter)=Qg;
    deltaTimeloop{iter}(:,:)=deltaQtime;
    best_in_history_Qgloop{iter}(:,:)=best_in_history_Qg;
    mmdrnodes_num(iter)=n-sum(canyujinhuaNode);
    mmdrnodes_percent(iter)=mmdrnodes_num(iter)/n;
    
    
    results(iter) = Qg;
    runtimes(iter) = toc;  %结果计时并存储运行时间
    
    
    


    cdlist=zeros(1,n);
    %检查groundtruth是否正确
    correct_node = 0;
    for i=1:n
        [max_value,max_index]=max(bestx(:,i));
        cdlist(i)=max_index;
        if max_index == groundtruth(i)
            correct_node=correct_node+1;
        end
    end

    real_cdlist=zeros(1,n);
    for i=1:n
        real_cdlist(i)=groundtruth(i);
    end

    nmi=NMI(real_cdlist,cdlist);
    NMIlist(iter)=nmi;
    
    %结构洞节点排序
    % 计算每一列的标准差 std
    std_values = std(bestx);
    % 按照标准差值从小到大进行排序
    [sorted_std_values, idx] = sort(std_values);
    % 定义前 n 个标准差的数量
    top_k = 9;
     % 获取前 n 个标准差的列序号
    selected_columns1 = [selected_columns1, idx(1:top_k)];
% 
%     % 计算每列的最小值，然后取平方根 k-overlap
%     min_values = sqrt(min(bestx));
%     % 排序结果
%     [sorted_values, sorted_indices] = sort(min_values, 'descend');
%     top_k = 20;
%    % 获取前 n 个列序号
%     selected_columns2 = [selected_columns2, sorted_indices(1:top_k)];
%     
% 
%     % 计算每一列中所有数相乘除以每一列数相乘后相加的结果  B指标
%     column_products = prod(bestx, 1);
%     sum_of_products = sum(column_products);
%     column_ratios1 = column_products ./ sum_of_products;
%     % 将结果按从大到小排序
%     [sorted_ratios1, indices] = sort(column_ratios1, 'descend');
%      top_k = 20;
%     %      % 获取前 n 个列序号
%     selected_columns3 = [selected_columns3, indices(1:top_k)];
% 
%     
%     % 获取每一列中次大值除以最大值的结果
%     column_max_values = max(bestx, [], 1);  % 每列的最大值
%     column_sorted = sort(bestx, 'descend');  % 对每列进行降序排序
%     column_second_max_values = column_sorted(2, :);  % 每列的次大值
%     column_ratios2 = column_second_max_values ./ column_max_values;
%     % 将结果按从大到小排序
%     [sorted_ratios2, sorted_idx] = sort(column_ratios2, 'descend');
%     top_k = 20;
%     % 获取前 n 个列序号
%     selected_columns4 = [selected_columns4, sorted_idx(1:top_k)];
    
%  save(['SOSFCD实验记录\',name],'Qglist','NMIlist','best_in_history_Qgloop','deltaTimeloop','mmdrnodes_num','mmdrnodes_percent');     
end
Qgavg=mean(Qglist(:));
StdQg=std(Qglist); %无偏标准差 n-1
MaxQg=max(Qglist(:));
MinQg=min(Qglist(:));
NMIavg=mean(NMIlist(:));

% % 绘制直方图
% % 创建画布，并绘制两个子图
% figure;
% 
% % 第一个子图：运行时间的直方图
% subplot(2, 1, 1);
% histogram(runtimes);
% title('Algorithm Runtime Distribution');
% xlabel('Runtime (s)');
% ylabel('Frequency');
% 
% % 第二个子图：结果的直方图
% subplot(2, 1, 2);
% histogram(results);
% title('Algorithm Results Distribution');
% xlabel('Results');
% ylabel('Frequency');

