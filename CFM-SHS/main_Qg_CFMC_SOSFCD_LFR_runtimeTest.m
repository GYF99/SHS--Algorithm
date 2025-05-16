%%%%% CFMC_Qg_LFR实验  runtimeTest  nmm只在后30%的迭代中使用  Qg用cpp版本
clc
clear
%运行次数
looptime=1000;
Qglist=zeros(1,looptime);
NMIlist=zeros(1,looptime);
best_in_history_Qgloop=cell(looptime,1);
deltaTimeloop=cell(looptime,1);
deltaTimeloopsum=zeros(1,looptime);
% 初始化结果和运行时间的数组
results = zeros(looptime, 1);
runtimes = zeros(looptime, 1);
for iter=1:looptime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%name
% name=['N100_LFR3_CFMC_1120_Qg_NP100_',num2str(iter)];
name=['football_CFMC_0310_Qg_NP100_',num2str(iter)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % % LFR网络
% % aaa=['3',num2str(iter-1)];
% aaa='5000';  %节点数
% name_net=['D:\研三上学期+寒假资料\实验代码-新\LFR网络\mu03-runtime\network',aaa,'-adj-',aaa,'.mat'];
% name_community=['D:\研三上学期+寒假资料\实验代码-新\LFR网络\mu03-runtime\community',aaa,'.dat'];
% load(name_net)
% community=load(name_community);
% adj=A;
% groundtruth=community(:,2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     load('Karate-adj-34.mat')
%     load Karate_groundtruth.txt  Karate_groundtruth -ASCII
%     groundtruth=Karate_groundtruth;   %记得修改
%     adj=A;
%     load('football-a-115.mat')
%     load football_groundtruth.txt  football_groundtruth -ASCII
%     groundtruth=football_groundtruth;   %记得修改
%     adj=A;
%     load('email-Eu-core-adj-1005.mat')
%     load emailEucore_groundtruth.txt emailEucore_groundtruth -ASCII
%     groundtruth=emailEucore_groundtruth;   %记得修改
%     adj=A;
%     load('NetScience-adj-1589.mat');
%     load netscience_groundtruth.txt  netscience_groundtruth -ASCII
%     groundtruth=netscience_groundtruth;   %记得修改    
%     adj=A;
%     load('Yeast-adj-2361.mat');
%     load Yeast_groundtruth.txt  Yeast_groundtruth -ASCII
%     groundtruth=Yeast_groundtruth;   %记得修改    
%     adj=A;
%     load('cora-adj-2708.mat');
%     load cora_groundtruth.txt  cora_groundtruth -ASCII
%     groundtruth=cora_groundtruth;   %记得修改
%     adj=A;
%     load('Dblp_adj.mat')
%     load Dblp_groundtruth.txt  Dblp_groundtruth-ASCII
%     groundtruth=Dblp_groundtruth;   %记得修改
%     adj=full(A);
tic;
load('LFR100-adj-100.mat')
load LFR100_groundtruth.txt  LFR100_groundtruth -ASCII
groundtruth=LFR100_groundtruth;   %记得修改
% adj=A;
adj=A;

    % 1. Initialization
    g=1;
    p=0;
    edge=numedges(adj);
    n=numnodes(adj);
    c=max(groundtruth(:));
    threshold_value=0.25;
    Gen = 500;
    NP = 100;
    mm=sum(adj,1);
    m2 = sum(sum(adj)); %||W||
    B=adj-(((mm').*(mm))./m2);
    deltaQtime=[];
    mututime=zeros(1,Gen);
    commtime=zeros(1,Gen);
    paratime=zeros(1,Gen);
    nmmtime=zeros(1,Gen);
    cfmctimemutu=zeros(1,NP);
    runqgtimemutu=zeros(1,NP);
    cfmctimecomm=zeros(1,NP);
    runqgtimecomm=zeros(1,NP);
    cfmctimepara=zeros(1,NP);
    runqgtimepara=zeros(1,NP);
    cfmctimenmm=zeros(1,NP);
    runqgtimenmm=zeros(1,NP);
    
    best_in_history_Qg=[];
    % 1.3 pop construction
    pop = inital_Brain_pop_consv4(n, c, NP,groundtruth,adj);
%     %初始化时把社区内部节点的增量模块度算出来
%     boundaryNode=culBoundaryNode(adj,groundtruth);    
%     Nodelist=(1:n);
%     IntraNode=setdiff(Nodelist,boundaryNode);
%     for r=1:NP
%         for i=IntraNode   
%             pop(:,i,r)=zeros(c,1);
%             pop(groundtruth(i),i,r)=1;
%         end 
%     end
    fit_Qg = fitness_Qgcpp(pop,NP,adj,n,c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%有偏操作
%     %2. the biased process used in the initialization
%     bias_pop = bias_init_brain_pop_consv4( pop, c, n, NP, adj, groundtruth);
%     bias_fit_Qg = fitness_Qgcpp(bias_pop,NP,adj,n,c);
    
%     win_index = find(bias_fit_Qg > fit_Qg); 
%     win_number = length(win_index);
%     pop(:,:,win_index) = bias_pop(:,:,win_index);
%     fit_Qg(win_index) = bias_fit_Qg(win_index);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    [best_Qg,index]=max(fit_Qg); 
    bestx=pop(:,:,randsrc(1,1,find(fit_Qg==max(fit_Qg))'));
    
    
    % 3. Main loop
    pop_best_history=zeros(c,n,Gen); 
    
    % do the evolution iteratively
    while p < 1
    totaltime=tic;
    mutu=tic;
    % 3.1 Mutualism
    Mutu_pop=pop;
    Mutu_fit=fit_Qg;
    better_number=0;
    for i=1:NP
        % (1)Xbest
        [bestfit,~]=max(Mutu_fit);   
        bestx=Mutu_pop(:,:,randsrc(1,1,find(Mutu_fit==max(Mutu_fit))'));
        seed=setdiff(randperm(NP),i); 
        j=randsrc(1,1,seed);
        mutual_vector=0.5.*(Mutu_pop(:,:,i)+Mutu_pop(:,:,j));
        BF1=round(1+rand); 
        BF2=round(1+rand);
        % (4)Xinew & Xjnew
        Xinew=Mutu_pop(:,:,i)+rand(c,n).*(bestx-BF1.*mutual_vector); 
        Xjnew=Mutu_pop(:,:,j)+rand(c,n).*(bestx-BF2.*mutual_vector);
        cfmcmutu=tic;
        Xinew=bound_brain_SOSFCD_consv4(Xinew,c,n,groundtruth,adj);
        Xjnew=bound_brain_SOSFCD_consv4(Xjnew,c,n,groundtruth,adj);
        
        if i>1
            cfmctimemutu(i)=cfmctimemutu(i-1)+toc(cfmcmutu);
        else
            cfmctimemutu(i)=toc(cfmcmutu);
        end
        runqgmutu=tic;
        differindex=1:n;
        Xinew_Qg_node=deltaQg_node(bestx,Xinew,adj,c,differindex,bestfit,m2,B);
        Xjnew_Qg_node=deltaQg_node(bestx,Xjnew,adj,c,differindex,bestfit,m2,B);
        Xinew_badnodes=find(Xinew_Qg_node<0);
        Xjnew_badnodes=find(Xjnew_Qg_node<0);
        Xinew(:,Xinew_badnodes)=bestx(:,Xinew_badnodes);
        Xjnew(:,Xjnew_badnodes)=bestx(:,Xjnew_badnodes);
        Xinew_Qg_node(Xinew_badnodes)=0;
        Xjnew_Qg_node(Xjnew_badnodes)=0;
        Xinew_Qg=bestfit+sum(Xinew_Qg_node);
        Xjnew_Qg=bestfit+sum(Xjnew_Qg_node);

        if i>1
            runqgtimemutu(i)=runqgtimemutu(i-1)+toc(runqgmutu);
        else
            runqgtimemutu(i)=toc(runqgmutu);
        end
        if Xinew_Qg > Mutu_fit(i)
           Mutu_pop(:,:,i)=Xinew;
           Mutu_fit(i)=Xinew_Qg;
           better_number=better_number+1;
        end
        if Xjnew_Qg > Mutu_fit(j)
           Mutu_pop(:,:,j)=Xjnew;
           Mutu_fit(j)=Xjnew_Qg;
           better_number=better_number+1;
        end  
    end
    mututime(g)=toc(mutu);
    comm=tic;
    % 3.2 Commensalism
    Comm_pop=Mutu_pop;
    Comm_fit=Mutu_fit;
    better_number=0;
    for i=1:NP
        % (1)Xbest
        [bestfit,~]=max(Comm_fit); 
        bestx=Comm_pop(:,:,randsrc(1,1,find(Comm_fit==max(Comm_fit))'));
        % (2)Xj~=Xi
        seed=setdiff(randperm(NP),i); 
        j=randsrc(1,1,seed);
        % (3)Xinew 
        Xinew=Comm_pop(:,:,i)+(rand(c,n)*2-1).*(bestx-Comm_pop(:,:,j));
        cfmccomm=tic;
        Xinew=bound_brain_SOSFCD_consv4(Xinew,c,n,groundtruth,adj);
        if i>1
            cfmctimecomm(i)=cfmctimecomm(i-1)+toc(cfmccomm);
        else
            cfmctimecomm(i)=toc(cfmccomm);
        end
        runqgcomm=tic;
        differindex=1:n;
        Xinew_Qg_node=deltaQg_node(bestx,Xinew,adj,c,differindex,bestfit,m2,B);
        Xinew_badnodes=find(Xinew_Qg_node<0);
        Xinew(:,Xinew_badnodes)=bestx(:,Xinew_badnodes);
        Xinew_Qg_node(Xinew_badnodes)=0;
        Xinew_Qg=bestfit+sum(Xinew_Qg_node);
        if i>1
            runqgtimecomm(i)=runqgtimecomm(i-1)+toc(runqgcomm);
        else
            runqgtimecomm(i)=toc(runqgcomm);
        end
        if Xinew_Qg > Comm_fit(i)
           Comm_pop(:,:,i)=Xinew;
           Comm_fit(i)=Xinew_Qg;
           better_number=better_number+1;
        end
    end 
    commtime(g)=toc(comm);
    % 3.5 Parasitism  / better_number = 4
    para=tic;
    Para_pop=Comm_pop;
    Para_fit=Comm_fit;
    better_number=0;
    for i=1:NP
        % (1)Xbest
        [bestfit,~]=max(Para_fit); 
        bestx=Para_pop(:,:,randsrc(1,1,find(Para_fit==max(Para_fit))'));
        % (1)Xj~=Xi
        seed=setdiff(randperm(NP),i); 
        j=randsrc(1,1,seed);
        % (2)Parasite Vector
        Para_vector=Para_pop(:,:,i);
        seed=randperm(n);
        pick=seed(1:ceil(rand*n));   
        Para_vector(:,pick) = inital_Brain_pop(length(pick), c, 1,groundtruth);
        cfmcpara=tic;
        Para_vector = bound_brain_SOSFCD_consv4(Para_vector,c,n,groundtruth,adj);
        if i>1
            cfmctimepara(i)=cfmctimepara(i-1)+toc(cfmcpara);
        else
            cfmctimepara(i)=toc(cfmcpara);
        end
        runqgpara=tic;
        differindex=1:n;
        Para_vector_Qg_node=deltaQg_node(bestx,Para_vector,adj,c,differindex,bestfit,m2,B);
        Para_vector_badnodes=find(Para_vector_Qg_node<0);
        Para_vector(:,Para_vector_badnodes)=bestx(:,Para_vector_badnodes);
        Para_vector_Qg_node(Para_vector_badnodes)=0;
        Para_vector_Qg=bestfit+sum(Para_vector_Qg_node);
        
        if i>1
            runqgtimepara(i)=runqgtimepara(i-1)+toc(runqgpara);
        else
            runqgtimepara(i)=toc(runqgpara);
        end
        %FE=FE+1;        
        % (3)updata Para_pop and Para_fit
        if Para_vector_Qg > Para_fit(j)
           Para_pop(:,:,j)=Para_vector;
           Para_fit(j)=Para_vector_Qg;
           better_number=better_number+1;
        end
    end 
    paratime(g)=toc(para);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3.6 Parasitism clean_uper  
    if g>(0.7*Gen)
        nmm=tic;
        clean_pop = Para_pop;
        clean_fit = Para_fit;
        clean_pop_before = crisp_fuzzypop(clean_pop, n, NP );  % NP个groundtruth？？
        clean=tic;
        clean_pop_after = clean_up_random( clean_pop_before, n, NP, adj, threshold_value );
        cleantime=toc(clean);

        for r=1:NP
            [bestfit,~]=max(clean_fit); 
            bestx=clean_pop(:,:,randsrc(1,1,find(clean_fit==max(clean_fit))'));
            for node = 1:n
                if clean_pop_after(r,node) ~= clean_pop_before(r,node)
                    k = clean_pop_after(r,node);
                    clean_pop(k,node,r) = clean_pop(k,node,r) + 0.5;
                end
            end
            cfmcnmm=tic;
            clean_pop(:,:,r) = bound_brain_SOSFCD_consv4(clean_pop(:,:,r),c,n,groundtruth,adj); %又修回
            if r>1
                cfmctimenmm(r)=cfmctimenmm(r-1)+toc(cfmcnmm);
            else
                cfmctimenmm(r)=toc(cfmcnmm);
            end
            runqgnmm=tic;
            differindex=1:n;
            clean_Qg_node=deltaQg_node(bestx,clean_pop(:,:,r),adj,c,differindex,bestfit,m2,B);
            clean_badnodes=find(clean_Qg_node<0);
            clean_pop(:,clean_badnodes,r)=bestx(:,clean_badnodes);
            clean_Qg_node(clean_badnodes)=0;
            clean_fit(r)=bestfit+sum(clean_Qg_node);
            if r>1
                runqgtimenmm(r)=runqgtimenmm(r-1)+toc(runqgnmm);
            else
                runqgtimenmm(r)=toc(runqgnmm);
            end
        end


        for r =1:NP
            if clean_fit(r) > Para_fit(r)
               Para_pop(:,:,r) = clean_pop(:,:,r);
               Para_fit(r) = clean_fit(r);
            end
        end
        nmmtime(g)=toc(nmm);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3.7 updata pop and fit
    pop = Para_pop;
    fit_Qg = Para_fit;
    % 3.8 recode bestfit and bestx
    [best_Qg,index]=max(fit_Qg); 
    bestx=pop(:,:,index);
    pop_best_history(:,:,g)=bestx;
    best_in_history_Qg(g)=best_Qg; 
    
    if g>1
        deltaQtime(g)=deltaQtime(g-1)+toc(totaltime);
    else
        deltaQtime(g)=toc(totaltime);
    end
    % 3.11 whether the loop stop?
    if ~mod(g,1)
       clc
       fprintf('%d generations completed\n',g);
    end
    
    if g>=Gen
       p=1;
   elseif g>30&&(abs(best_in_history_Qg(g)-best_in_history_Qg(g-30))<1e-8)
        p=1;
    end
    if p==0
       g=g+1;  
    end

    end % while p<1
    fprintf('\n');


    Qg = max(best_in_history_Qg(g))
    Qglist(iter)=Qg;
    deltaTimeloop{iter}(:,:)=deltaQtime;
    deltaTimeloopsum(iter)=deltaQtime(g-30);
    best_in_history_Qgloop{iter}(:,:)=best_in_history_Qg;
    cdlist=zeros(1,n);
    
    
    results(iter) = Qg;
    runtimes(iter) = toc;  %结果计时并存储运行时间
    
    
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
% save(['SOSFCD实验记录-runtime-improved\',name],'Qglist','NMIlist','best_in_history_Qgloop','deltaTimeloop','deltaTimeloopsum');    
end
Qgavg=mean(Qglist(:));
StdQg=std(Qglist); %无偏标准差 n-1
MaxQg=max(Qglist(:));
MinQg=min(Qglist(:));
NMIavg=mean(NMIlist(:));

% 绘制直方图
% 创建画布，并绘制两个子图
figure;

% 第一个子图：运行时间的直方图
subplot(2, 1, 1);
histogram(runtimes);
title('Algorithm Runtime Distribution');
xlabel('Runtime (s)');
ylabel('Frequency');

% 第二个子图：结果的直方图
subplot(2, 1, 2);
histogram(results);
title('Algorithm Results Distribution');
xlabel('Results');
ylabel('Frequency');


