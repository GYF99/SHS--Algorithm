%%%%CCFMM_Qg_LFR实验
clear all
% global net
looptime=1;
Qglist=zeros(1,looptime);
NMIlist=zeros(1,looptime);
mmdrnodes_num=zeros(1,looptime);
mmdrnodes_percent=zeros(1,looptime);
best_in_history_Qgloop=cell(looptime,1);
deltaTimeloop=cell(looptime,1);
% 初始化结果和运行时间的数组
results = zeros(looptime, 1);
runtimes = zeros(looptime, 1);
for iter=1:looptime 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%name
name=['LFR10000_d100_dmax500_1_CFMCMMDR_0501_Qg_NP100_',num2str(iter)];
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % % 网络
% % aaa=['5',num2str(iter-1)];
% aaa='1000';  %节点数
% name_net=['D:\研三上学期+寒假资料\实验代码-新\LFR网络\mu03-runtime\network',aaa,'-adj-',aaa,'.mat'];
% name_community=['D:\研三上学期+寒假资料\实验代码-新\LFR网络\mu03-runtime\community',aaa,'.dat'];
% load(name_net)
% community=load(name_community);
% adj=A;
% groundtruth=community(:,2);   
load('LFR10000_d100_dmax500-adj-10000.mat')
load LFR10000_d100_dmax500_groundtruth.txt  LFR10000_d100_dmax500_groundtruth -ASCII
groundtruth=LFR10000_d100_dmax500_groundtruth;   %记得修改
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
    mmdrnodesnum_gen=zeros(1,Gen);
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
    genbestx=zeros(c,n,Gen);
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
    mmdrnodesnum_gen(g)=n-sum(canyujinhuaNode);
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
    
    results(iter) = Qg; 
    runtimes(iter) = toc;  %结果计时并存储运行时间
    
 save(['SOSFCD实验记录-LFR\',name],'Qglist','NMIlist','best_in_history_Qgloop','deltaTimeloop','mmdrnodes_num','mmdrnodes_percent');     
end
Qgavg=mean(Qglist(:));
StdQg=std(Qglist); %无偏标准差 n-1
MaxQg=max(Qglist(:));
MinQg=min(Qglist(:));
NMIavg=mean(NMIlist(:));
