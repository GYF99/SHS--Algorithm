%%%%% SimpleMC+MMDR_Qg实验
clear all
% global net
looptime=1;
Qglist=zeros(1,looptime);
NMIlist=zeros(1,looptime);
best_in_history_Qgloop=cell(looptime,1);
deltaTimeloop=cell(looptime,1);
for iter=1:looptime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%name
% name=['Netscience_SimpleMCMMDR_1020_Qg_NP100_',num2str(iter)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % % 网络
% 需要操作的地方，把以下三行中的每个mu05改成mu06/07/08 对应μ=0.5/0.6/0.7/0.8  
%     load('LFR网络-操作有效性验证\mu06_N1000\LFRnetwork-adj-1000.mat')
%     community=load('LFR网络-操作有效性验证\mu06_N1000\community.dat');
%     name=['0806_SimpleMC_MMDR_Qg_mu06_',num2str(iter)];
%     groundtruth=community(:,2);   
%     adj=A;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%真实网络
    load('Karate-adj-34.mat')
    load Karate_groundtruth.txt  Karate_groundtruth -ASCII
    groundtruth=Karate_groundtruth;   %记得修改
    adj=A;
%     
%     load('brain47-adj-47.mat')
%     load brain47_groundtruth.txt  brain47_groundtruth -ASCII
%     groundtruth=brain47_groundtruth;   %记得修改
%     name=['0814_SimpleMC_MMDR_Qg_brain47_',num2str(iter)];
%     adj=A;

%     load('Dolphins-adj-62.mat')
%     load dolphins_groundtruth.txt  dolphins_groundtruth -ASCII
%     groundtruth=dolphins_groundtruth;   %记得修改
%     name=['0814_SimpleMC_MMDR_Qg_dolphins_',num2str(iter)];
%     adj=A;
%     
%     load('PolBooks-adj-105.mat')
%     load polbooks_groundtruth.txt  polbooks_groundtruth -ASCII
%     groundtruth=polbooks_groundtruth;   %记得修改
%     name=['0814_SimpleMC_MMDR_Qg_Polbooks_',num2str(iter)];
%     adj=A;

%     load('football-a-115.mat')
%     load football_groundtruth.txt  football_groundtruth -ASCII
%     groundtruth=football_groundtruth;   %记得修改
%     adj=A;

%     load('email-Eu-core-adj-1005.mat')
%     load emailEucore_groundtruth.txt emailEucore_groundtruth -ASCII
%     groundtruth=emailEucore_groundtruth;   %记得修改
%     adj=A;

%     load('cornell-adj-195.mat')
%     load cornell_groundtruth.txt  cornell_groundtruth -ASCII
%     groundtruth=cornell_groundtruth;   %记得修改
%     name=['0814_SimpleMC_MMDR_Qg_cornell_',num2str(iter)];
%     adj=A;
%     
%     load('NetScience-adj-1589.mat');
%     load netscience_groundtruth.txt  netscience_groundtruth -ASCII
%     groundtruth=netscience_groundtruth;   %记得修改   
%     adj=A;
%     
%     load('cora-adj-2708.mat');
%     load cora_groundtruth.txt  cora_groundtruth -ASCII
%     groundtruth=cora_groundtruth;   %记得修改
%     name=['0814_SimpleMC_MMDR_Qg_cora_',num2str(iter)];
%     adj=A;
%     
%     load wisconsin_adj.txt
%     adj=wisconsin_adj;
%     load wisconsin_groundtruth.txt  wisconsin_groundtruth -ASCII
%     groundtruth=wisconsin_groundtruth;   %记得修改
%     name=['0814_SimpleMC_MMDR_Qg_wisconsin_',num2str(iter)];
    
%     load ENZYMES_g163_adj.txt
%     adj=ENZYMES_g163_adj;
%     load ENZYMES_g163_groundtruth.txt  ENZYMES_g163_groundtruth -ASCII
%     groundtruth=ENZYMES_g163_groundtruth;   %记得修改
%     name=['0816_SimpleMC_MMDR_Qg_Enzymes_g163_',num2str(iter)];
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
    best_in_history_Qg=[];
    boundaryNode=culBoundaryNode(adj,groundtruth);
    Nodelist=(1:n);
    IntraNode=setdiff(Nodelist,boundaryNode);
    sigma=[1e-5,1e-3];
    canyujinhuaNode=ones(1,n);
    bubiancishu=zeros(1,n);
    % 1.3 pop construction
    pop = inital_Brain_pop(n, c, NP,groundtruth);
    fit_Qg = fitness_Qgcpp(pop,NP,adj,n,c);
    %2. the biased process used in the initialization
    bias_pop = bias_init_brain_pop( pop, c, n, NP, adj, groundtruth);
    bias_fit_Qg = fitness_Qgcpp(bias_pop,NP,adj,n,c);

    % 4  优秀的有偏个体保存到种群中
    win_index = find(bias_fit_Qg > fit_Qg); 
    win_number = length(win_index);
    % updata pop based on bias_pop
    pop(:,:,win_index) = bias_pop(:,:,win_index);
    fit_Qg(win_index) = bias_fit_Qg(win_index);
    %length(find(fit_Q>=bias_fit_Q))
    [bestfit,index]=max(fit_Qg); 
    bestx=pop(:,:,randsrc(1,1,find(fit_Qg==max(fit_Qg))'));

    pop_best_history=zeros(c,n,Gen); 
    % do the evolution iteratively
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

        Xinew=bound_brain_SOSFCD(Xinew,c,n,groundtruth);
        Xjnew=bound_brain_SOSFCD(Xjnew,c,n,groundtruth);
         canyu=find(canyujinhuaNode==1);
        differindex=canyu;
        differindex2=canyu;              
        Xinew_Qg=deltaQgcpp(Mutu_pop(:,:,i),Xinew,adj,c,differindex,Mutu_fit(i),m2,B);
        Xjnew_Qg = deltaQgcpp(Mutu_pop(:,:,j),Xjnew,adj,c,differindex2,Mutu_fit(j),m2,B);
            

        % (5)updata Mutu_pop and Mutu_fit
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
        Xinew=bound_brain_SOSFCD(Xinew,c,n,groundtruth);

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
        Para_vector = bound_brain_SOSFCD(Para_vector,c,n,groundtruth);

        canyu=find(canyujinhuaNode==1);
        differindex=canyu;
        Para_vector_Qg = deltaQgcpp(Para_pop(:,:,j),Para_vector,adj,c,differindex,Para_fit(j),m2,B);

        
        % (3)updata Para_pop and Para_fit
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
        clean_pop(:,:,r) = bound_brain_SOSFCD(clean_pop(:,:,r),c,n,groundtruth);

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
    Ubest_before=genbestx(:,:,g-1);
    canyu=find(canyujinhuaNode==1);
    for q=1:size(canyu,2)
        difference=std(genbestx(:,canyu(q),g));
        delta=(sigma(2)-sigma(1))*difference+sigma(1);
        if abs(Ubest_before(:,canyu(q))-genbestx(:,canyu(q),g))<delta
            bubiancishu(canyu(q))=bubiancishu(canyu(q))+1;
        elseif abs(Ubest_before(:,canyu(q))-genbestx(:,canyu(q),g))>=sigma(2)          
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
        for q=1:size(unchangeindex,2)
            i=unchangeindex(q);
            for r=1:NP
                pop(:,i,r)=bestx(:,i);
            end
            canyujinhuaNode(i)=0;
        end
        fit_Qg = fitness_Qgcpp(pop,NP,adj,n,c) ;     
    end
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
    
    end % while p<1
    fprintf('\n');
   


    Qg = max(best_in_history_Qg(g));
    Qglist(iter)=Qg;
    deltaTimeloop{iter}(:,:)=deltaQtime;
    best_in_history_Qgloop{iter}(:,:)=best_in_history_Qg;
    Genlist(iter)=g;

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
    mmdrnodes(iter)=n-sum(canyujinhuaNode);
    mmdrnodespercent(iter)=mmdrnodes(iter)/(n);
% save(['SOSFCD实验记录\',name],'Qglist','NMIlist','best_in_history_Qgloop','deltaTimeloop','mmdrnodespercent','mmdrnodes');        
end
Qgavg=mean(Qglist(:));
StdQg=std(Qglist); %无偏标准差 n-1
MaxQg=max(Qglist(:));
MinQg=min(Qglist(:));
NMIavg=mean(NMIlist(:));
