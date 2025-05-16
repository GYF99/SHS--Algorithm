%%%%% SimpleMC_Qg实验 
clc
clear all
%运行次数
looptime=1;
Qglist=zeros(1,looptime);
NMIlist=zeros(1,looptime);
best_in_history_Qgloop=cell(looptime,1);
deltaTimeloop=cell(looptime,1);
% 初始化结果和运行时间的数组
results = zeros(looptime, 1);
runtimes = zeros(looptime, 1);
for iter=1:looptime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%name
name=['LFR5000_D50_dmax250_SimpleMC_0428_Qg_d10_',num2str(iter)];
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % % 网络
% 需要操作的地方，把以下三行中的每个mu05改成mu06/07/08 对应μ=0.5/0.6/0.7/0.8  
%     load('LFR网络-操作有效性验证\mu08_N1000\LFRnetwork-adj-1000.mat')
%     community=load('LFR网络-操作有效性验证\mu08_N1000\community.dat');
%     groundtruth=community(:,2);    
%     name=['0803_SimpleMC_Qg_mu05_',num2str(iter)];
%     adj=A;
load('LFR5000_d50_dmax250-adj-5000.mat')
load LFR5000_d50_dmax250_groundtruth.txt  LFR5000_d50_dmax250_groundtruth -ASCII
groundtruth=LFR5000_d50_dmax250_groundtruth;   %记得修改
adj=A;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%真实网络
% load('LesMis-adj-9.mat')
% load LesMis_groundtruth.txt  LesMis_groundtruth -ASCII
% groundtruth=LesMis_groundtruth;
% adj=A;
%     load('Karate-adj-34.mat')
% 
%     load Karate_groundtruth.txt  Karate_groundtruth -ASCII
%     groundtruth=Karate_groundtruth;   %记得修改
%     adj=A;

%     load('Yeast-adj-2361.mat');
%     load Yeast_groundtruth.txt  Yeast_groundtruth -ASCII
%     groundtruth=Yeast_groundtruth;   %记得修改    
%     adj=A;

%     load('brain47-adj-47.mat')
%     load brain47_groundtruth.txt  brain47_groundtruth -ASCII
%     groundtruth=brain47_groundtruth;   %记得修改
%     name=['0814_SimpleMC_Qg_brain47_',num2str(iter)];
%     adj=A;
%     
%     load('Dolphins-adj-62.mat')
%     load dolphins_groundtruth.txt  dolphins_groundtruth -ASCII
%     groundtruth=dolphins_groundtruth;   %记得修改
%     name=['0814_SimpleMC_Qg_dolphins_',num2str(iter)];
%     adj=A;

%     load('PolBooks-adj-105.mat')
%     load polbooks_groundtruth.txt  polbooks_groundtruth -ASCII
%     groundtruth=polbooks_groundtruth;   %记得修改
%     name=['0814_SimpleMC_Qg_Polbooks_',num2str(iter)];
%     adj=A;
%     load('macaque.mat')
%     load macaque1_groundtruth.txt  macaque1_groundtruth -ASCII
%     groundtruth=macaque1_groundtruth;   %记得修改
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
%     name=['0814_SimpleMC_Qg_cornell_',num2str(iter)];
%     adj=A;
%     
%     load('NetScience-adj-1589.mat');
%     load netscience_groundtruth.txt  netscience_groundtruth -ASCII
%     groundtruth=netscience_groundtruth;   %记得修改 
%     adj=A;

%     load('blogs-adj-1224.mat');
%     load blogs_groundtruth.txt  blogs_groundtruth -ASCII
%     groundtruth=blogs_groundtruth;   %记得修改 
%     adj=A;

%     
%     load('cora-adj-2708.mat');
%     load cora_groundtruth.txt  cora_groundtruth -ASCII
%     groundtruth=cora_groundtruth;   %记得修改
%     name=['0814_SimpleMC_Qg_cora_',num2str(iter)];
%     adj=A;
%     
%     load wisconsin_adj.txt
%     adj=wisconsin_adj;
%     load wisconsin_groundtruth.txt  wisconsin_groundtruth -ASCII
%     groundtruth=wisconsin_groundtruth;   %记得修改
%     name=['0814_SimpleMC_Qg_wisconsin_',num2str(iter)];

%     load ENZYMES_g163_adj.txt
%     adj=ENZYMES_g163_adj;
%     load ENZYMES_g163_groundtruth.txt  ENZYMES_g163_groundtruth -ASCII
%     groundtruth=ENZYMES_g163_groundtruth;   %记得修改
%     name=['0816_SimpleMC_Qg_Enzymes_g163_',num2str(iter)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % 1. Initialization
    g=1;
    p=0;
    edge=numedges(adj);
    n=numnodes(adj);
    c=max(groundtruth(:));
    threshold_value=0.25;
    Gen = 500;
    NP = 100;
    deltaQtime=[];

    % 1.3 pop construction
%     pop = inital_Brain_pop(n, c, NP,groundtruth);
    pop = inital_Brain_pop_1(n, c, NP,groundtruth);
    fit_Qg = fitness_Qgcpp(pop,NP,adj,n,c);

    %2. the biased process used in the initialization
    bias_pop = bias_init_brain_pop( pop, c, n, NP, adj, groundtruth);
    bias_fit_Qg = fitness_Qgcpp(bias_pop,NP,adj,n,c);
    
    win_index = find(bias_fit_Qg > fit_Qg); 
    win_number = length(win_index);
    pop(:,:,win_index) = bias_pop(:,:,win_index);
    fit_Qg(win_index) = bias_fit_Qg(win_index);

    % 3. Main loop
    pop_best_history=zeros(c,n,Gen); 
    
    % do the evolution iteratively
    while p < 1
%     tic
    % 3.1 Mutualism
    Mutu_pop=pop;
    Mutu_fit=fit_Qg;
    better_number=0;
    for i=1:NP
        % (1)Xbest
        [bestfit,index]=max(Mutu_fit);   
        bestx=Mutu_pop(:,:,randsrc(1,1,find(Mutu_fit==max(Mutu_fit))'));
        seed=setdiff(randperm(NP),i); 
        j=randsrc(1,1,seed);
        mutual_vector=0.5.*(Mutu_pop(:,:,i)+Mutu_pop(:,:,j));
        BF1=round(1+rand); 
        BF2=round(1+rand);
        % (4)Xinew & Xjnew
        Xinew=Mutu_pop(:,:,i)+rand(c,n).*(bestx-BF1.*mutual_vector); 
        Xjnew=Mutu_pop(:,:,j)+rand(c,n).*(bestx-BF2.*mutual_vector);
        Xinew=bound_brain_SOSFCD(Xinew,c,n,groundtruth);
        Xjnew=bound_brain_SOSFCD(Xjnew,c,n,groundtruth);
        Xinew_Qg = fitness_Qgcpp(Xinew,1,adj,n,c);
        Xjnew_Qg = fitness_Qgcpp(Xjnew,1,adj,n,c);
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

    % 3.2 Commensalism
    Comm_pop=Mutu_pop;
    Comm_fit=Mutu_fit;
    better_number=0;
    for i=1:NP
        % (1)Xbest
        [bestfit,index]=max(Comm_fit); 
        bestx=Comm_pop(:,:,randsrc(1,1,find(Comm_fit==max(Comm_fit))'));
        % (2)Xj~=Xi
        seed=setdiff(randperm(NP),i); 
        j=randsrc(1,1,seed);
        % (3)Xinew 
        Xinew=Comm_pop(:,:,i)+(rand(c,n)*2-1).*(bestx-Comm_pop(:,:,j));
        Xinew=bound_brain_SOSFCD(Xinew,c,n,groundtruth);
        Xinew_Qg = fitness_Qgcpp(Xinew,1,adj,n,c);
        if Xinew_Qg > Comm_fit(i)
           Comm_pop(:,:,i)=Xinew;
           Comm_fit(i)=Xinew_Qg;
           better_number=better_number+1;
        end
    end 

    % 3.5 Parasitism  / better_number = 4
    Para_pop=Comm_pop;
    Para_fit=Comm_fit;
    better_number=0;
    for i=1:NP
        % (1)Xj~=Xi
        seed=setdiff(randperm(NP),i); 
        j=randsrc(1,1,seed);   
        % (2)Parasite Vector
        Para_vector=Para_pop(:,:,i);
        seed=randperm(n);
        pick=seed(1:ceil(rand*n));   
        Para_vector(:,pick) = inital_Brain_pop(length(pick), c, 1,groundtruth);
        Para_vector = bound_brain_SOSFCD(Para_vector,c,n,groundtruth);
        Para_vector_Qg = fitness_Qgcpp(Para_vector,1,adj,n,c);
        %FE=FE+1;        
        % (3)updata Para_pop and Para_fit
        if Para_vector_Qg > Para_fit(j)
           Para_pop(:,:,j)=Para_vector;
           Para_fit(j)=Para_vector_Qg;
           better_number=better_number+1;
        end
    end 


    % 3.6 Parasitism clean_uper  
    clean_pop = Para_pop;
    clean_fit = Para_fit;
    clean_pop_before = crisp_fuzzypop(clean_pop, n, NP );  
    clean_pop_after = clean_up_random( clean_pop_before, n, NP, adj, threshold_value );

    for r=1:NP
        for node = 1:n
            if clean_pop_after(r,node) ~= clean_pop_before(r,node)
                k = clean_pop_after(r,node);
                clean_pop(k,node,r) = clean_pop(k,node,r) + 0.5;
            end
        end
        clean_pop(:,:,r) = bound_brain_SOSFCD(clean_pop(:,:,r),c,n,groundtruth);
        clean_fit(r) =  fitness_Qgcpp(clean_pop(:,:,r),1,adj,n,c);
    end


    for r =1:NP
        if clean_fit(r) > Para_fit(r)
           Para_pop(:,:,r) = clean_pop(:,:,r);
           Para_fit(r) = clean_fit(r);
        end
    end

    % 3.7 updata pop and fit
    pop = Para_pop;
    fit_Qg = Para_fit;
%     fit_Qg = fitness_Qgcpp(pop,NP,adj,n,c);
    % 3.8 recode bestfit and bestx
    [best_Qg,index]=max(fit_Qg); 
    bestx=pop(:,:,index);
    pop_best_history(:,:,g)=bestx;
    best_in_history_Qg(g)=best_Qg; 
    
    if g>1
        deltaQtime(g)=deltaQtime(g-1)+toc;
    else
        deltaQtime(g)=toc;
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
    results(iter) = Qg;
    runtimes(iter) = toc;  %结果计时并存储运行时间
    deltaTimeloop{iter}(:,:)=deltaQtime;
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
    

    
save(['SOSFCD实验记录\',name],'Qglist','NMIlist','best_in_history_Qgloop','deltaTimeloop');    
end
Qgavg=mean(Qglist(:));
StdQg=std(Qglist); %无偏标准差 n-1
MaxQg=max(Qglist(:));
MinQg=min(Qglist(:));
NMIavg=mean(NMIlist(:));



