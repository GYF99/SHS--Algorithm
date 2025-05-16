%%%%% SimpleMC_Qg_LFR实验 
clc
clear all
%运行次数
looptime=3;
Qglist=zeros(1,looptime);
NMIlist=zeros(1,looptime);
best_in_history_Qgloop=cell(looptime,1);
deltaTimeloop=cell(looptime,1);
% 初始化结果和运行时间的数组
results = zeros(looptime, 1);
runtimes = zeros(looptime, 1);
for iter=1:looptime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%name
name=['LFR100_03_SimpleMC_0424_Qg_NP100_',num2str(iter)];
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % % 网络
% aaa=['5',num2str(iter-1)];
% name_net=['D:\研三上学期\实验代码-新\LFR网络\mu05\network',aaa,'-adj-1000.mat'];
% name_community=['D:\研三上学期\实验代码-新\LFR网络\mu05\community',aaa,'.dat'];
% load(name_net)
% community=load(name_community);
% adj=A;
% groundtruth=community(:,2); 
% load('LFR500_03-adj-500.mat')
% load LFR500_03_groundtruth.txt  LFR500_03_groundtruth -ASCII
% groundtruth=LFR500_03_groundtruth;   %记得修改
% adj=A;
load('LFR100_03-adj-100.mat')
load LFR100_03_groundtruth.txt  LFR100_03_groundtruth -ASCII
groundtruth=LFR100_03_groundtruth;   %记得修改
adj=A;
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
pop = inital_Brain_pop(n, c, NP,groundtruth);
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
% tic
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
clean_pop_before = crisp_fuzzypop(clean_pop, n, NP );  % NP个groundtruth？？
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
deltaTimeloop{iter}(:,:)=deltaQtime;
best_in_history_Qgloop{iter}(:,:)=best_in_history_Qg;
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
    
save(['SOSFCD实验记录-LFR\',name],'Qglist','NMIlist','best_in_history_Qgloop','deltaTimeloop');    
end
Qgavg=mean(Qglist(:));
StdQg=std(Qglist); %无偏标准差 n-1
MaxQg=max(Qglist(:));
MinQg=min(Qglist(:));
NMIavg=mean(NMIlist(:));



