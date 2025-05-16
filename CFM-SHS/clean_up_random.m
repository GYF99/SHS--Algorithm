
function pop = clean_up_random( pop, n, pop_size, adj, threshold_value )

for k = 1:pop_size
    
    all_node_index = 1:n;
    get_num = floor(rand*n)+1;
    use_node_index = [];
    
    for cu_i = 1:get_num
        cur_rand_index = randi([1,length(all_node_index)],1,1);     
        use_node_index = [use_node_index, all_node_index(cur_rand_index)];
        all_node_index(cur_rand_index) = [];      
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if length(find(21==use_node_index))==0
%        use_node_index = [use_node_index, 21]; 
%        get_num=get_num+1;
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     all_node_index = [];

    for rand_i = 1:get_num
        
        i = use_node_index(rand_i);  
        all_adj_node = [i];        
        for j = 1:n            
            if j~=i && adj(i,j)==1               
                all_adj_node = [ all_adj_node, j ];                
            end            
        end
        
        all_adj_comm = pop( k, all_adj_node ); 
       
        left_adj_comm = all_adj_comm( 2 : length(all_adj_comm) );
        
        comm_id_i = all_adj_comm(1);
        
        different_comm_id = find( left_adj_comm ~= comm_id_i );
        
        degree_i = sum( adj(i,:) );
        % CV
        cv_i = length(different_comm_id) / degree_i;
        
        if cv_i >= threshold_value
            
            temp_comm = left_adj_comm;
            
            comm_list=unique(temp_comm);  
            comm_num=[];
            
            while ~isempty(comm_list)                
                cur_comm = comm_list(1);                
                all_node = find( temp_comm==cur_comm );                
                cur_num = length( all_node );
                comm_num=[comm_num cur_num];
                         
                temp_comm( all_node ) = [];    
                comm_list(1)=[];
            end
            
            comm_list=unique(left_adj_comm);
         
            max_comm_num=max(comm_num);
            max_comm_list_index=find(comm_num==max_comm_num);
            max_comm_id=[];
            for i_index=1:length(max_comm_list_index)
                max_comm_id=[max_comm_id comm_list(max_comm_list_index(i_index))];
            end          
                     
            comm_num(max_comm_list_index)=[];
            comm_list(max_comm_list_index)=[];           
            
            if length(comm_list)>0
               second_max_comm_num=max(comm_num);
               second_max_comm_list_index=find(comm_num==second_max_comm_num);
               second_max_comm_id=[];
               for i_index=1:length(second_max_comm_list_index)
                   second_max_comm_id=[second_max_comm_id comm_list(second_max_comm_list_index(i_index))];
               end
            else
               second_max_comm_num=0;
               second_max_comm_id=[];
            end
            
            comm_id_for_choice=[max_comm_id second_max_comm_id];
            
            total_comm_num=length(max_comm_id).*max_comm_num+length(second_max_comm_id).*second_max_comm_num;
            for i_index=1:length(comm_id_for_choice)
                comm_id=comm_id_for_choice(i_index);
                if i_index==1 && length(find(comm_id==max_comm_id))>0
                    
                    pop( k, i ) = comm_id;
                elseif i_index>1 && length(find(comm_id==max_comm_id))>0
                    
                   gailv=0.5;
                   if rand<gailv
                      
                      pop( k, i ) = comm_id;
                   end
                elseif i_index>1 && length(find(comm_id==second_max_comm_id))>0            
                     
                     if rand<0.2
                        
                        pop( k, i ) = comm_id;
                     end 
                end
            end
   
        end  
        
    end  
    
end 

%end