%%%%%% parameters

weight = 0.5; % the weight to balance each dimension, when we set it as 1, it is a baseline using orienteering
k = 30; % number of new edges in a path
number_of_turns = 1;% number of turns, candidates
seeding_number = 5000;% number of candidate edges as initial paths, 174579 for nyc, 101774, and inserting all neighbors

max_len = k; % maximum number of edges in a path
max_iter = 100000;
area_lable = 5;% 1: Brooklyn, 2: Manhattan, 3: Statan, 4: Queens, 5 Bronx
area_indicator = 1;% only for NYC dataset

city_choice = 1; % 0:1;% test two cities, nyc (1) or chicago (0)
computation_choice = 0;%full bound computation or fast easy bound, fast (0) or slow (1)

objectives = [];%for logs

for test_nyc = city_choice
    optimization_choice = 0;%0:2; %0 for full optimization, 1 for all neighbors, 2 for not use domination table
    %% import the graph and dataset.
    stoplocations = load('/Users/sw160/Desktop/nyc/nyc_transit_stop_locations.mat');
    neighor_stops = load('/Users/sw160/Desktop/nyc/new_nyc_transit_neighbors0.5_0.2_pre.mat');
    K1_origin = load('/Users/sw160/Desktop/nyc/new_nyc_random_K1.mat');%for conn, randomly in the begining, use a fixed one
    C = load('/Users/sw160/Desktop/nyc/nyc_transit_weighted.mat');
    nodes_Number = 12340;

    if ~test_nyc
        area_indicator = 0;
        stoplocations = load('/Users/sw160/Desktop/transit_network/chicago/chi_transit_stop_locations.mat');
        neighor_stops = load('/Users/sw160/Desktop/transit_network/chicago/chi_transit_neighbors0.5_0.2.mat');
        K1_origin = load('/Users/sw160/Desktop/transit_network/chicago/new_chi_random_K1.mat');%for conn, randomly in the begining, use a fixed one
        C = load('/Users/sw160/Desktop/transit_network/chicago/chi_transit_weighted.mat');% this is the existing edges
    end

    stoplocations = stoplocations.arr;%it stores the location information of stops
    %neighor_stops = load('/Users/sw160/Desktop/nyc/nyc_transit_neighbors0.5.mat');
    B = neighor_stops.F;%
    %{
    if ~test_nyc
        
    else
        B = neighor_stops.B;%
    end
    %}
    %clean B, remove all the duplicates
    %clean_B = unique(B(:,1),'rows','stable');
    neighbors_node = containers.Map;%store it in the list
    neighbor_number = containers.Map;%store it in the list
    
    A = 0;
    cursor = 0;
    ranked_weight = [];%
    for i = 1:size(B,1)%store the neighbors' index
        if B(i)~=A
            neighbors_node(int2str(B(i))) = i;
            if i>1
                neighbor_number(int2str(B(i-1))) = i - cursor - 1;
            end
            cursor = i;
        end
        A = B(i);%the first item,
    end
    neighbor_number(int2str(B(i))) = i - cursor;%

    %% construct the transit graph
    C = C.arr;
    G = graph();
    tic % load the graph and insert to query
    for i = 1:size(C,1)
        if or(C(i,1)>nodes_Number, C(i,2)>nodes_Number)%for new york dataset
            continue;
        end
        G = addedge(G, C(i,1), C(i,2), i); % the weight is the index in D
        C(i,7) = 0;
        C(i,8) = 0;
        C(i,9) = 0;
    end
    
    toc
    A = adjacency(G);%adjacent matrix
    graph_dimension = size(A,1);
    %[SortC,idxC] = sort(C(:,6),'descend');%sort existing edges by the edge frequency

    if weight == 1
        D = B; % only consider the new edges
        G = graph();% empty the graph
    else
        D = [B;C];%combine possible edges and new edges, 
    end
    % tranfer the weight of edges to vertex.
    
    %% parameters and data structures for connectivity
    reps = 50;% number of random vectors
    iter = 10;% number of krylov iterations

    K1 =  K1_origin.K1;%
    disp(size(K1))
    disp(size(A))
    tic
    base = natural_connectivity(A, graph_dimension, K1, reps, iter);
    toc
    %% insert the heaviest unlinked edges into the queue 
    [SortA,idxA] = sort(D(:,6),'descend');%sort by the edge frequency
    bound_weight = SortA(1:max_len);%storing the highest max_len
%      disp(size(bound_weight));
    lb_weight_max = 0;
    for i=1:max_len
        lb_weight_max = lb_weight_max + SortA(i);% choose the maximum k edges' sum for normalization
    end

    %% upper bound with k edges
    [SortB,idx] = sort(B(:,7),'descend');%sort the edges by their connectivity increment
    bound_conn = SortB(1:k);%storing the highest max_len
    bottom_conn = SortB(k);%for bound estimation
    connectivity_ub = 0;%
    
    for i=1:k
        connectivity_ub = connectivity_ub + SortB(i);% choose the maximum k edges' sum to compute
    end
    lb_conn_max = connectivity_ub;
    
    tic
    musco_bound = path_upper_bound(k, A, graph_dimension, base) - base;
    toc
    %{
    tic
    erada_connectivity_ub = log((1+(exp(sqrt(size(C,1)+k+size(C,1)+k))-1)/graph_dimension))-base;% from Erada index
    toc

    tic
    musco_general_bound = general_upper_bound(k, A, graph_dimension, base) - base;
    toc
    
    disp("bound: "+erada_connectivity_ub+","+musco_general_bound+","+musco_bound+","+connectivity_ub);
    %}
    %% integrated equity increment
    for i=1:size(D)
        D(i,8) = weight*D(i,6)/lb_weight_max + (1-weight)*D(i,7)/lb_conn_max;
    end
    [SortA,idxA] = sort(D(:,8),'descend');%sort by the edge frequency
    bound_weight = SortA(1:max_len);
    bottom_weight = SortA(max_len);
    
    %{
    if area_indicator % revise the bound for small area, as the overall bound is too big,
        [lb_weight_max, connectivity_ub, DD] = combine_tables(B, C, stoplocations, area_lable, k);
        for i=1:size(DD)% integrated equity increment
          DD(i,8) = weight*DD(i,6)/lb_weight_max + (1-weight)*DD(i,7)/lb_conn_max;
        end
        SortA = sort(DD(:,8),'descend');%sort by the edge frequency
        bound_weight = SortA(1:max_len);
        bottom_weight = SortA(max_len);
    end
    disp(connectivity_ub);
    disp(lb_weight_max)
    %}
    %%running time estimation
    estimted_running_time = 0;
    time_cost_conn_compute = 0.09;
    time_cost_ub_compute = 0.13;
    if ~test_nyc
        time_cost_conn_compute = 0.05;
        time_cost_ub_compute = 0.08;
    end
    %[SortE,idxE] = sort(D(:,6),'descend');%sort by the edge frequency
    %[SortF,idxF] = sort(D(:,7),'descend');%sort by the edge frequency
    %disp(SortE(1:1000)), %present the top 50 edges
    %disp(SortF(1:1000))
    
    for full_computation = computation_choice
        if full_computation==1
            optimization_choice = 0;%no need to try other method.
        end
        %% the main algorithm part.
        for optimization_indicator = optimization_choice%test three optimization
            best_path = [];% store the current results
            min_equity = 0;% maximum value
            n = 1;
            q = PriorityQueue(1);% a priority queue to rank the candiates by the first column
            lb = 1;%the initial bound
            checked_edges = containers.Map;% for domination on candidates have same ends
            seeding_edges = containers.Map;
            counter = 0;
            for i=1:size(B,1)
                ix = idxA(i);
                newpath = [D(ix,1) D(ix,2)];
                
                if area_indicator % plan for an area
                    if or(area_lable ~= stoplocations(D(ix,1),4), area_lable~= stoplocations(D(ix,2),4))
                        continue;
                    end
                end
                
                if or(D(ix,1)>nodes_Number, D(ix,2)>nodes_Number) %for new york dataset
                    continue;
                end
                
                if or(isKey(seeding_edges, D(ix,1)+","+D(ix,2)), isKey(seeding_edges, D(ix,2)+","+D(ix,1)))%avoid repetitive
                    continue
                end
                
                seeding_edges(D(ix,1)+","+D(ix,2)) = 1;
                seeding_edges(D(ix,2)+","+D(ix,1)) = 1;
               % disp(newpath);
                %some edges are repetitive.
                bw = max_len;
                if D(ix,8) < bottom_weight
                    bw = bw - 1;%
                    lb = 1 - (bottom_weight-D(ix,8));
                end
                [frequency, conn, new_equity] = equity(A, newpath, graph_dimension, weight, connectivity_ub, lb_weight_max, K1, reps, iter, 0, D(ix,6),0, D(ix,7),0,base);
                
                if new_equity > min_equity
                    min_equity = new_equity;
                    best_path = newpath;
                    best_array = [lb*-1 bw 0 frequency conn D(ix,1) D(ix,2)];
                %    best_array(4) = newfre;
                %    best_array(5) = newconn;
                end
                
                                
                q.insert([lb*-1 bw 0 frequency conn D(ix,1) D(ix,2)]);%how to indicate new edges for later update, and the last two edges
                counter = counter+1;
                if counter>seeding_number
                    break
                end
            end

            bad = 0;
            tic
            disp("hello")
            while n < max_iter && q.size()>0 % the queue is not empty
                candidate = q.remove();
                lower_bound = candidate(1)*-1;% the potential
                
                if min_equity > lower_bound% current result already smaller than the rest potential
                    break
                end
                
                size_path = size(candidate);
                bw = candidate(2); % cursor for upper bound
                bc = candidate(3); % record number of turns
                frequency = candidate(4);% current frequency (demand)
                conn = candidate(5);% current conn
                path = candidate(6:size_path(2));%the path, new bound
                disp(n+" "+min_equity+" " + lower_bound+ " " +length(best_path)+" "+q.size()); %the lower bound 
                if mod(n,100) == 0
                    cc = n/100;
                    if test_nyc==0 && full_computation == 1
                        objectives(cc,4) = min_equity;
                    elseif test_nyc==1 && full_computation == 1
                        objectives(cc,8) = min_equity;
                    else
                        objectives(cc,optimization_indicator+1+test_nyc*4) = min_equity;
                    end
                    disp(min_equity);
                end
                last_node = candidate(size_path(2));
                first_node = candidate(6);
                ne = [];
                if weight < 1
                    N = transpose(neighbors(G,last_node));%get existing neighbor in the graph
                    for i=1:length(N)
                        ne = [ne, size(B) + G.Edges.Weight(findedge(G,last_node,N(i)))];%locate the edge
                    end
                end
                TF = isKey(neighbors_node, int2str(last_node));
                if TF==1%access the linked possible, 
                    for i = 0: neighbor_number(int2str(last_node))%get possible edge
                        ne = [ne, neighbors_node(int2str(last_node))+i];%add new the neighbor lists
                    end
                end
             %   disp(ne)% show 
                max_increment = 0;
                best_xx = 0;
                array = [];
                bcc = bc;
                best_neighbor = 0;
                for i = 1: length(ne)% check all the neighbors, choose the one with maximum increment
                    index_neighbor = ne(i);
                    if index_neighbor >= size(D,1)
                        continue;
                    end
                    xx = D(index_neighbor, 2);
                    
                    if D(index_neighbor, 2) == last_node %avoid the edge from graph
                        xx = D(index_neighbor, 1);
                    end
                    if area_indicator % plan for an area
                        if area_lable ~= stoplocations(xx,4)
                            continue;
                        end
                    end

                    if ismember(xx,path)%circle-free, not sure it i
                        continue
                    end

                    new_bc = bc; 
                    last_s_node = candidate(size_path(2)-1);
                    anlges = abs(angle(stoplocations(last_node, 3), stoplocations(last_node, 2), stoplocations(xx, 3), stoplocations(xx, 2), stoplocations(last_s_node, 3), stoplocations(last_s_node, 2))*180/pi);

                    if anlges < 90 % big turn, not allow, as it will go back
                        continue;
                    end
                    if anlges < 135 % small turn, can have some
                        new_bc = bc+1;
                    end
                    if new_bc > number_of_turns
                        continue
                    end
 
                    % optional
                    newpath = [path, xx];% add to the last end
                    K1 =  K1_origin.K1;%
                    [newfre, newconn, new_equity] = equity(A, newpath, graph_dimension, weight, connectivity_ub, lb_weight_max, K1, reps, iter, frequency, D(index_neighbor,6),conn, D(index_neighbor, 7), full_computation,base);
              %      disp(new_equity)
                    % domination table
                    %{
                    if xx < first_node
                        checking = xx+"_" + first_node; %% adding to checking lists
                    else
                        checking = first_node + "_"+ xx; %% adding to checking lists
                    end
                    if isKey(checked_edges, checking)==1
                        if checked_edges(checking) > new_equity
                            % a counter for pruned route
                            continue;
                        else
                            checked_edges(checking) = new_equity;
                        end
                    else
                        checked_edges(checking) = new_equity;
                    end
                    %}
                    estimted_running_time = estimted_running_time + time_cost_conn_compute;
                    if new_equity > min_equity
                        min_equity = new_equity;
                        best_path = newpath;
                        best_array = candidate;
                        best_array(4) = newfre;
                        best_array(5) = newconn;
                    end

                    %choose the one with maximum increment on equity
                    incrementcc = D(index_neighbor,8);
                    if full_computation
                        incrementcc = new_equity;
                    end

                    if incrementcc > max_increment
                        max_increment = incrementcc;
                        best_neighbor = index_neighbor;
                        best_xx = xx;
                        bcc = new_bc;
                    %   array = [lb*-1 new_bw bc start_angle D(index_neighbor,9) newfre newconn, newpath];
                    end
                end

                if max_increment>0
                    path = [path, best_xx];
                    bottom_weight = bound_weight(bw);
                    if D(best_neighbor,8) < bottom_weight
                        bw = bw - 1;%
                        lower_bound = lower_bound - (bottom_weight - D(best_neighbor,8));
                    end
                    frequency = frequency + D(best_neighbor,6)/lb_weight_max;
                    conn = conn + D(best_neighbor,7)/lb_conn_max;
                    bc = bcc;
                else
                    bad = bad+1;
                end

                if bw<1 % no more edges to be added, as it aready has k edges
                    continue
                end
                ne = [];
                if weight < 1
                    N = transpose(neighbors(G,first_node));%get existing neighbors in the graph
                    for i=1:length(N)
                        ne = [ne, size(B)+G.Edges.Weight(findedge(G,first_node,N(i)))];%locate the edge in D
                    end
                end
                TF = isKey(neighbors_node, int2str(first_node));
                if TF==1%access the linked possible, 
                    for i = 0: neighbor_number(int2str(first_node))%get possible edge
                        ne = [ne, neighbors_node(int2str(first_node))+i];%add new the neighbor lists
                    end
                end

                max_increment = 0;
                array = [];
                for i = 1: length(ne)
                    index_neighbor = ne(i);
                    if index_neighbor >= size(D,1)
                        continue;
                    end
                    xx = D(index_neighbor, 2);
                    if D(index_neighbor, 2) == last_node %avoid the edge from graph
                        xx = D(index_neighbor, 1);
                    end
                    if area_indicator % plan for an area
                        if area_lable ~= stoplocations(xx,4)
                            continue;
                        end
                    end
                    if ismember(xx,path)%circle-free
                        continue
                    end
                    new_bc = bc;
                    last_s_node = candidate(7);
                    anlges = abs(angle(stoplocations(first_node, 3), stoplocations(first_node, 2), stoplocations(xx, 3), stoplocations(xx, 2), stoplocations(last_s_node, 3), stoplocations(last_s_node, 2))*180/pi);
                    if anlges < 90 % not going back, pruned 
                        continue;
                    end
                    if anlges < 135 % one more turn in the route,
                        new_bc = bc+1;
                    end
                    if new_bc > number_of_turns
                        continue
                    end
                    

                    newpath = [xx, path];%add to the start position
                    K1 =  K1_origin.K1;%
                    [newfre, newconn, new_equity] = equity(A, newpath, graph_dimension, weight, connectivity_ub, lb_weight_max, K1, reps, iter, frequency, D(index_neighbor,6), conn, D(index_neighbor, 7),full_computation,base); 
                    estimted_running_time = estimted_running_time + time_cost_conn_compute;
                    if optimization_indicator ~= 2 % use domination table 
                        if xx < last_node
                            checking = xx +"_"+last_node; %% adding to checking lists
                        else
                            checking = last_node + "_"+ xx; %% adding to checking lists
                        end
                        if isKey(checked_edges, checking)==1
                            if checked_edges(checking) > new_equity
                                continue;
                            else
                                checked_edges(checking) = new_equity;
                            end
                        else
                            checked_edges(checking) = new_equity;
                        end
                    end

                    if new_equity > min_equity
                        min_equity = new_equity;
                        best_path = newpath;
                        best_array = candidate;
                        best_array(4) = newfre;
                        best_array(5) = newconn;
                    end

                    %choose the one with maximum increment on equity
                    incrementcc = D(index_neighbor,8);
                    if full_computation
                        incrementcc = new_equity;
                    end

                    if or(incrementcc > max_increment, optimization_indicator) % only use maximum, or use every one   
                        max_increment = incrementcc;
                        lb = lower_bound;
                        bottom_weight = bound_weight(bw);
                        if D(index_neighbor,8) < bottom_weight
                            new_bw = bw - 1;%
                            lb = lb - (bottom_weight-D(index_neighbor,8));
                        end
                        if new_bw >= 1 
                           array = [lb*-1 new_bw new_bc newfre newconn newpath];
                        end
                        if optimization_indicator == 1
                            if length(array) < (5 + max_len) && length(array) >= 7
                                q.insert(array); % insert every neighbor
                            end
                        end
                    end
                end

                if max_increment > 0 && optimization_indicator ~= 1
                    if length(array) < (5 + max_len) && length(array) >= 7
                        if -array(1) > min_equity
                            q.insert(array);% only insert the best 
                        end
                        estimted_running_time = estimted_running_time + time_cost_ub_compute;
                    end
                else
                    bad = bad+1;
                    a = 0;% no more edges into the queue, just skip this process
                end
                n = n+1;
            end
            disp("#bad:" + bad);
            disp("we found a path: ");
            disp(best_path);
            count_new_edges = 0;

            fileIDnew = fopen('./new_edge_b.txt','w');
            A_temp = A;
            content = 'geometry\n"{""type"": ""LineString"", ""coordinates"": [';
            teb = 0;
            for i=1:size(best_path,2)-1%count the number of new edges, we can also recompute the connectivity
                if A(best_path(i),best_path(i+1)) == 0
                    id = best_path(i);
                    count_new_edges = count_new_edges + 1;
                    fprintf(fileIDnew,'[%3.6f,%3.6f],',stoplocations(id,3), stoplocations(id,2));
                end
                A_temp(best_path(i),best_path(i+1)) = 1;
                A_temp(best_path(i+1),best_path(i)) = 1;
                for adsad = 1:size(D,1)
                    if best_path(i) == D(adsad, 1) && best_path(i+1) == D(adsad, 2)
                        teb = teb + D(adsad, 7)/connectivity_ub;
                    end
                end
            end
        %    fprintf(fileIDnew,'[%3.6f,%3.6f]]}"\n',stoplocations(length(best_path),3), stoplocations(length(best_path),2));

            disp("increased equity: "+min_equity);
            disp("#new edges: "+count_new_edges);
            base1 = natural_connectivity(A_temp, graph_dimension, K1, reps, iter);
            disp("real increase: "+(base1-base)/connectivity_ub);
            disp("estimated increase conn: "+best_array(5));
            disp("recompute: "+teb);
            disp("estimated increase demand: "+best_array(4));
            disp("estimated time on conn and ub: "+ estimted_running_time);
            %write file for mapv visualization
            fileID = fopen('/Users/sw160/Documents/mapv-2.0.12/travis/data/fairbus_k501.txt','w');
            fileID3 = fopen('/Users/sw160/Desktop/nyc/route_stops.txt','w');
            fileID2 = fopen('/Users/sw160/Desktop/nyc/busroutes_new.txt','w');

            if ~test_nyc
                fileID = fopen('/Users/sw160/Documents/mapv-2.0.12/travis/data/chi_fairbus_k501.txt','w');
                fileID3 = fopen('/Users/sw160/Desktop/transit_network/chicago/route_stops.txt','w');
                fileID2 = fopen('/Users/sw160/Desktop/transit_network/chicago/busroutes_new.txt','w');
            end


            fprintf(fileID,'geometry\n"{""type"": ""LineString"", ""coordinates"": [');
            for i=1:length(best_path)
                id = best_path(i);
                fprintf(fileID2,'%d\n',id);
                fprintf(fileID3,'%3.6f,%3.6f\n',stoplocations(id,3), stoplocations(id,2));%for finding related routes in the edge
                fprintf(fileID,'[%3.6f,%3.6f]',stoplocations(id,3), stoplocations(id,2));
                if i~=length(best_path)
                    fprintf(fileID,',');
                end
            end
            %write the stops location, and new edges,
            fprintf(fileID,']}"\n');% \n is essential, it will not work if without it.
            fclose(fileID);
            fclose(fileID2);
            fclose(fileID3);
            q.clear();
        end
    end
end
toc

logs = "/Users/sw160/Desktop/nyc/exp_logs";%give every d
%{
if test_nyc
    if area_indicator
        if area_lable == 1
            logs = logs+"_brooklyn";
        elseif area_lable == 2
            logs = logs+"_manhattan";
        elseif area_lable == 3
            logs = logs+"_statan";
        elseif area_lable == 4
            logs = logs+"_queens";
        elseif area_lable == 5
            logs = logs+"_bronx";
        end
    else
        logs = logs+"_nyc";
    end
else
    logs = logs+"_chicago";
end
%}
logs = logs+"_"+k+"_"+weight+"_"+number_of_turns+"_"+seeding_number;
logs = logs+".xls";
delete(logs);
writematrix(objectives, logs);