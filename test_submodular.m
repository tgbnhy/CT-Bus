%%%% test the submodular by randomly choosing various number of edges.
test_times = 50;
test_maximum_edge = 50;
percent = zeros(test_maximum_edge-1,test_times);
for number_edge = 2:test_maximum_edge
    counterexp = 0;
    for i=1:test_times
        %generate k edges and 
        A_temp = A;%construct the new edge by check the path
        r = randi([1 95304],1,number_edge); %160790 for nyc
        sum = 0;
        for j=1:number_edge
            A_temp(B(r(j),1),B(r(j),2)) = 1;
            A_temp(B(r(j),2),B(r(j),1)) = 1;
            sum = sum + B(r(j),7);
        end
        
        %{
        ix = idx(i);
        iy = idx(i+1);
        A_temp(B(ix,1),B(ix,2)) = 1;%update for a new edge
        A_temp(B(ix,2),B(ix,1)) = 1;
        A_temp(B(iy,1),B(iy,2)) = 1;%update for a new edge
        A_temp(B(iy,2),B(iy,1)) = 1;
        if SortB(i) == SortB(i+1) % two equal scores
            continue;
        end
        sum = SortB(i) + SortB(i+1);
        %}
        K1 =  K1_origin.K1;%
        gap = natural_connectivity(A_temp, graph_dimension, K1, reps, iter) - base;
    %    disp("a"+sum);
    %    disp("b"+gap);
        if sum < gap
           counterexp = counterexp + 1;% get the number of counter example
        end
        percent(number_edge-1,i) = (gap - sum)/sum;
        if gap == 0
            percent(number_edge-1,i) = 0;
        end
        disp(number_edge+","+(gap - sum)/sum);
    end
 %   disp(counterexp)
end
%disp(percent);

