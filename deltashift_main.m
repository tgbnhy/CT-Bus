%arr = Problem.A;
C = load('./nyc_data/nyc_transit_weighted.mat');
C = C.arr;
G = graph();
nodes_Number = 12340;
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

arr = adjacency(G);

n = size(arr,1); %Dim of graph
m = 50;         %Seq length
l2 = 50;         %Number of mat-vec for deltas (A2-A1)
iter = 10;       %Iteration for Lanczos algorithm
exp = 1;        %Number of experiments
%l_true = 500;   %No. of RV to approximate true trace

l1 = 2*l2;       %Mat-vec for hutchinson
l3 = l1;       %Mat-vec for first matrix in the sequence

trace_truth = zeros([exp,m]);
trace_hutch = zeros([exp,m]);
trace_ds_hutch = zeros([exp,m]);
trace_ds_combined = zeros([exp,m]);


for e = 1:exp
    disp(e);
    
    %Approximation to true trace
    A1 = full(arr);
    
    
    %G_T = 2*randi(2,n,l_true)-3;
    %A1G_T = lanczos(A1, G_T, @exp, iter);
    %trace_truth(e,1) = trace(G_T' * A1G_T)/l_true;
    trace_truth(e,1) = trace(expm(A1));
    
    
    
    %Trace A1 using hutchinsons
    G = 2*randi(2,n,l3)-3;
    A1G = lanczos(A1, G, @exp, iter);
    trace_hutch(e,1) = trace(G' * A1G)/l3;          %Hutch
    
    trace_ds_hutch(e,1) = trace_hutch(e,1);         %DS(Hutch)
    trace_ds_combined(e,1) = trace_hutch(e,1);      %Hutch + DS(Hutch++)
    prev_var_hutch = (2/l3^2) * trace(A1G' * A1G);  %Var for DS(Hutch)
    prev_var_combined = prev_var_hutch;             %Var for DS(Hutch++)
    
    
    
    for z = 2:m
        disp(z)
        %Choose two nodes at random and add edge
        delta = zeros(n,n);
        no_nodes = 2;
        nodes = randperm(n, no_nodes);
        edges = nchoosek(nodes,2);
        for i = 1:size(edges,1)
            delta(edges(i,1), edges(i,2)) = 1;
            delta(edges(i,2), edges(i,1)) = 1;
        end
        
        A2 = A1 + delta;
        
        %True value
        %G_T = 2*randi(2,n,l_true)-3;
        %A2G_T = lanczos(A2, G_T, @exp, iter);
        %A2G = A2_exp * G;
        %trace_truth(e,z) = trace(G_T' * A2G_T)/l_true;
      %  trace_truth(e,z) = trace(expm(A2));
        
        %Hutch
     %   G = 2*randi(2,n,l1)-3;
    %    A2G = lanczos(A2, G, @exp, iter);
    %    trace_hutch(e,z) = trace(G' * A2G)/l1;
        
        %DeltaShift - Hutchinson
    %    [trace_ds_hutch(e,z), prev_var_hutch]= deltashift_hutch(A1, A2, trace_ds_hutch(e,z-1), prev_var_hutch, l2, iter);
        
       
        %DeltaShift - Hutch++
        tic
        [trace_ds_combined(e,z), prev_var_combined] = deltashift_hutchpp(A1, A2, trace_ds_combined(e,z-1), prev_var_combined, l2, iter);
        toc
        disp(trace_ds_combined(e,z))
        A1 = A2;
    end

end


x = 1:m;

plot(x,mean(abs(trace_truth-trace_hutch)./trace_truth,1),x,mean(abs(trace_truth-trace_ds_hutch)./trace_truth,1),x,mean(abs(trace_truth-trace_ds_combined)./trace_truth,1));
legend("Hutchinson","DS-Hutch","DS-Hutch++");


