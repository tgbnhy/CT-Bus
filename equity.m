%input the path and graph, generate the equity score composed of 
%it computes all the exact score, and lower bound for the compuatiton
function [frequency, conn, equity] = equity(A, path, n, weight, connectivity_ub, distance_ub, b, rpes, iter, current_fre, newfre, current_conn, newconn,full_computation,base)
    A_temp = A;%construct the new edge by check the path
    if size(path,2)>1
        for i=1:size(path,2)-1
            A_temp(path(i),path(i+1)) = 1;%update for a new edge
            A_temp(path(i+1),path(i)) = 1;
        end
    end
    conn = current_conn + newconn/connectivity_ub; %in case 
    if full_computation
        conn = (natural_connectivity(A_temp, n, b, rpes, iter) - base)/connectivity_ub;
    end
    frequency = current_fre + newfre/distance_ub;
    equity = (1 - weight)*(conn) + weight*frequency;
end
