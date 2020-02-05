%compute the angle
function an = angle(x0, y0, x1, y1, x2, y2)
    P0 = [x0, y0];
    P1 = [x1, y1];
    P2 = [x2, y2];
    n1 = (P2 - P0) / norm(P2 - P0);  % Normalized vectors
    n2 = (P1 - P0) / norm(P1 - P0);
    an = atan2(norm(det([n2; n1])), dot(n1, n2));
end