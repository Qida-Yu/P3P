%Author: Yu, Qida
function y = polishRoots(U,s_root)
% use Gauss-Newton's method to polish roots for a quartic equation
iteraterations = 10;
for i = 1 : iteraterations
    for j = 1 : length(s_root)
        err = (((U(1,1) * s_root(j,1) + U(1,2)) * s_root(j,1) + U(1,3)) * s_root(j,1) + U(1,4)) * s_root(j,1) + U(1,5); 
        derivative = ((4 * U(1,1) * s_root(j,1) + 3 * U(1,2)) * s_root(j,1) + 2 * U(1,3)) * s_root(j,1) + U(1,4);
        s_root(j,1) = s_root(j,1) - err / derivative;
    end
end
y = s_root;
end

