%Author: Yu, Qida
function q = getRandomUnitVector( k )
%get k-dimension random unit vector
q = randn(1,k);
if (sum(abs(q)) < 1e-9)
    q = getRandomUnitVector( k );
end
q = q  / norm(q);
end

