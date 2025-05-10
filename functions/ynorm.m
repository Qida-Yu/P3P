% normalize each row
function V1 = ynorm( V0 )

V1 = V0 .* kron(ones(1,3), 1./sqrt(sum(V0.^2,2)));


