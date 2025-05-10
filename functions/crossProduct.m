function w = crossProduct(u,v)
for i = 1:3
    w(1,i) = u(1,ROTF(i)) * v(1,ROTB(i)) - u(1,ROTB(i)) * v(1,ROTF(i));
end
end

