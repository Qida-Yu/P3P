function v = multMatrixVector(m,u)
%UNTITLED9 此处显示有关此函数的摘要
%   此处显示详细说明
for i = 1 : 3
    v(i)  = 0;
    for j = 1 :3
        v(i) = v(i) + m(i,j) * u(j);
    end
end

end

