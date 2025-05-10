function mi = invertMatrix(m)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
det = m(1,1)*m(2,2)*m(3,3) + m(1,2)*m(2,3)*m(3,1) + m(1,3)*m(2,1)*m(3,2) - m(1,1)*m(2,3)*m(3,2) - m(1,2)*m(2,1)*m(3,3) - m(1,3)*m(2,2)*m(3,1);
for i = 1 : 3
    for j = 1 : 3
        mi(i,j) = (m(ROTF(j),ROTF(i))*m(ROTB(j),ROTB(i))- m(ROTF(j),ROTB(i))*m(ROTB(j),ROTF(i))) / det;
    end
end

end

