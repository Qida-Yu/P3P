function v = multMatrixVector(m,u)
%UNTITLED9 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
for i = 1 : 3
    v(i)  = 0;
    for j = 1 :3
        v(i) = v(i) + m(i,j) * u(j);
    end
end

end

