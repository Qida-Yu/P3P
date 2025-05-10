function [e1,e2] = root2real(b,c)
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
v = SQR(b) - 4*c;
if (v < 0)
    e1 = (0.5)*b;
    e2 = (0.5)*b;
else
    y = sqrt(v);
    if (b < 0)
        e1 = 0.5*(-b + y);
        e2 = 0.5*(-b - y);
    else
        e1 = 2*c / (-b + y);
        e2 = 2*c / (-b - y);
    end
end

end

