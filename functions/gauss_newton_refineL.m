function L = gauss_newton_refineL(L,a12,a13,a23,b12,b13,b23)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
refinement_iterations = 2;
for i = 1 : refinement_iterations
    l1 = L(1,1);
    l2 = L(1,2);
    l3 = L(1,3);
    r1 = l1*l1 + l2*l2 +b12*l1*l2 -a12;
    r2 = l1*l1 + l3*l3 +b13*l1*l3 -a13;
    r3 = l2*l2 + l3*l3 +b23*l2*l3 -a23;
    if ((abs(r1) + abs(r2) + abs(r3)) < 1e-10)
        break;
    end
    dr1dl1 = 2*l1 + b12*l2;
    dr1dl2 = 2*l2 + b12*l1;
    
    dr2dl1 = 2*l1 + b13*l3;
    dr2dl3 = 2*l3 + b13*l1;
    
    dr3dl2 = 2*l2 + b23*l3;
    dr3dl3 = 2*l3 + b23*l2;
    
    r = [r1,r2,r3];
    v0 = dr1dl1;
    v1 = dr1dl2;
    v3 = dr2dl1;
    v5 = dr2dl3;
    v7 = dr3dl2;
    v8 = dr3dl3;
    det = ((-v0)*v5*v7 - v1*v3*v8);
    Ji = [(-v5)*v7,(-v1)*v8,v1*v5;(-v3)*v8,v0*v8,(-v0)*v5;v3*v7,(-v0)*v7,(-v1)*v3];
    r_tmp = (r*Ji.') / det;
    L1 = L - r_tmp;
    
    l1=L1(1,1);
    l2=L1(1,2);
    l3=L1(1,3);
    r11 = SQR(l1) + SQR(l2) +b12*l1*l2 -a12;
    r12 = SQR(l1) + SQR(l3) +b13*l1*l3 -a13;
    r13 = SQR(l2) + SQR(l3) +b23*l2*l3 -a23;
    if ((abs(r11) + abs(r12) + abs(r13)) > (abs(r1) + abs(r2) + abs(r3)))
        break;
    else
        L = L1;
    end
end
    

end

