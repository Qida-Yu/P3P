%Author: Yu, Qida
function R = getRotationMatrix( q )
%get rotation martix from a quaternion
   aa = q(1,1) * q(1,1);
   ab = q(1,1) * q(1,2);
   ac = q(1,1) * q(1,3);
   ad = q(1,1) * q(1,4);
   bb = q(1,2) * q(1,2);
   bc = q(1,2) * q(1,3);
   bd = q(1,2) * q(1,4);
   cc = q(1,3) * q(1,3);
   cd = q(1,3) * q(1,4);
   dd = q(1,4) * q(1,4);
   
   a=aa + bb - cc - dd;
   b=(2.0) * (ad + bc);
   c=(2.0) * (bd - ac);
   d=(2.0) * (bc - ad);
   e=aa - bb + cc - dd;
   f=(2.0) * (ab + cd);
   g=(2.0) * (ac + bd);
   h=(2.0) * (cd - ab);
   i=aa - bb - cc + dd;
   
   R=[a, d, g;b,e,h;c,f,i];

end

