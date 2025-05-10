%Author: Yu, Qida
function q = getRotationQuaternion(R)
% get quaternion from a rotation martix
tr = trace(R) + 1;
if (tr > 1e-5)
    S = 0.5 / sqrt(tr);
    q(1,1) = 0.25/S;
    q(1,2) = ( R(3,2) - R(2,3) ) * S;
    q(1,3) = ( R(1,3) - R(3,1) ) * S;
    q(1,4) = ( R(2,1) - R(1,2) ) * S;
else
    if ( (R(1,1) > R(2,2)) && (R(1,1) > R(3,3)) )
        S  = sqrt( 1.0 + R(1,1) - R(2,2) - R(3,3) ) * 2.0;
        q(1,1) = (R(3,2) - R(2,3) ) / S;
        q(1,2) = 0.25 * S;
        q(1,3) = (R(2,1) + R(1,2) ) / S;
        q(1,4) = (R(1,3) + R(3,1) ) / S;
    else
        if ( R(2,2) > R(3,3) )
            S  = sqrt( 1.0 + R(2,2) - R(1,1) - R(3,3) ) * 2.0;
            q(1,1) = (R(1,3) - R(3,1) ) / S;
            q(1,2) = (R(2,1) + R(1,2) ) / S;
            q(1,3) = 0.25 * S;
            q(1,4) = (R(3,2) + R(2,3) ) / S;
        else
            S  = sqrt( 1.0 + R(3,3) - R(1,1) - R(2,2) ) * 2.0;
            q(1,1) = (R(2,1) - R(1,2) ) / S;
            q(1,2) = (R(1,3) + R(3,1) ) / S;
            q(1,3) = (R(3,2) + R(2,3) ) / S;
            q(1,4) = 0.25 * S;
        end
    end
end
%q = q / norm(q);        
        
end

