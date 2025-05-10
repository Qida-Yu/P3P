function [res_out,position_error,rotation_error]=cal_pose_error_noise(res,R,t)
%find min error with noise
min_error=inf;
Quaternion_standard = getRotationQuaternion(R);
for i=1:length(res)
    Quaternion_temp = getRotationQuaternion(res{i}.R);
    temp=norm((Quaternion_standard-Quaternion_temp),2);
    if temp<min_error
        min_error=temp;
        res_out.R=res{i}.R;
        res_out.t=res{i}.T;
    end
end 
Quaternion_estimation= getRotationQuaternion(res_out.R);
rotation_error = norm((Quaternion_standard-Quaternion_estimation),2)/pi*180;
position_error=(norm(res_out.t-t,2)/norm(t,2))*100;
end