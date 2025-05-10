function [res_out,position_error,rotation_error]=cal_pose_error(res,R,t)
%find min error without noise
min_error=inf;
for i=1:length(res)
    temp=norm(res{i}.R-R,2);
    if temp<min_error
        min_error=temp;
        res_out.R=res{i}.R;
        res_out.t=res{i}.T;
    end
end
rotation_error=norm(res_out.R-R,2);
position_error=norm(res_out.t-t,2);
end