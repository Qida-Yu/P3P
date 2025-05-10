function [res_out,position_error,rotation_error] = pose_error(res,R,t)

min_error=inf;
for i=1:length(res)
    temp=norm(res{i}.R-R,2);
    if temp<min_error
        min_error=temp;
        res_out.R=res{i}.R;
        res_out.t=res{i}.T;
    end
end

Angle_standard = Angle(R);
Angle_estimation= Angle(res_out.R);
rotation_error=abs((Angle_estimation-Angle_standard));
position_error=abs((res_out.t-t));
end

