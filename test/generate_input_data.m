function [input_data, gt] = generate_input_data(R_gt, t_gt, solver_name)
    % Number of points
    num_points = 3; % P3P 求解器需要至少 3 个点

    % 生成 3D 世界坐标点
    Xw = randn(3, num_points); % 随机生成 3D 世界坐标

    % 投影到相机坐标系
    Xc = R_gt * Xw + t_gt;

    % 归一化 2D 图像点（无噪声）
    xn = Xc(1:2, :) ./ Xc(3, :); % 转换为归一化图像平面

    % 调试信息
    fprintf('Generating input data for solver: %s (No noise added)\n', solver_name);

    % 根据不同的求解器，调整输入数据格式
    switch solver_name
        case 'Gao'
            input_data.Xw = Xw; % Gao 求解器需要 Xw 和 xn
            input_data.xn = xn;
        case 'Kneip'
            input_data.worldPoints = Xw; % Kneip 求解器需要 worldPoints 和 xn
            input_data.xn = xn;
        case {'Banno', 'Ke', 'Persson', 'Yu', 'Jrrr'}
            input_data.P = Xw; % 这些求解器需要 P 和 xn
            input_data.xn = xn;
        otherwise
            error('Unsupported solver: %s', solver_name);
    end

    % Ground truth
    gt.R = R_gt; % 旋转矩阵
    gt.t = t_gt; % 平移向量
end
