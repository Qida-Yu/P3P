function [valid, unique, duplicate, ground_match, incorrect_sol] = evaluate_solutions(solutions, gt, error_threshold, duplicate_threshold)
    % 初始化指标
    valid = 0;
    unique = 0;
    duplicate = 0;
    ground_match = 0;
    incorrect_sol = 0;

    % 检查解是否为空
    if isempty(solutions)
        fprintf('No solutions provided for evaluation.\n');
        return;
    end

    % 打印解的数量
    num_solutions = numel(solutions);
    fprintf('Number of solutions to evaluate: %d\n', num_solutions);

    % 遍历解并评估
    for i = 1:num_solutions
        solution = solutions{i};
        fprintf('Evaluating solution %d:\n', i);

        % 打印解的内容
        disp('Solution content:');
        disp(solution);

        % 检查是否包含 R 和 T（字段名大写）
        if ~isfield(solution, 'R') || ~isfield(solution, 'T')
            fprintf('Solution %d is missing required fields (R or T).\n', i);
            incorrect_sol = incorrect_sol + 1;
            continue;
        end

        % 提取 R 和 T
        R = solution.R;
        t = solution.T;

        % 计算与地面真值的误差
        R_error = norm(R - gt.R, 'fro'); % Frobenius 范数
        t_error = norm(t - gt.t);       % 欧几里得距离

        fprintf('Solution %d Rotation error: %f\n', i, R_error);
        fprintf('Solution %d Translation error: %f\n', i, t_error);

        % 判断是否为有效解
        if R_error < error_threshold && t_error < error_threshold
            valid = valid + 1;
            if R_error < error_threshold && t_error < error_threshold
    ground_match = ground_match + 1;
            end

        else
            incorrect_sol = incorrect_sol + 1;
        end

        % 检测重复解
        is_duplicate = false;
        for j = 1:i-1
            prev_R = solutions{j}.R;
            prev_t = solutions{j}.T;
            if norm(prev_R - R, 'fro') < duplicate_threshold && norm(prev_t - t) < duplicate_threshold
                is_duplicate = true;
                duplicate = duplicate + 1;
                fprintf('Solution %d is a duplicate of solution %d.\n', i, j);
                break;
            end
        end

        if ~is_duplicate
            unique = unique + 1;
            fprintf('Solution %d is unique.\n', i);
        end
    end
end
