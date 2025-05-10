function [R_gt, t_gt] = generate_test_cases(num_trials)
    % Generate ground truth rotation matrices and translation vectors
    R_gt = zeros(3, 3, num_trials);
    t_gt = zeros(3, num_trials);
    
    for i = 1:num_trials
        % Generate a random rotation matrix
        R_gt(:, :, i) = generate_random_rotation();
        
        % Generate a random translation vector (normalized)
        t_gt(:, i) = randn(3, 1); 
        t_gt(:, i) = t_gt(:, i) / norm(t_gt(:, i));
    end
end

function R = generate_random_rotation()
    % Generate a random rotation matrix
    [Q, ~] = qr(randn(3)); % QR decomposition to get orthogonal matrix
    R = Q * diag([1, 1, det(Q)]); % Ensure it's a proper rotation matrix
end
