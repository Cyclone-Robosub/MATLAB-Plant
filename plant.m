function plant_simulation()
    % Constants
    rev_pulse = 1100 * 1000;
    stop_pulse = 1500 * 1000;
    fwd_pulse_raw = 1900 * 1000;
    rev_adj = 0.97;
    fwd_pulse = int32(fwd_pulse_raw * rev_adj);
    frequency = 100;
    
    zero_set = zeros(1, 8);
    stop_set = stop_pulse * ones(1, 8);
    fwd_set = [stop_pulse * ones(1, 4), fwd_pulse * ones(1, 4)];
    crab_set = [stop_pulse * ones(1, 4), fwd_pulse, rev_pulse, rev_pulse, fwd_pulse];
    down_set = [rev_pulse * ones(1, 4), stop_pulse * ones(1, 4)];
    
    % Thruster positions
    thruster_positions = [
        0.2535, -0.2035, 0.042;
        0.2535, 0.2035, 0.042;
        -0.2545, -0.2035, 0.042;
        -0.2545, 0.2035, 0.042;
        0.1670, -0.1375, -0.049;
        0.1670, 0.1375, -0.049;
        -0.1975, -0.1165, -0.049;
        -0.1975, 0.1165, -0.049
    ];
    
    % Thruster directions
    sin45 = sin(pi / 4);
    thruster_directions = [
        0, 0, 1;
        0, 0, 1;
        0, 0, 1;
        0, 0, 1;
        -sin45, -sin45, 0;
        -sin45, sin45, 0;
        -sin45, sin45, 0;
        -sin45, -sin45, 0
    ];
    
    % Compute thruster torques
    thruster_torques = cross(thruster_positions, thruster_directions);
    
    % Compute wrench matrix
    wrench_matrix_transposed = [thruster_directions, thruster_torques]';
    wrench_matrix = wrench_matrix_transposed';
    
    % Simulation loop setup
    time_steps = 0:1/frequency:10; % 10 seconds simulation
    position_log = zeros(length(time_steps), 6);
    velocity_log = zeros(length(time_steps), 6);
    
    % Initial conditions
    position = zeros(1, 6);
    velocity = zeros(1, 6);
    
    for i = 1:length(time_steps)
        % Compute forces and update dynamics
        force = compute_thruster_force(stop_set, wrench_matrix);
        acceleration = force / 5.51; % Assuming mass = 5.51 kg
        
        velocity = velocity + acceleration / frequency;
        position = position + velocity / frequency;
        
        % Log data
        position_log(i, :) = position;
        velocity_log(i, :) = velocity;
    end
    
    % Plot results
    figure;
    plot(time_steps, position_log);
    title('Position Over Time');
    xlabel('Time (s)');
    ylabel('Position');
    legend('X', 'Y', 'Z', 'Roll', 'Pitch', 'Yaw');
    
    figure;
    plot(time_steps, velocity_log);
    title('Velocity Over Time');
    xlabel('Time (s)');
    ylabel('Velocity');
    legend('X', 'Y', 'Z', 'Roll', 'Pitch', 'Yaw');
end

function force = compute_thruster_force(pwm_set, wrench_matrix)
    % Define a sample force conversion function
    force_scalars = arrayfun(@pwm_force_scalar, pwm_set);
    force = force_scalars * wrench_matrix;
end

function force = pwm_force_scalar(pwm)
    pwm = pwm / 1000;
    if pwm >= 1100 && pwm < 1460
        force = (-1.24422882971549e-8) * pwm^3 + (4.02057100632393e-5) * pwm^2 - 0.0348619861030835 * pwm + 3.90671429105423;
    elseif pwm >= 1460 && pwm <= 1540
        force = 0;
    elseif pwm > 1540 && pwm <= 1900
        force = (-1.64293565374284e-8) * pwm^3 + (9.45962838560648e-5) * pwm^2 - 0.170812079190679 * pwm + 98.7232373648272;
    else
        error('PWM value out of valid range (1100-1900)');
    end
end
