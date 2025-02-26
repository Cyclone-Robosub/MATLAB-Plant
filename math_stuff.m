function plot_pos_vel(V_i, S_i, m, C, T, dt, t)
    % Constants
    RHO = 1000;
    GRAVITY = 9.81;
    
    % Time steps
    time_steps = 0:dt:t;
    positions = zeros(size(time_steps));
    velocities = zeros(size(time_steps));
    
    S = S_i;
    V = V_i;
    
    for i = 1:length(time_steps)
        [S, V] = pos_vel(V, S, m, C, T, dt);
        positions(i) = S;
        velocities(i) = V;
    end
    
    % Plot results
    figure;
    yyaxis left;
    plot(time_steps, positions, 'b');
    ylabel('Position');
    
    yyaxis right;
    plot(time_steps, velocities, 'r--');
    ylabel('Velocity');
    xlabel('Time (s)');
    title('Position and Velocity Over Time');
end

function R = rotation_matrix(roll, pitch, yaw)
    % Rotation matrices
    Rx = [1, 0, 0; 0, cos(roll), -sin(roll); 0, sin(roll), cos(roll)];
    Ry = [cos(pitch), 0, sin(pitch); 0, 1, 0; -sin(pitch), 0, cos(pitch)];
    Rz = [cos(yaw), -sin(yaw), 0; sin(yaw), cos(yaw), 0; 0, 0, 1];
    R = Rz * Ry * Rx;
end

function angle = wrap_to_pi(angle)
    % Wrap an angle to (-pi, pi)
    angle = mod(angle + pi, 2*pi) - pi;
    if angle == -pi
        angle = pi;
    end
end

function result = wrap_array_to_pi(angles)
    result = arrayfun(@wrap_to_pi, angles);
end

function F_rotated = weight_force(roll, pitch, yaw, mass, g)
    if nargin < 5, g = 9.81; end
    g_vector = [0; 0; -mass * g];
    R = rotation_matrix(roll, pitch, yaw);
    F_rotated = R * g_vector;
end

function output = buoyant_force(roll, pitch, yaw, volume, d, rho, g)
    if nargin < 6, rho = 1000; end
    if nargin < 7, g = 9.81; end
    d = d(:);
    b_vector = [0; 0; rho * g * volume];
    R = rotation_matrix(roll, pitch, yaw);
    buoyant_force_body = R * b_vector;
    torque = cross(-d, buoyant_force_body);
    output = [buoyant_force_body; torque];
end

function [S, V] = pos_vel(V_i, S_i, m, C, T, dt)
    if C < 0 || m <= 0 || dt <= 0
        error('Invalid input');
    elseif T == 0
        [S, V] = pos_vel_case_1(V_i, S_i, m, C, T, dt);
    elseif V_i == 0
        [S, V] = pos_vel_case_2(V_i, S_i, m, C, T, dt);
    elseif (V_i > 0 && T > 0) || (V_i < 0 && T < 0)
        [S, V] = pos_vel_case_2(V_i, S_i, m, C, T, dt);
    elseif (V_i < 0 && T > 0) || (V_i > 0 && T < 0)
        [S, V] = pos_vel_case_3(V_i, S_i, m, C, T, dt);
    else
        error('Invalid input');
    end
end

function [S, V] = pos_vel_case_1(V_i, S_i, m, C, T, dt)
    if T == 0 && V_i == 0
        S = S_i;
        V = 0;
    elseif T == 0 && V_i > 0
        S = S_i - (m / C) * log(m / abs(-C * V_i * dt - m));
        V = - (-1 / V_i - C * dt / m) ^ (-1);
    elseif T == 0 && V_i < 0
        S = -S_i - (m / C) * log(m / abs(C * V_i * dt - m));
        V = - (1 / V_i - C * dt / m) ^ (-1);
        S = -S;
        V = -V;
    else
        error('Invalid input');
    end
end

function [S, V] = pos_vel_case_2(V_i, S_i, m, C, T, dt)
    a = sqrt(T / C);
    if V_i == a
        S = S_i + dt * V_i;
        V = V_i;
    else
        A = (a + V_i) / (a - V_i);
        B = 2 * a * C / m;
        D = abs(A * exp(B * dt) + 1);
        E = 2 * log(D) / B - dt;
        F = 2 * log(abs(A + 1)) / B;
        S = S_i + a * E - a * F;
        num3 = A * exp(B * dt) - 1;
        den3 = A * exp(B * dt) + 1;
        V = a * num3 / den3;
    end
end

function [S, V] = pos_vel_case_3(V_i, S_i, m, C, T, dt)
    a = sqrt(C / m);
    A = atan(V_i / a);
    B = (C * a * dt / m) + A;
    D = m * log(abs(1 / cos(B))) / C;
    E = m * log(abs(1 / cos(A))) / C;
    if T > 0 && V_i < 0
        S = S_i + D - E;
        V = a * tan(A + a * C * dt / m);
    elseif T < 0 && V_i > 0
        [S, V] = pos_vel_case_3(-V_i, -S_i, m, C, -T, dt);
        S = -S;
        V = -V;
    else
        error('Invalid input');
    end
end
