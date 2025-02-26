classdef Plant < handle
    %% A MATLAB class replicating all functionality of the original Python code.
    %% This class includes:
    %  - Thruster constants (rev_pulse, stop_pulse, etc.).
    %  - Data arrays (e.g., stop_set, crab_set, etc.).
    %  - The full simulation logic with weight, buoyant, drag, and thruster forces.
    %  - The pos_vel, rotation_matrix, and wrap methods from the Python math_stuff.
    %  - Logging states in a structure array (state_log) for each time step.
    %  - Graphing methods for acceleration, velocity, position, forces, and PWM signals.

    properties (Constant)
        %% Thruster PWM Constants
        rev_pulse       = 1100 * 1000;
        stop_pulse      = 1500 * 1000;
        fwd_pulse_raw   = 1900 * 1000;
        rev_adj         = 0.97;  % Because forward is stronger than reverse
        frequency       = 100;   % Not used directly, but kept for reference
        pwm_file        = 'pwm_file.csv';  % Unused in code, but kept for completeness
    end

    properties
        %% Derived Thruster Constants
        fwd_pulse       % Forward pulse factoring in rev_adj
        
        zero_set        % 8-element array of zeros
        stop_set        % 8-element array of stop_pulse
        fwd_set         % Four stop, four forward
        crab_set        % Four stop, then forward, reverse, reverse, forward
        down_set        % Four reverse, four stop
        barrel_set      % [rev, fwd, rev, fwd], then four stop
        summer_set      % [rev, rev, fwd, fwd], then four stop
        spin_set        % Four stop, then [fwd, rev, fwd, rev]
        torpedo_set     % [rev, fwd, rev, fwd], then four forward

        %% Plant Physical Constants and States
        default_frequency = 100;      % Simulation steps per second
        mass = 5.51;                 % In kilograms
        height = 0.3;                % Not used in code, but kept
        mass_moment_of_inertia = [1, 1, 1];  % 3D inertia
        six_axis_mass;               % [m, m, m, Ixx, Iyy, Izz]

        volume_inches = 449.157;     % Volume in cubic inches
        volume;                      % Volume in cubic meters
        rho_water = 1000;           % Density of water (kg/m^3)

        combined_drag_coefs = [0.041, 0.05, 0.125, 0.005, 0.005, 0.005];

        mass_center_inches = [0, 0, 0];
        mass_center;                % In meters
        volume_center = [0, 0, 0.1]; % Dist from COM to COB

        %% Thrusters
        thruster_positions;         % 8x3
        thruster_directions;        % 8x3
        thruster_torques;           % 8x3 cross product of positions, directions
        wrench_matrix;              % 6x8 net directions + torques

        %% Current PWM Commands
        next_pwm;                   % 1x8 array for next time step

        %% Log of states, as an array of structures. Each entry:
        %   time, position(6), velocity(6), acceleration(6),
        %   totalForces(6), weightForces(6), buoyantForces(6), pwm(1x8)
        state_log;                  
    end

    methods
        %% Constructor
        function obj = Plant()
            % Derived forward pulse factoring in rev_adj
            obj.fwd_pulse = floor(obj.fwd_pulse_raw * obj.rev_adj);

            % Expand combined_drag_coefs by factor of 10, as in Python code
            obj.combined_drag_coefs = obj.combined_drag_coefs * 10;

            % Initialize the sets used in Python
            obj.zero_set    = zeros(1, 8, 'int32');
            obj.stop_set    = ones(1, 8, 'int32') * obj.stop_pulse;
            obj.fwd_set     = [ones(1, 4, 'int32') * obj.stop_pulse, ones(1, 4, 'int32') * obj.fwd_pulse];
            obj.crab_set    = [ones(1, 4, 'int32') * obj.stop_pulse, int32([obj.fwd_pulse, obj.rev_pulse, obj.rev_pulse, obj.fwd_pulse])];
            obj.down_set    = [ones(1, 4, 'int32') * obj.rev_pulse, ones(1, 4, 'int32') * obj.stop_pulse];
            obj.barrel_set  = [int32([obj.rev_pulse, obj.fwd_pulse, obj.rev_pulse, obj.fwd_pulse]), ones(1, 4, 'int32') * obj.stop_pulse];
            obj.summer_set  = [int32([obj.rev_pulse, obj.rev_pulse, obj.fwd_pulse, obj.fwd_pulse]), ones(1, 4, 'int32') * obj.stop_pulse];
            obj.spin_set    = [ones(1, 4, 'int32') * obj.stop_pulse, int32([obj.fwd_pulse, obj.rev_pulse, obj.fwd_pulse, obj.rev_pulse])];
            obj.torpedo_set = [int32([obj.rev_pulse, obj.fwd_pulse, obj.rev_pulse, obj.fwd_pulse]), ones(1, 4, 'int32') * obj.fwd_pulse];

            % Compute volume in m^3
            obj.volume = obj.volume_inches * (0.0254)^3;

            % The 6-axis mass = [mx, my, mz, Ixx, Iyy, Izz]
            obj.six_axis_mass = zeros(1, 6);
            obj.six_axis_mass(1:3) = obj.mass;  % same mass in x,y,z
            obj.six_axis_mass(4:6) = obj.mass_moment_of_inertia;  % inertias

            % Convert center of mass to meters
            obj.mass_center = obj.mass_center_inches * 0.0254;

            % Thruster positions
            obj.thruster_positions = [
                0.2535, -0.2035, 0.042;
                0.2535,  0.2035, 0.042;
               -0.2545, -0.2035, 0.042;
               -0.2545,  0.2035, 0.042;
                0.1670, -0.1375,-0.049;
                0.1670,  0.1375,-0.049;
               -0.1975, -0.1165,-0.049;
               -0.1975,  0.1165,-0.049
            ];

            % Thruster directions
            sin45 = sin(pi/4);
            obj.thruster_directions = [
                0,      0,      1;
                0,      0,      1;
                0,      0,      1;
                0,      0,      1;
               -sin45, -sin45,  0;
               -sin45,  sin45,  0;
               -sin45,  sin45,  0;
               -sin45, -sin45,  0
            ];

            % Thruster torques
            obj.thruster_torques = zeros(8,3);
            for i=1:8
                obj.thruster_torques(i,:) = cross(obj.thruster_positions(i,:), obj.thruster_directions(i,:));
            end

            % Combine directions + torques => 6x8 matrix
            % Actually, we store as (8 x 6) then transpose.
            big = [obj.thruster_directions, obj.thruster_torques]'; % shape (6x8)
            obj.wrench_matrix = big;  % (6x8)

            % Start with stop commands
            obj.next_pwm = obj.stop_set;

            % Initialize state_log as an array of structs. We'll store the first entry.
            s.time          = 0;
            s.position      = zeros(1, 6);
            s.velocity      = zeros(1, 6);
            s.acceleration  = zeros(1, 6);
            s.totalForces   = zeros(1, 6);
            s.weightForces  = zeros(1, 6);
            s.buoyantForces = zeros(1, 6);
            s.pwm           = obj.stop_set;

            obj.state_log   = s;
        end

        %% Accessor methods (Python-like)
        function t = time(obj)
            t = obj.state_log(end).time;
        end
        function r = roll(obj)
            r = obj.state_log(end).position(4);
        end
        function p = pitch(obj)
            p = obj.state_log(end).position(5);
        end
        function y = yaw(obj)
            y = obj.state_log(end).position(6);
        end
        function pos = position(obj)
            pos = obj.state_log(end).position;
        end
        function vel = velocity(obj)
            vel = obj.state_log(end).velocity;
        end
        function acc = acceleration(obj)
            acc = obj.state_log(end).acceleration;
        end
        function tf = total_forces(obj)
            tf = obj.state_log(end).totalForces;
        end
        function wf = weight_forces(obj)
            wf = obj.state_log(end).weightForces;
        end
        function bf = buoyant_forces(obj)
            bf = obj.state_log(end).buoyantForces;
        end
        function p = pwms(obj)
            p = obj.state_log(end).pwm;
        end

        %% PWM force scalar function
        function force = pwm_force_scalar(~, x)
            % x in microseconds from 1100 to 1900, e.g. 1500e3 => 1500000.
            % In Python code, x was divided by 1000, so let's do same.
            val = x / 1000.0; % typical range 1100..1900 => 1.1..1.9k?
            if (val >= 1100 && val < 1460)
                force = (-1.24422882971549e-8) * val^3 + (4.02057100632393e-5) * val^2 - 0.0348619861030835 * val + 3.90671429105423;
            elseif (val >= 1460 && val <= 1540)
                force = 0;
            elseif (val > 1540 && val <= 1900)
                force = (-1.64293565374284e-8) * val^3 + (9.45962838560648e-5) * val^2 - 0.170812079190679 * val + 98.7232373648272;
            else
                error('PWM value %f out of valid range (1100..1900)', val);
            end
        end

        %% Summed thruster force across all thrusters
        function force = pwm_force(obj, pwm_set)
            thruster_forces = arrayfun(@(p) obj.pwm_force_scalar(p), pwm_set);
            force = thruster_forces * obj.wrench_matrix';  % thruster_forces(1x8) * (8x6) => (1x6)
            % But in Python code, we had (thruster_forces) as shape (8,) and wrench_matrix as shape (6,8).
            % So the result was shape (6,). We'll just transpose to get (1x6).
        end

        %% Set PWM commands for next iteration
        function obj = set_pwm(obj, pwm_set)
            obj.next_pwm = pwm_set;
        end

        %% rotation_matrix (roll, pitch, yaw)
        function R = rotation_matrix(~, roll, pitch, yaw)
            Rx = [1, 0, 0; 0, cos(roll), -sin(roll); 0, sin(roll), cos(roll)];
            Ry = [cos(pitch), 0, sin(pitch); 0, 1, 0; -sin(pitch), 0, cos(pitch)];
            Rz = [cos(yaw), -sin(yaw), 0; sin(yaw), cos(yaw), 0; 0, 0, 1];
            R = Rz * Ry * Rx;
        end

        %% wrap_to_pi(angle)
        function ang = wrap_to_pi(~, angle)
            ang = mod(angle + pi, 2*pi) - pi;
            % If exactly -pi, return pi
            if abs(ang + pi) < 1e-12
                ang = pi;
            end
        end

        function arr = wrap_array_to_pi(obj, angles)
            arr = zeros(size(angles));
            for i = 1:length(angles)
                arr(i) = obj.wrap_to_pi(angles(i));
            end
        end

        %% weight_force(roll, pitch, yaw, mass, g=9.81)
        function w = weight_force(obj)
            g = 9.81;
            mass_ = obj.mass;
            g_vector = [0; 0; -mass_*g];
            R = obj.rotation_matrix(obj.roll(), obj.pitch(), obj.yaw());
            F_rotated = R * g_vector;  % shape (3,1)
            w = [F_rotated' 0 0 0];    % shape (1,6), only linear forces
        end

        %% buoyant_force(roll, pitch, yaw, volume, d, rho=1000, g=9.81)
        function b = buoyant_force(obj)
            rho = obj.rho_water;
            g = 9.81;
            vol = obj.volume;
            d = -obj.volume_center;  % from code -> buoyant_force(self.roll(), self.pitch(), self.yaw(), self.volume, -1 * self.volume_center)

            b_vector = [0; 0; rho*g*vol];
            R = obj.rotation_matrix(obj.roll(), obj.pitch(), obj.yaw());
            buoyant_force_body = R * b_vector;  % (3,1)
            torque = cross(d, buoyant_force_body');

            % shape: combine force + torque => (1,6)
            b = [buoyant_force_body' torque];
        end

        %% drag_force
        function d = drag_force(obj)
            d = zeros(1,6);
            v = obj.velocity();
            C = obj.combined_drag_coefs;
            % For i in range(6): drag_force[i] = - velocity[i]*abs(velocity[i]) * drag_coefs[i]
            for i=1:6
                d(i) = -v(i)*abs(v(i))*C(i);
            end
        end

        %% total_force
        function tf = total_force(obj)
            w = obj.weight_force();
            b = obj.buoyant_force();
            t = obj.pwm_force(obj.pwms());
            d = obj.drag_force();
            tf = w + b + t + d;
        end

        %% step() => updates the state by one dt
        function obj = step(obj)
            dt = 1 / obj.default_frequency;

            new_state = obj.state_log(end);
            new_state.time = new_state.time + dt;

            pos_old = obj.state_log(end).position;
            vel_old = obj.state_log(end).velocity;

            % We need to do per-axis integrator using the pos_vel approach from math_stuff.
            % That code calls pos_vel(V_i, S_i, m, C, T, dt). We'll replicate it below.

            new_position = zeros(1,6);
            new_velocity = zeros(1,6);
            new_acc      = zeros(1,6);

            totalF = obj.total_force();
            wF     = obj.weight_force();
            bF     = obj.buoyant_force();

            for i = 1:6
                Si = pos_old(i);
                Vi = vel_old(i);
                m  = obj.six_axis_mass(i);
                C  = obj.combined_drag_coefs(i);
                T  = totalF(i);

                % We call pos_vel(Vi, Si, m, C, T, dt)
                [S, V] = obj.pos_vel(Vi, Si, m, C, T, dt);

                new_position(i) = S;
                new_velocity(i) = V;
                new_acc(i)      = T/m;
            end

            new_state.position = new_position;
            new_state.velocity = new_velocity;
            new_state.acceleration = new_acc;
            new_state.totalForces = totalF;
            new_state.weightForces = wF;
            new_state.buoyantForces = bF;
            new_state.pwm = obj.next_pwm;

            % Append to the log
            obj.state_log(end+1) = new_state;
        end

        %% run_pwm
        function obj = run_pwm(obj, pwm_set, time_sec)
            obj.set_pwm(pwm_set);
            steps = round(time_sec * obj.default_frequency);
            for i=1:steps
                obj = obj.step();
            end
        end

        %% print_dictionary => print the entire state_log
        function print_dictionary(obj)
            for i=1:length(obj.state_log)
                st = obj.state_log(i);
                fprintf('Time: %.3f\n', st.time);
                fprintf('Position: ['); fprintf(' %.3f', st.position); fprintf(' ]\n');
                fprintf('Velocity: ['); fprintf(' %.3f', st.velocity); fprintf(' ]\n');
                fprintf('Acceleration: ['); fprintf(' %.3f', st.acceleration); fprintf(' ]\n');
                fprintf('Total Forces: ['); fprintf(' %.3f', st.totalForces); fprintf(' ]\n');
                fprintf('Weight Forces: ['); fprintf(' %.3f', st.weightForces); fprintf(' ]\n');
                fprintf('Buoyant Forces: ['); fprintf(' %.3f', st.buoyantForces); fprintf(' ]\n');
                fprintf('PWM: ['); fprintf(' %d', st.pwm); fprintf(' ]\n\n');
            end
        end

        %% Graphing methods
        function graph_acceleration(obj)
            % gather data from log
            dt = 1/obj.default_frequency;
            time_vals = dt * (0:(length(obj.state_log)-1));

            accel_mat = reshape([obj.state_log.acceleration], 6, [])';
            figure;
            subplot(2,1,1);
            plot(time_vals, accel_mat(:,1), 'b', time_vals, accel_mat(:,2), 'r', time_vals, accel_mat(:,3), 'g');
            xlabel('Time (s)'); ylabel('m/s^2');
            legend('Ax','Ay','Az'); title('Linear Acceleration');

            subplot(2,1,2);
            plot(time_vals, accel_mat(:,4), 'b', time_vals, accel_mat(:,5), 'r', time_vals, accel_mat(:,6), 'g');
            xlabel('Time (s)'); ylabel('rad/s^2');
            legend('Alpha Roll','Alpha Pitch','Alpha Yaw'); title('Rotational Acceleration');
        end

        function graph_velocity(obj)
            dt = 1/obj.default_frequency;
            time_vals = dt * (0:(length(obj.state_log)-1));

            vel_mat = reshape([obj.state_log.velocity], 6, [])';
            figure;
            subplot(2,1,1);
            plot(time_vals, vel_mat(:,1), 'b', time_vals, vel_mat(:,2), 'r', time_vals, vel_mat(:,3), 'g');
            xlabel('Time (s)'); ylabel('m/s');
            legend('Vx','Vy','Vz'); title('Linear Velocity');

            subplot(2,1,2);
            plot(time_vals, vel_mat(:,4), 'b', time_vals, vel_mat(:,5), 'r', time_vals, vel_mat(:,6), 'g');
            xlabel('Time (s)'); ylabel('rad/s');
            legend('Roll Rate','Pitch Rate','Yaw Rate'); title('Rotational Velocity');
        end

        function graph_position(obj)
            dt = 1/obj.default_frequency;
            time_vals = dt * (0:(length(obj.state_log)-1));

            pos_mat = reshape([obj.state_log.position], 6, [])';
            figure;
            subplot(2,1,1);
            plot(time_vals, pos_mat(:,1), 'b', time_vals, pos_mat(:,2), 'r', time_vals, pos_mat(:,3), 'g');
            xlabel('Time (s)'); ylabel('m');
            legend('X','Y','Z'); title('Linear Position');

            subplot(2,1,2);
            plot(time_vals, pos_mat(:,4), 'b', time_vals, pos_mat(:,5), 'r', time_vals, pos_mat(:,6), 'g');
            xlabel('Time (s)'); ylabel('rad');
            legend('Roll','Pitch','Yaw'); title('Rotational Position');
        end

        function graph_total_forces(obj)
            dt = 1/obj.default_frequency;
            time_vals = dt * (0:(length(obj.state_log)-1));

            tf_mat = reshape([obj.state_log.totalForces], 6, [])';
            figure;
            subplot(2,1,1);
            plot(time_vals, tf_mat(:,1), 'b', time_vals, tf_mat(:,2), 'r', time_vals, tf_mat(:,3), 'g');
            xlabel('Time (s)'); ylabel('N');
            legend('Fx','Fy','Fz'); title('Linear Total Forces');

            subplot(2,1,2);
            plot(time_vals, tf_mat(:,4), 'b', time_vals, tf_mat(:,5), 'r', time_vals, tf_mat(:,6), 'g');
            xlabel('Time (s)'); ylabel('N*m');
            legend('Mx','My','Mz'); title('Rotational Total Forces');
        end

        function graph_weight_forces(obj)
            dt = 1/obj.default_frequency;
            time_vals = dt * (0:(length(obj.state_log)-1));

            wf_mat = reshape([obj.state_log.weightForces], 6, [])';
            figure;
            subplot(2,1,1);
            plot(time_vals, wf_mat(:,1), 'b', time_vals, wf_mat(:,2), 'r', time_vals, wf_mat(:,3), 'g');
            xlabel('Time (s)'); ylabel('N');
            legend('Wx','Wy','Wz'); title('Linear Weight Forces');

            subplot(2,1,2);
            plot(time_vals, wf_mat(:,4), 'b', time_vals, wf_mat(:,5), 'r', time_vals, wf_mat(:,6), 'g');
            xlabel('Time (s)'); ylabel('N*m');
            legend('Mx','My','Mz'); title('Rotational Weight Forces');
        end

        function graph_buoyant_forces(obj)
            dt = 1/obj.default_frequency;
            time_vals = dt * (0:(length(obj.state_log)-1));

            bf_mat = reshape([obj.state_log.buoyantForces], 6, [])';
            figure;
            subplot(2,1,1);
            plot(time_vals, bf_mat(:,1), 'b', time_vals, bf_mat(:,2), 'r', time_vals, bf_mat(:,3), 'g');
            xlabel('Time (s)'); ylabel('N');
            legend('Bx','By','Bz'); title('Linear Buoyant Forces');

            subplot(2,1,2);
            plot(time_vals, bf_mat(:,4), 'b', time_vals, bf_mat(:,5), 'r', time_vals, bf_mat(:,6), 'g');
            xlabel('Time (s)'); ylabel('N*m');
            legend('Mx','My','Mz'); title('Rotational Buoyant Forces');
        end

        function graph_pwm_signals(obj)
            dt = 1/obj.default_frequency;
            time_vals = dt * (0:(length(obj.state_log)-1));

            % Each log entry has pwm(1x8). We'll plot all 8 thrusters vs time.
            all_pwm = reshape([obj.state_log.pwm], 8, [])';

            figure;
            hold on;
            for thr = 1:8
                plot(time_vals, all_pwm(:,thr));
            end
            hold off;
            xlabel('Time (s)'); ylabel('PWM (microseconds)');
            legend('T1','T2','T3','T4','T5','T6','T7','T8'); title('PWM Signals');
        end
    end

    %% Private or helper methods for the pos_vel logic
    methods (Access = private)
        function [S, V] = pos_vel(obj, V_i, S_i, m, C, T, dt)
            %% Direct port of pos_vel from Python
            if C < 0 || m <= 0 || dt <= 0
                error('Invalid input to pos_vel');
            elseif T == 0
                [S, V] = obj.pos_vel_case_1(V_i, S_i, m, C, T, dt);
            elseif V_i == 0
                [S, V] = obj.pos_vel_case_2(V_i, S_i, m, C, T, dt);
            elseif (V_i > 0 && T > 0) || (V_i < 0 && T < 0)
                [S, V] = obj.pos_vel_case_2(V_i, S_i, m, C, T, dt);
            elseif (V_i < 0 && T > 0) || (V_i > 0 && T < 0)
                [S, V] = obj.pos_vel_case_3(V_i, S_i, m, C, T, dt);
            else
                error('Invalid input in pos_vel');
            end
        end

        function [S, V] = pos_vel_case_1(~, V_i, S_i, m, C, T, dt)
            if T == 0 && abs(V_i) < 1e-12
                % V_i == 0
                S = S_i;
                V = 0;
            elseif T == 0 && V_i > 0
                S = S_i - (m / C) * log(m / abs(-C * V_i * dt - m));
                V = - (-1 / V_i - C * dt / m)^(-1);
            elseif T == 0 && V_i < 0
                S_pos = - S_i - (m / C) * log(m / abs(C * V_i * dt - m));
                V_temp = - (1 / V_i - C * dt / m)^(-1);
                % then we flip sign to replicate python logic
                S = -S_pos;
                V = -V_temp;
            else
                error('Invalid input in pos_vel_case_1');
            end
        end

        function [S, V] = pos_vel_case_2(~, V_i, S_i, m, C, T, dt)
            if T == 0
                % This can appear if V_i == 0 & T=0 => we handle it in case_1 or above
                S = S_i;
                V = V_i;
                return;
            end
            if (T > 0 && V_i >= 0)
                a = sqrt(T / C);
                if abs(V_i - a) < 1e-12
                    % at terminal velocity already
                    S = S_i + dt * V_i;
                    V = V_i;
                else
                    A = (a + V_i) / (a - V_i);
                    B = 2 * a * C / m;
                    D_ = abs(A * exp(B * dt) + 1);
                    E_ = 2 * log(D_) / B - dt;
                    F_ = 2 * log(abs(A + 1)) / B;
                    S = S_i + a * E_ - a * F_;
                    num3 = A * exp(B * dt) - 1;
                    den3 = A * exp(B * dt) + 1;
                    V = a * (num3 / den3);
                end
            elseif (T < 0 && V_i <= 0)
                % sign flip
                [S_temp, V_temp] = pos_vel_case_2([], -V_i, -S_i, m, C, -T, dt);
                S = -S_temp;
                V = -V_temp;
            else
                error('Invalid input in pos_vel_case_2');
            end
        end

        function [S, V] = pos_vel_case_3(obj, V_i, S_i, m, C, T, dt)
            a = sqrt(C / m);
            A_ = atan(V_i / a);
            B_ = (C * a * dt / m) + A_;
            D_ = m * log(abs(1 / cos(B_))) / C;
            E_ = m * log(abs(1 / cos(A_))) / C;

            if (T > 0 && V_i < 0)
                S = S_i + D_ - E_;
                V = a * tan(A_ + a * C * dt / m);
            elseif (T < 0 && V_i > 0)
                [S_temp, V_temp] = obj.pos_vel_case_3(-V_i, -S_i, m, C, -T, dt);
                S = -S_temp;
                V = -V_temp;
            else
                error('Invalid input in pos_vel_case_3');
            end
        end
    end
end
