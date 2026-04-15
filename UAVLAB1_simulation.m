function UAV_work1()
clc
close all
%% =========================================================
%  STATE VECTOR
%
% x = [p; v; lambda; omega]
%
% p      = [x y z]'         position in inertial frame
% v      = [u v w]'         velocity in body frame
% lambda = [phi theta psi]' Euler angles (rotation order - ZYX)
% omega  = [p q r]'         angular velocity in body frame
%
%  4   ^    1
%      |           Motor positions
%  3        2
%
%% =========================================================


%% =========================================================

% Crazyflie 2.1 physical parameters

%% =========================================================

P.m = 0.027; % mass [kg]

P.g = 9.81; % gravity [m/s^2]

P.l = 0.046; % arm length [m]

P.J = diag([1.40e-5 1.40e-5 2.17e-5]); % inertia matrix (NOT CHECKED)

P.D = diag([0.015 0.015 0.030]); % aerodynamic drag     (NOT CHECKED)

P.kQ = 0.001; % rotor drag constant                     (NOT CHECKED)

P.e3 = [0 0 1]';

%% =========================================================

% Crazyfly Motor Constants ?? (BUSCAR FONTES PORRA)
P.Ct = 2.1302e-08;     % Thrust coefficient [N/(rad/s)^2]
P.Cq = 1.95e-09;       % Torque coefficient [Nm/(rad/s)^2]

%% =========================================================
%  Input: Motor Speeds (omega)
%% =========================================================
% Calculate the speed required for a steady hover
% T_hover = m*g / 4  =>  Ct * omega^2 = m*g / 4
omega_hover = sqrt((P.m * P.g) / (4 * P.Ct));

% Define speeds for the 4 motors (rad/s)
% Adding a tiny difference to motor 1 to see some rotation
w1 = omega_hover * 1.01; 
w2 = omega_hover * 1.01;
w3 = omega_hover * 1.00;
w4 = omega_hover * 1.00;

omega_rot = [w1; w2; w3; w4]; 

%% =========================================================
% Simulation
%% =========================================================
x0 = zeros(12,1);
Tf = 2;
tspan = [0 Tf];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);

% We pass omega_rot (input)
[t,x] = ode45(@(t,x) crazyflie_model(t,x,omega_rot,P), tspan, x0, opts);




%% =========================================================
% Plots
%% =========================================================
figure('Name','UAV States')
subplot(3,1,1)
plot(t,x(:,1:3),'LineWidth',1.5)
grid on; ylabel('Position [m]'); legend('x','y','z'); title('Inertial Position')

subplot(3,1,2)
plot(t,x(:,4:6),'LineWidth',1.5)
grid on; ylabel('Velocity [m/s]'); legend('u','v','w'); title('Body Velocity')

subplot(3,1,3)
phi   = wrap180(rad2deg(x(:,7)));
theta = wrap180(rad2deg(x(:,8)));
psi   = wrap180(rad2deg(x(:,9)));
plot(t,phi, t,theta, t,psi,'LineWidth',1.5)
ylim([-180 180]); grid on; ylabel('Euler [deg]'); 
legend('\phi','\theta','\psi'); title('Attitude (Wrapped)')

end




%% =========================================================
% NONLINEAR DYNAMICS
%% =========================================================
function dx = crazyflie_model(~,x,omega_rot,P)
    v      = x(4:6);
    lambda = x(7:9);
    omega  = x(10:12);
    
    R = rotm_zyx(lambda);
    Q = euler_rate_matrix(lambda);
    [fg, fa, fp, np] = force_moment_breakdown(x, omega_rot, P);
    
    pdot = R*v;
    vdot = -skew(omega)*v + (fg + fa + fp)/P.m;
    lambdadot = Q*omega;
    omegadot = P.J \ (-skew(omega)*P.J*omega + np);
    
    dx = [pdot; vdot; lambdadot; omegadot];
end

%% =========================================================
% FORCE AND MOMENT MODEL (SPEED AS INPUT)
%% =========================================================
function [fg, fa, fp, np] = force_moment_breakdown(x, omega_rot, P)
    v = x(4:6);
    R = rotm_zyx(x(7:9));
    
    % Gravity & Drag
    fg = -P.m*P.g*(R')*P.e3;
    fa = -P.D*v;
    
    % Motor Positions (X-configuration)
    l = P.l;
    
    % Cos/Sin(pi/4) = 0.7071
    p1 = [ l*0.7071; -l*0.7071; 0 ];
    p2 = [-l*0.7071; -l*0.7071; 0 ];
    p3 = [-l*0.7071;  l*0.7071; 0 ];
    p4 = [ l*0.7071;  l*0.7071; 0 ];
    
    % Calculate Individual Thrusts: T = Ct * w^2
    T = P.Ct .* (omega_rot.^2);
    
    f1 = [0; 0; T(1)];
    f2 = [0; 0; T(2)];
    f3 = [0; 0; T(3)];
    f4 = [0; 0; T(4)];
    
    % Calculate Individual Torques: Tau = Cq * w^2
    % (Note: Motors 2 and 4 spin clockwise, 1 and 3 counter-clockwise)
    n1 = [0; 0;  P.Cq * omega_rot(1)^2];
    n2 = [0; 0; -P.Cq * omega_rot(2)^2];
    n3 = [0; 0;  P.Cq * omega_rot(3)^2];
    n4 = [0; 0; -P.Cq * omega_rot(4)^2];
    
    % Total Propulsion Force
    fp = f1 + f2 + f3 + f4;
    
    % Total Propulsion Moment
    np = (n1 + cross(p1,f1)) + (n2 + cross(p2,f2)) + ...
         (n3 + cross(p3,f3)) + (n4 + cross(p4,f4));
end


%% =========================================================
% HELPER FUNCTIONS
%% =========================================================
function R = rotm_zyx(lambda)
    phi = lambda(1); theta = lambda(2); psi = lambda(3);
    Rx = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
    Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    Rz = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
    R = Rz*Ry*Rx;
end



function Q = euler_rate_matrix(lambda)
    phi = lambda(1); theta = lambda(2);
    Q = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
         0, cos(phi),           -sin(phi);
         0, sin(phi)/cos(theta), cos(phi)/cos(theta)];
end


function S = skew(a)
    S = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
end



function a = wrap180(a)
    for i = 1:numel(a)
        while a(i) > 180
            a(i) = a(i) - 360;
        end
        while a(i) < -180
            a(i) = a(i) + 360;
        end
    end
end