%% ================================================
%  1.8, 1.9, 1.10, 1.11
%% ================================================
clear; clc; close all;

%% ---- PARAMETERS ----
m   = 0.027;
g   = 9.81;
Ixx = 1.4730e-5;
Iyy = 1.4797e-5;
Izz = 2.8476e-5;
rho = 1.225;
Cd  = 0.3;
S   = 0.0081;
kd  = 0.5*rho*Cd*S;
ve  = 10;

Z  = zeros(3,3);
I3 = eye(3);

%% ================================================
%% 1.8 — EIGENVALUES
%% ================================================

%% ---- OP.1: HOVER ----
grav_hover = [0, -g,  0;
              g,  0,  0;
              0,  0,  0];

A_hover = [Z,  I3,         Z,          Z;
           Z,  Z,          grav_hover, Z;
           Z,  Z,          Z,          I3;
           Z,  Z,          Z,          Z];

eig_hover = eig(A_hover);

%% ---- OP.2: FORWARD FLIGHT ----
theta_e = atan((kd*ve^2)/(m*g));
ct = cos(theta_e);
st = sin(theta_e);

fprintf('Forward flight equilibrium:\n');
fprintf('  theta_e = %.4f rad = %.2f deg\n', theta_e, rad2deg(theta_e));
fprintf('  T_e     = %.4f N\n', sqrt((kd*ve^2)^2 + (m*g)^2));

% R(lambda_e) — pure pitch
R_e = [ct,  0,  st;
        0,  1,   0;
      -st,  0,  ct];

% dR/dlambda * ve
dRdtheta_mat = [-st,  0,  ct;
                  0,  0,   0;
                -ct,  0, -st];
v_e = [ve; 0; 0];
dRdlambda_ve = [zeros(3,1), dRdtheta_mat*v_e, zeros(3,1)];

% gravity in ENU
gI = [0; 0; -g];

% 1/m * dR^T/dlambda * gI
dRTdphi = [0,            0,           0;
           0,            0,  cos(theta_e);
           0, -cos(theta_e),           0];

dRTdtheta = [-st,  0, -ct;
               0,  0,   0;
              ct,  0, -st];

dRTdpsi = [0, -1,  0;
           1,  0,  0;
           0,  0,  0];

grav_fwd = (1/m) * [dRTdphi*gI, dRTdtheta*gI, dRTdpsi*gI];

% -S(ve)
mSve = [0,    0,    0;
        0,    0,   ve;
        0,  -ve,    0];

% drag linearization
drag = -(2*kd*ve/m)*I3;

% Q(lambda_e)
Q_e = [1,  0,  st;
       0,  1,   0;
       0,  0,  ct];

A_fwd = [Z,  R_e,   dRdlambda_ve,  Z;
         Z,  drag,  grav_fwd,      mSve;
         Z,  Z,     Z,             Q_e;
         Z,  Z,     Z,             Z];

eig_fwd = eig(A_fwd);

%% ---- EIGENVALUE PLOT ----
figure('Name','1.8 — Eigenvalues');
hold on; grid on;
plot(real(eig_hover), imag(eig_hover), 'bx', ...
    'MarkerSize', 14, 'LineWidth', 2.5, ...
    'DisplayName', 'OP.1 Hover');
plot(real(eig_fwd), imag(eig_fwd), 'ro', ...
    'MarkerSize', 10, 'LineWidth', 2.5, ...
    'DisplayName', 'OP.2 Forward Flight');
xline(0, '--k', 'LineWidth', 1.2);
yline(0, '--k', 'LineWidth', 1.2);
xlabel('Real part', 'FontSize', 13);
ylabel('Imaginary part', 'FontSize', 13);
title('1.8 — Eigenvalues: OP.1 Hover vs OP.2 Forward Flight', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 12);
axis tight;

fprintf('\n--- OP.1 Hover eigenvalues ---\n');      disp(eig_hover);
fprintf('\n--- OP.2 Forward Flight eigenvalues ---\n'); disp(eig_fwd);






%% ================================================
%% 1.9 — CONTROLLABILITY, OBSERVABILITY, STABILITY
%% ================================================

% B matrix (same for both OPs)
B = [zeros(3,1),    zeros(3,3);
     (1/m)*[0;0;1], zeros(3,3);
     zeros(3,1),    zeros(3,3);
     zeros(3,1),    diag([1/Ixx, 1/Iyy, 1/Izz])];

% Output matrix — full state
Cmat = eye(12);

fprintf('\n================================================\n');
fprintf('1.9 — CONTROLLABILITY, OBSERVABILITY, STABILITY\n');
fprintf('================================================\n');

for op = 1:2
    if op == 1
        A_op = A_hover;
        ev   = eig_hover;
        name = 'OP.1 HOVER';
    else
        A_op = A_fwd;
        ev   = eig_fwd;
        name = 'OP.2 FORWARD FLIGHT';
    end

    fprintf('\n--- %s ---\n', name);

    % Controllability
    Mc      = ctrb(A_op, B);
    rk_Mc   = rank(Mc);
    fprintf('Controllability rank: %d / 12 --> %s\n', rk_Mc, ...
        ternary(rk_Mc==12, 'CONTROLLABLE', 'NOT controllable'));

    % Observability
    Mo      = obsv(A_op, Cmat);
    rk_Mo   = rank(Mo);
    fprintf('Observability rank:   %d / 12 --> %s\n', rk_Mo, ...
        ternary(rk_Mo==12, 'OBSERVABLE', 'NOT observable'));

    % Stability
    max_re = max(real(ev));
    if max_re < -1e-10
        stab = 'EXPONENTIALLY STABLE';
    elseif max_re <= 1e-10
        [V, ~] = eig(A_op);
        if rank(V) < 12
            stab = 'MARGINALLY UNSTABLE';
        else
            stab = 'MARGINALLY STABLE';
        end
    else
        stab = 'UNSTABLE (positive real eigenvalue)';
    end
    fprintf('Max Re(lambda): %.4f --> %s\n', max_re, stab);
end






%% ================================================
%% 1.10 — TRANSFER FUNCTIONS
%% ================================================

fprintf('\n================================================\n');
fprintf('1.10 — TRANSFER FUNCTIONS\n');
fprintf('================================================\n');

% Inner loops
G_nx_phi   = tf(1, [Ixx 0 0]);
G_ny_theta = tf(1, [Iyy 0 0]);
G_nz_psi   = tf(1, [Izz 0 0]);

% Outer loops
G_phi_py   = tf(g,  [1 0 0]);
G_theta_px = tf(-g, [1 0 0]);
G_T_pz     = tf(1,  [m 0 0]);

fprintf('\nG_nx_phi   = phi(s)/nx(s):\n');    zpk(G_nx_phi)
fprintf('G_ny_theta = theta(s)/ny(s):\n');    zpk(G_ny_theta)
fprintf('G_nz_psi   = psi(s)/nz(s):\n');      zpk(G_nz_psi)
fprintf('G_phi_py   = py(s)/phi(s):\n');      zpk(G_phi_py)
fprintf('G_theta_px = px(s)/theta(s):\n');    zpk(G_theta_px)
fprintf('G_T_pz     = pz(s)/T(s):\n');        zpk(G_T_pz)

% Bode plots
figure('Name','1.10 — Bode Plots');
subplot(3,2,1); bode(G_nx_phi);   title('G_{n_x,\phi}(s)');   grid on;
subplot(3,2,2); bode(G_phi_py);   title('G_{\phi,p_y}(s)');   grid on;
subplot(3,2,3); bode(G_ny_theta); title('G_{n_y,\theta}(s)'); grid on;
subplot(3,2,4); bode(G_theta_px); title('G_{\theta,p_x}(s)'); grid on;
subplot(3,2,5); bode(G_nz_psi);   title('G_{n_z,\psi}(s)');   grid on;
subplot(3,2,6); bode(G_T_pz);     title('G_{T,p_z}(s)');      grid on;
sgtitle('1.10 — Bode Plots: Inner and Outer Transfer Functions');





%% ================================================
%% 1.11 — ROOT LOCUS / CLOSED-LOOP STABILITY
%% ================================================

figure('Name','1.11 — Root Locus');

subplot(1,2,1);
rlocus(G_T_pz);
title('Root Locus — G_{T,p_z}(s) = 1/(ms^2)', 'FontSize', 12);
xline(0, '--k', 'LineWidth', 1);
grid on;

subplot(1,2,2);
rlocus(G_phi_py);
title('Root Locus — G_{\phi,p_y}(s) = g/s^2', 'FontSize', 12);
xline(0, '--k', 'LineWidth', 1);
grid on;


%% ================================================
%% HELPER FUNCTION
%% ================================================
function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end