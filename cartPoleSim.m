clear
clc

%% Simulation settings
% Initial conditions
sim.q0 = [-0.3; 4.0];
sim.p0 = [0; 0];

% Simulation time
sim.Tend = 15;

% Option to generate level-set plots (slow)
generateLevelSets = 0;

% Code to verify Hd at initial condition (use after running simulation)
% nonlin.Hd(sim.q0,sim.p0)

%% Nonlinear model
% Define degrees of freedom
sys.N = length(sim.q0);

% symbolic variables
syms q_sym [sys.N,1]
syms p_sym [sys.N,1]

% model parameters
sys.l = 1;
sys.mp = 1;
sys.mc = 1;
sys.g = 1;
sys.d1 = 1;
sys.d2 = 1;

% Mass matrix
sys.m11 = sys.mp*sys.l^2;
sys.m12 = @(q) sys.mp*sys.l*cos(q(1));
sys.m22 = sys.mc + sys.mp;
sys.M = @(q) [sys.m11, sys.m12(q); sys.m12(q), sys.m22];

% Damping matrix
sys.D = [sys.d1, 0; 0, sys.d2];

% Input matrix
sys.G = [0; 1];
sys.Gp = [1, 0];

% Potential energy
sys.k = sys.mp*sys.g*sys.l;
sys.V = @(q) sys.k*cos(q(1));

% Kinetic energy
sys.T = @(q,p) 0.5*p.'*(sys.M(q)\p);

% System energy
sys.H = @(q,p) sys.T(q,p) + sys.V(q);

% energy gradients
sys.dVdq = matlabFunction(jacobian(sys.V(q_sym),q_sym).','vars',{q_sym});
sys.dTdq = matlabFunction(jacobian(sys.T(q_sym,p_sym),q_sym).','vars',[{q_sym}, {p_sym}]);
sys.dHdq = matlabFunction(jacobian(sys.H(q_sym,p_sym),q_sym).','vars',[{q_sym}, {p_sym}]);
sys.dHdp = @(q,p) sys.M(q)\p;

% Compute the C matrix
sys.dTdqdp = matlabFunction(jacobian(sys.dTdq(q_sym,p_sym),p_sym),'vars',[{q_sym}, {p_sym}]);
sys.C = @(q,p) 0.5*sys.dTdqdp(q,p)*sys.M(q);
% Test the definition of C -- should equal zero
% simplify(sys.C(q_sym,p_sym)*(sys.M(q_sym)\p_sym) - sys.dTdq(q_sym,p_sym))

% Interconnection and damping structure
sys.F = [zeros(sys.N) eye(sys.N); -eye(sys.N) -sys.D];

% system ODE
sys.dx = @(q,p,u) sys.F*[sys.dHdq(q,p); sys.dHdp(q,p)] + [zeros(sys.N,1); sys.G]*u;

%% Linearised system and controller
% Define linearised system
lin.M = sys.M(zeros(2,1));
lin.D = sys.D;
lin.K = double(subs(jacobian(sys.dVdq(q_sym),q_sym), q_sym, zeros(2,1)));
lin.A = [zeros(sys.N) eye(sys.N)/lin.M; -lin.K -sys.D/lin.M];
lin.B = [zeros(sys.N,1); sys.G];

% Construct a linear control law
[lin.K_ctrl,~,~] = lqr(lin.A,lin.B,diag([1 1 1 1]),1);
lin.ctrl = @(x) -lin.K_ctrl*x;

% Construct the linear closed-loop dynamics
lin.Ad = lin.A - lin.B*lin.K_ctrl;

% Select a damping structure for the linear system
lin.Rd = diag([1 1 1 1]);

% Solve the Riccatti equation to find the closed-loop Hessian
lin.Qdi = lyap(lin.Ad, 2*lin.Rd);
lin.Qd  = inv(lin.Qdi);

% Construct the closed-loop interconncetion and damping matrix for teh
% linear system
lin.Fd = lin.Ad*lin.Qdi;

% Verify controller definition
%K = -(B.'*B)\B.'*(Fd*Qd - A)

% Break Qd into components
lin.Qd11 = lin.Qd(1:2,1:2);
lin.Qd12 = lin.Qd(1:2,3:4);
lin.Qd22 = lin.Qd(3:4,3:4);

% Construct the linear energy function terms
lin.Mdi = lin.Qd22;
lin.Md = inv(lin.Mdi);
lin.Ja = lin.Qd22\lin.Qd12.';
lin.Kd = lin.Qd11 - lin.Qd12*lin.Md*lin.Qd12.';
%lin.Fd11 = -lin.M\(lin.Ja/lin.Kd);

% Break Fd into components. This vaue for Fd11 should agree with above
lin.Fd11 = lin.Fd(1:2,1:2);
lin.Fd12 = lin.Fd(1:2,3:4);
lin.Fd21 = lin.Fd(3:4,1:2);
lin.Fd22 = lin.Fd(3:4,3:4);

%% Construct nonlinear controller
% % Choose Fd11
nonlin.kd = diag([11, 11]);
nonlin.Fd11 = @(q) lin.Fd11 - (1-cos(q(1)))*nonlin.kd;

% Select closed-loop mass
nonlin.Mdi = @(q) lin.Mdi;

% Select Vd to simplify matching condition
nonlin.Vd = @(q) 0.5*q.'*lin.Kd*q;
nonlin.dVddq = @(q) lin.Kd*q;
Vhess0 = double(subs(hessian(sys.V(q_sym), q_sym), ...
                     q_sym, [0; 0]));
nonlin.k = -Vhess0(1,1);
nonlin.SV = @(q) [-nonlin.k*localSinc(q(1)/pi), 0];
nonlin.Phi1_2 = @(q) -nonlin.SV(q)/lin.Kd;

% Construct the vactor alpha and jacobian
nonlin.alpha = @(q) -sys.M(q)*nonlin.Fd11(q)*nonlin.dVddq(q);
nonlin.dalphadq = matlabFunction(jacobian(nonlin.alpha(q_sym),q_sym),'vars',{q_sym});

% Compute the Cd matrix
nonlin.bTd = @(q,p) 0.5*p.'*(nonlin.Mdi(q)*p);
nonlin.dbTddq = matlabFunction(jacobian(nonlin.bTd(q_sym,p_sym),q_sym).','vars',[{q_sym}, {p_sym}]);
nonlin.dbTddqdp = matlabFunction(jacobian(nonlin.dbTddq(q_sym,p_sym),p_sym),'vars',[{q_sym}, {p_sym}]);
nonlin.Cd1 = @(q,p) 0.5*nonlin.dbTddqdp(q,p)/nonlin.Mdi(q);
nonlin.Cd0 = @(q) 0.5*nonlin.dbTddqdp(q,nonlin.alpha(q))/nonlin.Mdi(q);
nonlin.Cd = @(q,p) nonlin.Cd0(q) + nonlin.Cd1(q,p);

% Construct Fd12
nonlin.Fd12 = @(q,p) (sys.M(q)\eye(2))/nonlin.Mdi(q) - nonlin.Fd11(q)*(nonlin.dalphadq(q).' + nonlin.Cd(q,p));

% Construct T matrix
nonlin.T = @(q) [sys.Gp; sys.G.'];
nonlin.Ti = @(q) [sys.Gp.'/(sys.Gp*sys.Gp.'), sys.G/(sys.Gp*sys.Gp.')];

% Choose Phi1
nonlin.Phi1_1 = @(q,p) -sys.Gp*(sys.C(q,p) + sys.D)*nonlin.Fd11(q);
% nonlin.Phi1_2 = @(q) -(1/nonlin.beta)*[1, -nonlin.Lv(1,1)/nonlin.Lv(1,2)];
nonlin.Phi1 = @(q,p) nonlin.Phi1_1(q,p) + nonlin.Phi1_2(q);

% % Verify matching condition. Should be equal to 0
% nonlin.Phi1_2(q_sym)*nonlin.dVddq(q_sym) + sys.Gp*sys.dVdq(q_sym)

% Choose Phi 2. Take it to make free component of Fd21 skew symmetric with
% Fd12
nonlin.Phi2 = @(q,p) [0, 1]*nonlin.T(q)*(-nonlin.Fd12(q,p).');

% Construct Phi matrix
nonlin.Phi = @(q,p) [nonlin.Phi1(q,p); nonlin.Phi2(q,p)];

% Construct Fd21
nonlin.Fd21 = @(q,p) nonlin.Ti(q)*nonlin.Phi(q,p);

% Set Gamma
% nonlin.T([0;0])*lin.Fd22*nonlin.T([0;0]).'
nonlin.Gamma11_12 = @(q,p) -sys.Gp*(((sys.C(q,p) + sys.D)/sys.M(q))/nonlin.Mdi(q) + nonlin.Fd21(q,p)*(nonlin.Cd(q,p) + nonlin.dalphadq(q).'))*nonlin.T(q).';
nonlin.Gamma21_22 = @(q,p) [-nonlin.Gamma11_12(q,p)*[0; 1], -1];
nonlin.Gamma = @(q,p) [nonlin.Gamma11_12(q,p);
                        nonlin.Gamma21_22(q,p)];
nonlin.Gamma11 = @(q,p) nonlin.Gamma11_12(q,p)*[1; 0];

% Construct Fd22
nonlin.Fd22 = @(q,p) nonlin.Ti(q)*nonlin.Gamma(q,p)*nonlin.Ti(q).';

% Construct the full Fd matrix
nonlin.Fd = @(q,p) [nonlin.Fd11(q), nonlin.Fd12(q,p);
                    nonlin.Fd21(q,p), nonlin.Fd22(q,p)];

% Define closed-loop energy function and gradient
nonlin.Hd = @(q,p) 0.5*(p + nonlin.alpha(q)).'*(nonlin.Mdi(q)*(p + nonlin.alpha(q))) + nonlin.Vd(q);
nonlin.dHdp = matlabFunction(jacobian(nonlin.Hd(q_sym,p_sym),p_sym).','vars',[{q_sym}, {p_sym}]);
nonlin.dHdq = matlabFunction(jacobian(nonlin.Hd(q_sym,p_sym),q_sym).','vars',[{q_sym}, {p_sym}]);

% COnstruct the nonlinear control law
nonlin.g = [0;0;sys.G];
u = @(q,p) ((nonlin.g.'*nonlin.g)\nonlin.g.')*(nonlin.Fd(q,p)*[nonlin.dHdq(q,p); nonlin.dHdp(q,p)] - sys.F*[sys.dHdq(q,p); sys.dHdp(q,p)]);

% Construct Lambda matrix
nonlin.Lambda11 = @(q,p) [nonlin.Fd11(q), nonlin.Fd12(q,p)*sys.Gp.';
                nonlin.Phi1(q,p), nonlin.Gamma11(q,p)];

% Test the negativity of Fd
q_test = [-0.6; 0.0];
p_test = [0; 0];
eig(nonlin.Fd(q_test,p_test) + nonlin.Fd(q_test,p_test).')
eig(nonlin.Lambda11(q_test,p_test) + nonlin.Lambda11(q_test,p_test).')

%% Run simulation
sim.u_wrap = @(x) u(x(1:2),x(3:4));
sim.sys_wrap = @(t,x) sys.dx(x(1:2),x(3:4),sim.u_wrap(x));
% sim.sys_wrap = @(t,x) sys.dx(x(1:2),x(3:4),lin.ctrl(x));
sim.tspan = [0 sim.Tend];
sim.x0 = [sim.q0; sim.p0];

sim.options = odeset('RelTol',1e-6);
[res.t, res.x] = ode45(sim.sys_wrap, sim.tspan, sim.x0, sim.options);

%% Plot results
% Seperate q and p
res.q = res.x(:,1:2);
res.p = res.x(:,3:4);

% Evaluate the Hamiltonian functions and verify the damping is negative
% along trajectory
res.H = zeros(length(res.t), 1);
res.Hd = zeros(length(res.t), 1);
res.lambda_max = zeros(length(res.t), 1);
for i=1:length(res.t)
    res.H(i) = sys.H(res.q(i,:).',res.p(i,:).');
    res.Hd(i) = nonlin.Hd(res.q(i,:).',res.p(i,:).');

    % Maximum eigenvalue of symmetric Lambda11
    Lambda = nonlin.Lambda11(res.q(i,:).',res.p(i,:).');
    res.lambda_max(i) = max(eig(Lambda + Lambda.'));
end
% Check that the damping was positive along the trajectory
% plot(res.lambda_max)

% Plot state trajectories and closed-loop hamiltonian
figure(1)
subplot(2,1,1)
plot(res.t, res.q,'LineWidth',1.5)
legend('q_1','q_2')
grid on
subplot(2,1,2)
plot(res.t, res.p,'LineWidth',1.5)
legend('p_1','p_2')
grid on
figure(2)
plot(res.t, log(res.Hd),'LineWidth',1.5)
legend('log(H_d)')
grid on

%% Generate configuration level set plot
if(generateLevelSets)
    % Configuration grid
    q1_vals = linspace(-0.7, 0.7, 200);
    q2_vals = linspace(-10, 10, 200);
    HamiltonianLimit = 0.112;
    
    % Evaluate the hamiltonian and eigenvales of Lambda11
    lambda_max = zeros(length(q1_vals), length(q2_vals));
    Hd_val = zeros(length(q1_vals), length(q2_vals));
    for i = 1:length(q1_vals)
        for j = 1:length(q2_vals)
            % Test points
            q_test = [q1_vals(i); q2_vals(j)];
            p_test = [0; 0];
    
            % Maximum eigenvalue of symmetric Lambda11
            Lambda = nonlin.Lambda11(q_test,p_test);
            lambda_max(i,j) = max(eig(Lambda + Lambda.'));
    
            % Desired Hamiltonian
            Hd_val(i,j) = nonlin.Hd(q_test,p_test);
        end
    end
    
    % Define regions
    % 0 = negative semidefinite
    % 1 = positive eigenvalue
    region = double(lambda_max > 0);
    
    % Mask Hamiltonian outside dissipative region
    % and above the selected energy level
    Hd_plot = Hd_val;
    Hd_plot(lambda_max > 0) = NaN;
    Hd_plot(Hd_plot > HamiltonianLimit) = NaN;
    
    % Create meshgrid
    [Q1,Q2] = meshgrid(q1_vals,q2_vals);
    
    % Create figure
    figure

    % AXIS 1: White/grey dissipativity region
    ax1 = axes;
    ax1.FontSize = 16;
    hold(ax1,'on');
    contourf(ax1,Q1,Q2,region.', [-0.5 0.5 1.5], 'LineColor','none');
    colormap(ax1,[1 1 1; 0.8 0.8 0.8]);
    
    % Boundary where Lambda11 loses negative semidefiniteness
    contour(ax1,Q1,Q2,lambda_max.',[0 0], 'k','LineWidth',1.5);
    
    % Labels
    xlabel(ax1,'$q_1$','Interpreter','latex','FontSize',20);
    ylabel(ax1,'$q_2$','Interpreter','latex','FontSize',20);

    % AXIS 2: Coloured Hamiltonian sublevel sets
    ax2 = axes;
    ax2.FontSize = 16;
    hold(ax2,'on');
    
    % Make background transparent
    ax2.Color = 'none';
    
    % Coloured Hamiltonian contours
    contourf(ax2,Q1,Q2,Hd_plot.', linspace(0,HamiltonianLimit,20), 'LineColor','none');
    colormap(ax2,turbo);
    clim(ax2,[0 0.1]);
    contour(ax2,Q1,Q2,Hd_val.',[HamiltonianLimit HamiltonianLimit], 'r','LineWidth',1.5);
    
    % Axes
    xlabel(ax2,'$q_1$','Interpreter','latex','FontSize',20);
    ylabel(ax2,'$q_2$','Interpreter','latex','FontSize',20);
end

%% Generate momentum level set plot
if(generateLevelSets)
    % Momentum grid
    p1_vals = linspace(-1, 1, 200);
    p2_vals = linspace(-2, 2, 200);
    % HamiltonianLimit = 0.07;
    
    % Evaluate the hamiltonian and eigenvales of Lambda11
    lambda_max = zeros(length(p1_vals), length(p2_vals));
    Hd_val = zeros(length(p1_vals), length(p2_vals));
    for i = 1:length(p1_vals)
        for j = 1:length(p2_vals)
            % Test points
            q_test = [0; 0];
            p_test = [p1_vals(i); p2_vals(j)];
    
            % Maximum eigenvalue of symmetric Lambda11
            Lambda = nonlin.Lambda11(q_test,p_test);
            lambda_max(i,j) = max(eig(Lambda + Lambda.'));
    
            % Desired Hamiltonian
            Hd_val(i,j) = nonlin.Hd(q_test,p_test);
        end
    end
    
    % Define regions
    % 0 = negative semidefinite
    % 1 = positive eigenvalue
    region = double(lambda_max > 0);
    
    % Mask Hamiltonian outside dissipative region
    % and above the selected energy level
    Hd_plot = Hd_val;
    
    Hd_plot(lambda_max > 0) = NaN;
    Hd_plot(Hd_plot > HamiltonianLimit) = NaN;
    
    % Create meshgrid
    [Q1,Q2] = meshgrid(p1_vals,p2_vals);
    
    % Create figure
    figure
    
    % AXIS 1: White/grey dissipativity region
    ax1 = axes;
    ax1.FontSize = 16;
    hold(ax1,'on');
    xlabel(ax1,'$p_1$','Interpreter','latex','FontSize',20);
    ylabel(ax1,'$p_2$','Interpreter','latex','FontSize',20);
    contourf(ax1,Q1,Q2,region.', [-0.5 0.5 1.5], 'LineColor','none');
    colormap(ax1,[1 1 1; 0.8 0.8 0.8]);
    
    % Boundary where Lambda11 loses negative semidefiniteness
    contour(ax1,Q1,Q2,lambda_max.',[0 0], 'k','LineWidth',1.5);

    % AXIS 2: Coloured Hamiltonian sublevel sets
    ax2 = axes;
    ax2.FontSize = 16;
    hold(ax2,'on');
    % Make background transparent
    ax2.Color = 'none';
    
    % Coloured Hamiltonian contours
    contourf(ax2,Q1,Q2,Hd_plot.', linspace(0,HamiltonianLimit,20), 'LineColor','none');
    colormap(ax2,turbo);
    clim(ax2,[0 0.1]);
    contour(ax2,Q1,Q2,Hd_val.',[HamiltonianLimit HamiltonianLimit], 'r','LineWidth',1.5);
    xlabel(ax2,'$p_1$','Interpreter','latex','FontSize',20);
    ylabel(ax2,'$p_2$','Interpreter','latex','FontSize',20);
end

%% Functions
function y = localSinc(x)
    if abs(x) < 1e-6
        % Taylor expansion near zero for numerical accuracy
        y = 1 - x^2/6 + x^4/120;
    else
        y = sin(x)/x;
    end
end