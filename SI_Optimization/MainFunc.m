%==========================================================================
% Main Function
% Author: Tina Gui
% Developed in MATLAB R2016a
%
% Objective functions:
% F1 - Sphere
% F2 - Step
% F3 - Schwefel
% F4 - Griewank
% F5 - Cigar
% F6 - APHE (Axis Parallel Hyper-Ellipsoid)
%==========================================================================

function MainFunc(swarm_no,max_iter,dim,f)
% swarm_no - number of search agents
% max_iter - maximum number of generations
% dim - number of your variables
% lb_no - lower bound value
% ub_no - upper bound value
% CostFunction - objective function

if (nargin < 1)
    swarm_no = 30;
    max_iter = 50;
    dim = 30;
    f = 'F1';
end

switch(f)
    case 'F1'
        lb_no = -5.12;
        ub_no = 5.12;
        CostFunction = @Sphere; % Objective function
        title('F1 - Sphere');
    case 'F2'
        lb_no = -100;
        ub_no = 100;
        CostFunction = @Sphere; % Objective function
        title('F2 - Step');
    case 'F3'
        lb_no = -500;
        ub_no = 500;
        CostFunction = @Sphere; % Objective function
        title('F3 - Schwefel');
    case 'F4'
        lb_no = -600;
        ub_no = 600;
        CostFunction = @Sphere; % Objective function
        title('F4 - Griewank');
    case 'F5'
        lb_no = -10;
        ub_no = 10;
        CostFunction = @Sphere; % Objective function
        title('F5 - Cigar');
    case 'F6'
        lb_no = -5.12;
        ub_no = 5.12;
        CostFunction = @APHE; % Objective function
        title('F6 - Axis Parallel Hyper-Ellipsoid');
    otherwise
        fprintf('Invalid objective function\n' );
end
GWO2(swarm_no,max_iter,lb_no,ub_no,dim,CostFunction);
WOA2(swarm_no,max_iter,lb_no,ub_no,dim,CostFunction);
SSO2(swarm_no,max_iter,lb_no,ub_no,dim,CostFunction);
MOA2(swarm_no,max_iter,lb_no,ub_no,dim,CostFunction);
hold on;

% Display legend (optional)
h = zeros(4, 1);
h(1) = plot(NaN,NaN,':k','LineWidth', 2); % GWO
h(2) = plot(NaN,NaN,'--b','LineWidth', 1.5); % WOA
h(3) = plot(NaN,NaN,'-.r','LineWidth', 1.5);  % SSO
h(4) = plot(NaN,NaN,'-m','LineWidth', 1.5);  % MOA
legend(h, 'GWO','WOA','SSO', 'MOA');