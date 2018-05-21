%==========================================================================
% The Monkey Optimization Algorithm
% Author: Tina Gui (tgui@go.olemiss.edu)
%
% Developed in MATLAB R2016a
%
% Reference: Bansal et al.
%            Spider monkey optimization algorithm for numerical optimization 
%            Memetic computing 6, no. 1 (2014): 31-47.
%
% 
%==========================================================================

function MOA2(swarm_no,max_iter,lb_no,ub_no,dim,CostFunction)
% CostFunction - @YourCostFunction
% dim - number of your variables
% max_iter - maximum number of generations
% swarm_no - number of search agents
% lb_no - lower bound value
% ub_no - upper bound value


%% Parameters setting
ub = ub_no * ones(1, dim); % upper bound
lb = lb_no * ones(1, dim); % lower bound
pr = 0.5; % perturbation rate 
iter = 0;% iteration counter


%% Population initialization
Local_pos = zeros(1, dim);
Local_score = inf; 
Positions = initialization(swarm_no, dim, ub, lb);
Convergence_curve = zeros(1, max_iter);


%% Main Loop
while iter < max_iter
    
    for i = 1:size(Positions,1)
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub = Positions(i,:)>ub;
        Flag4lb = Positions(i,:)<lb;
        Positions(i,:) = (Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate objective function for each search agent
        fitness = CostFunction(Positions(i,:), dim);
        All_fitness(1,i) = fitness;
        
        % Update the leader
        if fitness < Local_score % Change this to > for maximization problem
            Local_score = fitness; % Update alpha
            Local_pos = Positions(i,:);
        end
    end

    % Local Leader Phase 
    for i = 1:size(Positions,1)
        p = rand();
        r1 = rand(); % r1 is a random number in [0,1]
        a = -1; b = 1;
        r2 = a + (b-a).*rand(); % r2 is a random number in [-1,1]
        
        for j = 1:size(Positions,2) % j in each dimension
            if p > pr % p greater than perturbation rate 
                SM_min = min(Positions(:,j));
                SM_max = min(Positions(:,j));
                Positions(i,j) = SM_min + rand() * (SM_max - SM_min);
                
            elseif p <= pr % p less than and equal to perturbation rate 
                LLkj = r1 * (Local_pos(j) - Positions(i,j));
                r_index = floor(swarm_no*rand() + 1);
                SMr = Positions(r_index, :);
                SMrj = r2 * (SMr(j) - Positions(i,j));
                % SMnewij = SMij + R(0,1)*(LLkj-SMij) + R(-1,1)*(SMrj-SMij;
                % r != i
                Positions(i,j) = Positions(i,j) + LLkj + SMrj;       
            end
        end
    end
    
    iter = iter+1;
    Convergence_curve(iter) = Local_score;

    if iter > 2
        line([iter-1 iter], [Convergence_curve(iter-1) Convergence_curve(iter)], 'Color','Magenta','LineStyle','-', 'LineWidth', 1.5);
%         title('Monkey Optimization Algorithm');
        xlabel('Number of Iteration');
        ylabel('Best Score Obtained');        
        drawnow
    end
 
    results{1,1} = iter;
    results{1,2} = Local_score;

end




