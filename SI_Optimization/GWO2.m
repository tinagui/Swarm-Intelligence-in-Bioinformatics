%==========================================================================
% The Grey Wolf Optimization Alogrithm
% Author: Tina Gui (tgui@go.olemiss.edu)
% Modified on Seyedali Mirjalili's original code
% (seyedali.mirjalili@griffithuni.edu.au)
%
% Developed in MATLAB R2016a
%
% Reference: http://www.alimirjalili.com                                   
%            Grey Wolf Optimizer, Advances in Engineering        
%            Software, Volume 69, March 2014, Pages 46-61,       
%            http://dx.doi.org/10.1016/j.advengsoft.2013.12.007
%
%==========================================================================

function GWO2(swarm_no,max_iter,lb,ub,dim,CostFunction)
% CostFunction - @YourCostFunction
% dim - number of your variables
% max_iter - maximum number of generations
% swarm_no - number of search agents
% lb_no - lower bound value
% ub_no - upper bound value

%% Parameters setting
ub = 5.12 * ones(1, dim); % upper bound
lb = -5.12 * ones(1, dim); % lower bound

Alpha_pos=zeros(1,dim);
Alpha_score= inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score= inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score= inf; %change this to -inf for maximization problems

%Initialize the positions of search agents
Positions=initialization(swarm_no,dim,ub,lb);

%Convergence_curve=zeros(1,max_iter);

l=0;% Loop counter

%% Main Loop
while l<max_iter
    
    for i=1:size(Positions,1)             
        % Calculate objective function for each search agent
        fitness=CostFunction(Positions(i,:), dim);
        All_fitness(1,i)=fitness;
        
        % Update Alpha, Beta, and Delta
        if fitness<Alpha_score 
            Alpha_score=fitness; % Update alpha
            Alpha_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness<Beta_score 
            Beta_score=fitness; % Update beta
            Beta_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score 
            Delta_score=fitness; % Update delta
            Delta_pos=Positions(i,:);
        end
    end
    
    a=2-l*((2)/max_iter); % a decreases linearly fron 2 to 0
    
    % Update the Position of search agents including omegas
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)     
                       
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % Equation (3.3)
            C1=2*r2; % Equation (3.4)
            
            D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % Equation (3.5)-part 1
            X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
                       
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2; % Equation (3.4)
            
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (3.5)-part 2
            X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
            
            r1=rand();
            r2=rand(); 
            
            A3=2*a*r1-a; % Equation (3.3)
            C3=2*r2; % Equation (3.4)
            
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
            X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
            
            Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)
            
        end
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    end
    l=l+1;    
    Convergence_curve(l)=Alpha_score;
    if l>1
        line([l-1 l], [Convergence_curve(l-1) Convergence_curve(l)],'Color','Black','LineStyle',':','LineWidth', 2)
%         title('Grey Wolf Optimization');
        xlabel('Number of Iteration');
        ylabel('Best Score Obatained');        
        drawnow
    end

    results{1,1} = l;
    results{1,2} = Alpha_score;

end



