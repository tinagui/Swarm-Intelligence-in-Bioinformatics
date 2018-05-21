%==========================================================================
% The Whale Optimization Algorithmm
% Author: Tina Gui (tgui@go.olemiss.edu)
% Modified on Seyedali Mirjalili's original code
% (seyedali.mirjalili@griffithuni.edu.au)
%
% Developed in MATLAB R2016a
%
% Reference: http://www.alimirjalili.com                                   
%            The Whale Optimization Algorithm,                         
%            Advances in Engineering Software , in press, 2016         
%            DOI: http://dx.doi.org/10.1016/j.advengsoft.2016.01.008  
%
%==========================================================================

function WOA2(swarm_no,max_iter,lb_no,ub_no,dim,CostFunction)
% CostFunction - @YourCostFunction
% dim - number of your variables
% max_iter - maximum number of generations
% swarm_no - number of search agents
% lb_no - lower bound value
% ub_no - upper bound value


%% Parameters setting
ub = ub_no * ones(1, dim); % upper bound
lb = lb_no * ones(1, dim); % lower bound

% initialize position vector and score for the leader
Leader_pos=zeros(1,dim);
Leader_score= inf; %change this to -inf for maximization problems

%Initialize the positions of search agents
Positions=initialization(swarm_no,dim,ub,lb);

Convergence_curve=zeros(1,max_iter);

t=0;% Loop counter

%% Main Loop
while t<max_iter
    
    for i=1:size(Positions,1) 
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate objective function for each search agent
        fitness=CostFunction(Positions(i,:), dim);
        All_fitness(1,i)=fitness;
        
        % Update the leader
        if fitness<Leader_score % Change this to > for maximization problem
            Leader_score=fitness; % Update alpha
            Leader_pos=Positions(i,:);
        end  
    end
    
    a=2-t*((2)/max_iter); % a decreases linearly fron 2 to 0 in Eq. (2.3)
    
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+t*((-1)/max_iter);
    
    % Update the Position of search agents 
    for i=1:size(Positions,1)
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper
        
        
        b=1;               %  parameters in Eq. (2.5)
        l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
        
        p = rand();        % p in Eq. (2.6)
        
        for j=1:size(Positions,2)
            
            if p<0.5   
                if abs(A)>=1
                    rand_leader_index = floor(swarm_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand=abs(C*X_rand(j)-Positions(i,j)); % Eq. (2.7)
                    Positions(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
                    
                elseif abs(A)<1
                    D_Leader=abs(C*Leader_pos(j)-Positions(i,j)); % Eq. (2.1)
                    Positions(i,j)=Leader_pos(j)-A*D_Leader;      % Eq. (2.2)
                end
                
            elseif p>=0.5
              
                distance2Leader=abs(Leader_pos(j)-Positions(i,j));
                % Eq. (2.5)
                Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);         
            end     
        end
    end
    
    t=t+1;
    Convergence_curve(t)=Leader_score;
    
    if t>2
        line([t-1 t], [Convergence_curve(t-1) Convergence_curve(t)], 'Color','Blue','LineStyle','--','LineWidth', 1.5 )
%         title('Whale Optimization Algorithm');
        xlabel('Number of Iteration');
        ylabel('Best Score Obtained');        
        drawnow
    end
 
    results{1,1} = t;
    results{1,2} = Leader_score;
    
end



