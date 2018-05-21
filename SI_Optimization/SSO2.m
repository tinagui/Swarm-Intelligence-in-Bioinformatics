%==========================================================================
% The Whale Optimization Algorithmm
% Author: Tina Gui (tgui@go.olemiss.edu)
% Modified on Erik Cuevas' original code
%
% Developed in MATLAB R2016a
%
% Reference:  Cuevas, E., Cienfuegos, M., Zaldívar, D., Pérez-Cisneros, M. A 
%   swarm optimization algorithm inspired in the behavior of the social-spider, 
%   Expert Systems with Applications, 40 (16), (2013), pp. 6374-6384 
%   http://arxiv.org/abs/1406.3282
%
%==========================================================================

function SSO2(swarm_no,max_iter,lb_no,ub_no,dim,CostFunction)
%SCO Summary of this function goes here
%   Detailed explanation goes here
    %% Preliminares
	for i=1:dim
        lb(i,:)=lb_no;
        ub(i,:)=ub_no;
    end
    %% Initial Parameters
    rand('state',0');  % Reset the random generator
    % Define the poblation of females and males
    fpl = 0.65;     % Lower Female Percent
    fpu = 0.9;      % Upper Female Percent
    fp = fpl+(fpu-fpl)*rand;	% Aleatory Percent
    fn = round(swarm_no*fp);   % Number of females
    mn = swarm_no-fn;          % Number of males
	%Probabilities of attraction or repulsion
	% Proper tunning for better results
    pm = exp(-(0.1:(3-0.1)/(max_iter-1):3));
    % Initialization of vectors
    fsp = zeros(fn,dim);   % Initlize females
    msp = zeros(mn,dim);   % Initlize males
    fefit = zeros(fn,1);    % Initlize fitness females
    mafit = zeros(mn,1);    % Initlize fitness males
    spwei = zeros(swarm_no,1); % Initlize weigth spiders
    fewei = zeros(fn,1); % Initlize weigth spiders
    mawei = zeros(mn,1); % Initlize weigth spiders

    %% Population Initialization
    % Generate Females
    for i=1:fn
        fsp(i,1:dim)=lb(1)+rand(1,dim).*(ub(1)-lb(1));
    end
    % Generate Males
    for i=1:mn
        msp(i,1:dim)=lb(1)+rand(1,dim).*(ub(1)-lb(1));
    end
    %% **** Evaluations *****
	% Evaluation of function for females
    for i=1:fn
        fefit(i)=CostFunction(fsp(i,:),dim);
    end
	% Evaluation of function for males
    for i=1:mn
        mafit(i)=CostFunction(msp(i,:),dim);
    end
    %% ***** Assign weigth or sort ***********
	% Obtain weight for every spider
    spfit = [fefit' mafit']';   % Mix Females and Males
    bfitw = min(spfit);          % best fitness
    wfit = max(spfit);          % worst fitness
    for i=1:swarm_no
        spwei(i) = 0.001+((spfit(i)-wfit)/(bfitw-wfit));
    end
    fewei = spwei(1:fn);      % Separate the female mass
    mawei = spwei(fn+1:swarm_no);% Separate the male mass
    %% Memory of the best
    % Check the best position
    [~,Ibe] = max(spwei);
    % Check if female or male
    if Ibe > fn
        % Is Male
        spbest=msp(Ibe-fn,:);   % Asign best position to spbest
        bfit = mafit(Ibe-fn);      % Get best fitness for memory
    else
        % Is Female
        spbest=fsp(Ibe,:);      % Asign best position to spbest
        bfit = fefit(Ibe);      % Get best fitness for memory
    end
    %% Start the iterations
    for i=1:max_iter
        %% ***** Movement of spiders *****
        % Move Females
        [fsp] = FeMove(swarm_no,fn,fsp,msp,spbest,Ibe,spwei,dim,lb,ub,pm(i));
        % Move Males
        [msp] = MaMove(fn,mn,fsp,msp,fewei,mawei,dim,lb,ub,pm(i));
        %% **** Evaluations *****
        % Evaluation of function for females
        for j=1:fn
            fefit(j)=CostFunction(fsp(j,:),dim);
        end
        % Evaluation of function for males
        for j=1:mn
            mafit(j)=CostFunction(msp(j,:),dim);
        end
        %% ***** Assign weigth or sort ***********
        spfit = [fefit' mafit']';   % Mix Females and Males
        bfitw = min(spfit);          % best fitness
        wfit = max(spfit);          % worst fitness
        % Obtain weight for every spider
        for j=1:swarm_no
            spwei(j) = 0.001+((spfit(j)-wfit)/(bfitw-wfit));
        end
        fewei = spwei(1:fn);      % Separate the female mass
        mawei = spwei(fn+1:swarm_no);% Separate the male mass
        %% Mating Operator
        [ofspr] = Mating(fewei,mawei,fsp,msp,dim);
        %% Selection of the Mating
        if isempty(ofspr)
%             % Do nothing
        else
            [fsp,msp,fefit,mafit] = Survive(fsp,msp,ofspr,fefit,mafit,spfit,CostFunction,fn,dim);
            % ***** Recalculate the weigth or sort ***********
            spfit = [fefit' mafit']';   % Mix Females and Males
            bfitw = min(spfit);          % best fitness
            wfit = max(spfit);          % worst fitness
            % Obtain weight for every spider
            for j=1:swarm_no
                spwei(j) = 0.001+((spfit(j)-wfit)/(bfitw-wfit));
            end
            fewei = spwei(1:fn);      % Separate the female mass
            mawei = spwei(fn+1:swarm_no);% Separate the male mass
        end
        %% Memory of the best
        % Check if best position belongs to male or female
        [~,Ibe2] = max(spwei);
        if Ibe2 > fn
            % Is Male
            spbest2=msp(Ibe2-fn,:);      % Asign best position to spbest
            bfit2 = mafit(Ibe2-fn);      % Get best fitness for memory
        else
            % Is Female
            spbest2 = fsp(Ibe2,:);  % Asign best position to spbest
            bfit2 = fefit(Ibe2);    % Get best fitness for memory
        end
        %% Global Memory
        if bfit<=bfit2
            bfit = bfit;
            spbest = spbest;      % Asign best position to spbest
            befit(i) = bfit;
        else
            bfit = bfit2;
            spbest = spbest2;      % Asign best position to spbest
            befit(i) = bfit;
        end
        spbesth(i,:)=spbest;
%         %% Plot Results
%         plot(fsp(:,1),fsp(:,2),'r.',msp(:,1),msp(:,2),'bx',spbest(:,1),spbest(:,2),'go');
%         hold on
%         axis([lb(1) ub(1) lb(2) ub(2)])
%         drawnow
%         hold off

        Convergence_curve(i) = bfit;
        if i > 2
            line([i-1 i], [Convergence_curve(i-1) Convergence_curve(i)], 'Color','Red','LineStyle','-.', 'LineWidth', 1.5);
%           title('Social Spider Optimization');
            xlabel('Number of Iteration');
            ylabel('Best Score Obtained');        
        drawnow
        end
    end
    %% Display of results
%     figure
%     plot(befit)
end