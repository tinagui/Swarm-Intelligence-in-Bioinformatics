
% Griewank
function Cost = Griewank(x, d)
    % Griewank Function
    % [-600 600]
    % The global minima: x* =  (0, …, 0), f(x*) = 0.
    n = d;
    fr = 4000;
    s = 0;
    p = 1;
    for i = 1:n;
        s = s+x(i)^2; 
    end
    for j = 1:n;
        p = p*cos(x(j)/sqrt(j)); 
    end
    Cost = s/fr-p+1;
end
