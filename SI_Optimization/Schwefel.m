
% Schwefel
function Cost = Schwefel( x, d )
    for i = 1:d;
        Cost=sum(x(i)*sin(sqrt(abs(x(i))))); 
    end
end
