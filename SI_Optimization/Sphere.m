
% Sphere
function Cost = Sphere( x, d )
    for i = 1:d;
        Cost=sum(x(i).^2); 
    end
end

