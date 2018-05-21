% Step
function Cost = Step( x, d )
    for i = 1:d;
        Cost=sum((x(i)+0.5).^2); 
    end
end
