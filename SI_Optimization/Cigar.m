
% Cigar
function Cost = Cigar( x, d )
    x_0 = x(1).^2;
    for i = 1:d;
        Cost=x_0 + sum(d*x(i).^2); 
    end
end
