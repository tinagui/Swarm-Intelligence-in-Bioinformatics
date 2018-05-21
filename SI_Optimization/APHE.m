
% Axis Parallel Hyper-Ellipsoid
function Cost = APHE( x, d )
    for i = 1:d;
        Cost= sum(d*x(i).^2); 
    end
end

