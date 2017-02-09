function K = nomographK_sidesway_inhibited(Ga,Gb)

% Check for known solutions
if Ga == Inf && Gb == Inf
    K = 1.0;
    return
elseif Ga == 0 && Gb == 0
    K = 0.5;
    return
end

% Change Inf to a large number that wont return an error
if Ga == Inf; Ga = 1e10; end
if Gb == Inf; Gb = 1e10; end

% Solve for K
options.Display = 'off';
options.TolX = 1e-6;

if Ga == 0
    [K,~,exitflag] = fsolve(@(K)fcn_with_one_zero(K,Gb),0.7,options);
elseif Gb == 0
    [K,~,exitflag] = fsolve(@(K)fcn_with_one_zero(K,Ga),0.7,options);
else
    [K,~,exitflag] = fsolve(@(K)fcn(K,Ga,Gb),0.7,options);
end
    
if (exitflag < 1)
    error('Could not deterine K'); 
end

end


function num = fcn(K,Ga,Gb)
pK = pi/K;
num = (Ga*Gb/4)*pK^2 + ...
    ((Ga+Gb)/2)*(1-pK/tan(pK)) + ...
    2*tan(pK/2)/(pK)-1;
end

function num = fcn_with_one_zero(K,G)
pK = pi/K;
num = (G/2)*(1-pK/tan(pK)) + ...
    2*tan(pK/2)/(pK)-1;
end
