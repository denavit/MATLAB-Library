function K = nomographK_sidesway_uninhibited(Ga,Gb)

% Check for known solutions
if Ga == Inf && Gb == Inf
    K = Inf;
    return
elseif (Ga == Inf && Gb == 0) || (Ga == 0 && Gb == Inf)
    K = 2.0;
    return
elseif Ga == 0 && Gb == 0
    K = 1.0;
    return
end

% Change Inf to a large number that wont return an error
if Ga == Inf; Ga = 1e10; end
if Gb == Inf; Gb = 1e10; end

% Solve for K
options.Display = 'off';
options.TolX = 1e-6;
   
if Ga == 0
    [K,~,exitflag] = fsolve(@(K)fcn_with_one_zero(K,Gb),2.0,options);
elseif Gb == 0
    [K,~,exitflag] = fsolve(@(K)fcn_with_one_zero(K,Ga),2.0,options);
else
    [K,~,exitflag] = fsolve(@(K)fcn(K,Ga,Gb),2.0,options);
end

if (exitflag < 1)
    error('Could not determine K'); 
end

if K < 1
    error('K found by fsolve is not correct');
    % If this becomes a problem then try using fmincon
end

end


function num = fcn(K,Ga,Gb)
pK = pi/K;
num = (Ga*Gb*pK^2 - 36)/(6*(Ga+Gb)) - pK/tan(pK);
end

function num = fcn_with_one_zero(K,G)
pK = pi/K;
num = -6/G - pK/tan(pK);
end
