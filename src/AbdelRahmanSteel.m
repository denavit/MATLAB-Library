function stress = AbdelRahmanSteel(strain,Fy,Es,type,Fu,ri_over_t)
if nargin < 4
    type = 'flat';
end

% Check Data
assert(Fy > 0 ,'Fy should be greater than zero');
assert(Es > 0 ,'Es should be greater than zero');

% Initilize Data
stress = zeros(size(strain));

switch type
    case 'flat'
        
    case 'corner'
        assert(Fu > Fy,'Fu should be greater than Fy');
        assert(ri_over_t >= 0,'r_over_t should be positive');
        Bc = 3.69*(Fu/Fy) - 0.819*(Fu/Fy)^2 - 1.79;
        m  = 0.192*(Fu/Fy) - 0.068;
        Fy = Fy + 0.6*(Bc/(ri_over_t)^m - 1.0)*Fy;
    otherwise
        error('Unknown type: %s',type);
end

% Material Properties
Fp = 0.75*Fy;
Fym = 0.875*Fy;

Et1 = Es/2;
Et2 = Es/10;
Et3 = Es/200;

e1 = Fp/Es;
e2 = e1+(Fym-Fp)/Et1;
e3 = e2+(Fy-Fym)/Et2;

% Compute Stress
ind =  e3 <= strain;
stress(ind) =  Fy  + Et3*(strain(ind)-e3);
ind =  e2 <= strain & strain <  e3;
stress(ind) =  Fym + Et2*(strain(ind)-e2);
ind =  e1 <= strain & strain <  e2;
stress(ind) =  Fp  + Et1*(strain(ind)-e1);
ind = -e1 <= strain & strain <  e1;
stress(ind) =  Es*strain(ind);
ind = -e2 <= strain & strain < -e1;
stress(ind) = -Fp  + Et1*(strain(ind)+e1);
ind = -e3 <= strain & strain < -e2;
stress(ind) = -Fym + Et2*(strain(ind)+e2);
ind = strain < -e3;
stress(ind) = -Fy  + Et3*(strain(ind)+e3);

end 
