function phi = ACI_phi(type,et,ety)

switch lower(type)
    case 'spiral'
        phicc = 0.75;
    case 'other'
        phicc = 0.65;
    otherwise
        error('Unknown type %s',type)
end

if strcmpi(ety,'Grade60')
    ety = 0.002;
end

phitc = 0.9;

phi = zeros(size(et));
% Tension Controlled
ind = 0.005 <= et;
phi(ind) = phitc;
% Transition
ind = ety < et & et < 0.005;
phi(ind) = phicc + (phitc-phicc)*(et(ind)-ety)/(0.005-ety);
% Compression Controlled
ind = et <= ety;
phi(ind) = phicc;

end