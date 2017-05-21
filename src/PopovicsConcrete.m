function stress = PopovicsConcrete(strain,fcc,ecc,Ec)

% Check Data
assert(fcc < 0,'fcc should be negative');
assert(ecc < 0,'ecc should be negative');
assert( Ec > 0,'Ec should be positive');

% Initilize Data
stress = zeros(size(strain));

% Compute Stress (Compression)
ind = strain < 0;
x = strain(ind)/ecc;
n = Ec*ecc/fcc;
r = n/(n-1);
y = (r*x)./(r-1+x.^r);
stress(ind) = fcc*y;

end 
