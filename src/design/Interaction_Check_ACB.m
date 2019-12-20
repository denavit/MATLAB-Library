function ratio = Interaction_Check_ACB(Pr,Mrz,Mry,Pa,Pc,Mcz,Mcy,alpha)
    assert(all(Pr >= 0.0,'all'),'Pr should be positive, indicating compression')
    assert(and(Pa > Pc,Pc > 0.0),'Need Pa > Pc > 0 (Pa = %g, Pc = %g)',Pa,Pc)
    assert(Mcz > 0.0,'Mcz should be positive (Mcz = %g)',Mcz)
    assert(Mcy > 0.0,'Mcy should be positive (Mcy = %g)',Mcy)
    if nargin < 8
        alpha = 1;
    end
    ind = Pr > Pc;
    ratio(ind) = (Pr(ind)-Pc)/(Pa-Pc)+((abs(Mrz(ind))/Mcz).^alpha + (abs(Mry(ind))/Mcy).^alpha).^(1/alpha);
    ind = Pr < Pc;
    ratio(ind) = ((abs(Mrz(ind))/Mcz).^alpha + (abs(Mry(ind))/Mcy).^alpha).^(1/alpha); 
end