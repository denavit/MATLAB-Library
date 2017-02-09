function PoOverPe = AISC_inverse_column_curve(PnOverPo)

PoOverPe = nan(size(PnOverPo));

ind = PnOverPo <= 1.0 & PnOverPo > 0.658^2.25;
PoOverPe(ind) = log(PnOverPo(ind))/log(0.658);

ind = PnOverPo <= 0.658^2.25 & PnOverPo > 0.0;
PoOverPe(ind) = 0.877./PnOverPo(ind);

ind = PnOverPo == 0.0;
PoOverPe(ind) = Inf;

end