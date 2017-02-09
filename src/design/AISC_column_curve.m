function PnOverPo = AISC_column_curve(PoOverPe)

PnOverPo = nan(size(PoOverPe));

ind = PoOverPe >= 0.0 & PoOverPe < 2.25;
PnOverPo(ind) = 0.658.^PoOverPe(ind);

ind = PoOverPe >= 2.25;
PnOverPo(ind) = 0.877./PoOverPe(ind);

end