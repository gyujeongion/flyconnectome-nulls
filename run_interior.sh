#!/bin/bash
# is the recurrent interior functionally necessary? generation-500 samples of the corrected grid
cd "$(dirname "$0")"
for s in 0 1 2 3 4 5 6 7 8 9; do
  for cfg in "0.0 0.0" "0.2 0.0" "0.0 1.0" "0.2 1.0"; do
    pv=${cfg%% *}; pt=${cfg##* }; run=runs/ecofix_P${pv}_T${pt}_s${s}
    [ -f $run/interior.flag ] && continue
    ./venv/bin/python evo_int.py --homeo global --pop 2560 --conds A,N1,N2,N4,N5 --seed $s --assay $run --assay_gen 500 --assay_set interior \
      --gain_bad 0.25 --pred_fixed $pv --ptox0 $pt >> ${run}_interior.log 2>&1 && touch $run/interior.flag && echo "OK $run" || echo "FAIL $run"
  done
done
echo INTERIOR_DONE
