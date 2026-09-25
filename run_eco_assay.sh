#!/bin/bash
# cross-ecology assay on generation-500 elites; starts only after the evolution grid finishes
cd ~/Development/260917_fly_evolution
until grep -q ECO_DONE eco.log; do sleep 120; done
for s in 0 1 2 3 4 5 6 7 8 9; do
  for cfg in "0.0 0.0" "0.2 0.0" "0.0 1.0" "0.2 1.0"; do
    pv=${cfg%% *}; pt=${cfg##* }; run=runs/eco_P${pv}_T${pt}_s${s}
    [ -f $run/assay_eco.flag ] && continue
    ./venv/bin/python evo_assay.py --homeo global --seed $s --assay $run --assay_gen 500 --assay_set eco \
      --gain_bad 0.25 --pred_fixed $pv --ptox0 $pt > ${run}_assay.log 2>&1 && touch $run/assay_eco.flag && echo "OK $run" || echo "FAIL $run"
  done
done
echo ASSAY_DONE
