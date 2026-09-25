#!/bin/bash
# PROTOCOL_DOSE.md: graded olfactory shortcut share + rewiring-matched sham, two ecologies, 10 seeds, 600 generations
cd ~/Development/260917_fly_evolution
for s in 0 1 2 3 4 5 6 7 8 9; do
  for cfg in "0.0 0.0" "0.2 0.0"; do
    pv=${cfg%% *}; pt=${cfg##* }; run=dose_P${pv}_T${pt}_s${s}
    [ -f runs/$run/done.flag ] && continue
    rm -rf runs/$run
    ./venv/bin/python evo_dose.py --homeo global --pop 3072 --conds A,AS1,AS3,AS5,AS10,SHAM --seed $s --gens 600 --probe_every 50 --probe_reps 8 \
      --ptox0 $pt --probe_ptox $pt --pred_fixed $pv --gain_bad 0.25 --tox_curriculum 0 --migrate_every 0 \
      --run $run > runs_$run.log 2>&1 && touch runs/$run/done.flag && echo "OK $run" || echo "FAIL $run"
  done
  echo "SEED_DONE $s"
done
echo DOSE_DONE
