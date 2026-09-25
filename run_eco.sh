#!/bin/bash
# paper-v2 ecology grid: seed-major order so an early stop still leaves a balanced 2x2
cd ~/Development/260917_fly_evolution
for s in 0 1 2 3 4 5 6 7 8 9; do
  for cfg in "0.0 0.0" "0.2 0.0" "0.0 1.0" "0.2 1.0"; do
    pv=${cfg%% *}; pt=${cfg##* }
    run=eco_P${pv}_T${pt}_s${s}
    [ -f runs/$run/done.flag ] && continue
    rm -rf runs/$run
    ./venv/bin/python evo_paper.py --homeo global --seed $s --gens 600 --probe_every 50 --probe_reps 8 \
      --ptox0 $pt --probe_ptox $pt --pred_fixed $pv --gain_bad 0.25 --tox_curriculum 0 --migrate_every 0 \
      --run $run > runs_$run.log 2>&1 && touch runs/$run/done.flag && echo "OK $run" || echo "FAIL $run"
  done
  echo "SEED_DONE $s"
done
echo ECO_DONE
