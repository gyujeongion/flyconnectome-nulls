#!/bin/bash
# After the campaign driver finishes, extend two registered grids to seeds 10-19 with each
# grid's own frozen engine and its original launch flags, copied from run_fix.sh and
# run_dose.sh. The first attempt at the corrected-grid extension went through the campaign
# driver with the newer engine, left --probe_ptox at that engine's default of 0.5 and so could
# not be pooled with seeds 0-9 (research log 45-2); those runs are kept as runs/ecofix_pp05_*.
# Order: corrected grid first (it answers the per-ecology equivalence question), then dose.
cd "$(dirname "$0")"
log() { printf "%s %s\n" "$(date "+%F %T")" "$*" >> followup.log; }

log "waiting for the campaign driver to finish"
until grep -q "CAMPAIGN COMPLETE" chuseok.log 2>/dev/null; do sleep 120; done
log "driver finished; starting extensions"

for grid in fix dose; do
  for s in 10 11 12 13 14 15 16 17 18 19; do
    for cfg in "0.0 0.0" "0.2 0.0"; do
      pv=${cfg%% *}; pt=${cfg##* }
      if [ "$grid" = fix ]; then
        run=ecofix_P${pv}_T${pt}_s${s}; eng=evo_fix.py; conds=A,N1,N2,N4,N5; pop=2560
      else
        run=dose_P${pv}_T${pt}_s${s}; eng=evo_dose.py; conds=A,AS1,AS3,AS5,AS10,SHAM; pop=3072
      fi
      [ -f "runs/$run/done.flag" ] && { log "SKIP $run"; continue; }
      rm -rf "runs/$run"
      log "RUN $run ($eng)"
      if ./venv/bin/python "$eng" --homeo global --pop "$pop" --conds "$conds" --seed "$s" --gens 600 \
           --probe_every 50 --probe_reps 8 --ptox0 "$pt" --probe_ptox "$pt" --pred_fixed "$pv" \
           --gain_bad 0.25 --tox_curriculum 0 --migrate_every 0 \
           --run "$run" > "runs_$run.log" 2>&1; then
        touch "runs/$run/done.flag"; log "OK $run"
      else
        log "FAIL $run"
      fi
    done
  done
  log "GRID_DONE $grid"
done
log "FOLLOWUP_DONE"
