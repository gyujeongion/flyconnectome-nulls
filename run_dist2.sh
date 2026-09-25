#!/bin/bash
# PROTOCOL_DIST2.md: the corrected grid's four ecologies, seeds 0-9, re-run with distribution-matched
# calibration. Same frozen engine and flags as run_fix.sh; only --homeo changes (global -> distmatch).
# Starts after run_followup.sh finishes. Never starts a run after the holiday window closes, and waits
# while a digital audio workstation is open so the machine is not loaded while someone is working on it.
cd "$(dirname "$0")"
log() { printf "%s %s\n" "$(date "+%F %T")" "$*" >> dist2.log; }
# To resume in a later idle window: DIST2_DEADLINE="YYYY-MM-DD HH:MM" ./run_dist2.sh
DEADLINE=$(date -j -f "%Y-%m-%d %H:%M" "${DIST2_DEADLINE:-2026-09-28 05:15}" +%s)
ENGINE_HASH=c0adad2f709e93b0

[ "$(shasum -a 256 evo_fix.py | cut -c1-16)" = "$ENGINE_HASH" ] || { log "ABORT engine hash differs from PROTOCOL_DIST2"; exit 1; }

log "waiting for run_followup.sh to finish"
until grep -q "FOLLOWUP_DONE" followup.log 2>/dev/null; do sleep 120; done
log "follow-up finished; starting distmatch grid"

for s in 0 1 2 3 4 5 6 7 8 9; do
  for cfg in "0.0 0.0" "0.2 0.0" "0.0 1.0" "0.2 1.0"; do
    pv=${cfg%% *}; pt=${cfg##* }; run=ecodist_P${pv}_T${pt}_s${s}
    [ -f "runs/$run/done.flag" ] && { log "SKIP $run"; continue; }
    while pgrep -xq "$DAW_PROCESS"; do
      [ "$(date +%s)" -ge "$DEADLINE" ] && break
      log "a digital audio workstation is open; waiting"; sleep 300
    done
    if [ "$(date +%s)" -ge "$DEADLINE" ]; then log "DEADLINE reached before $run; stopping (resume later with the same script)"; exit 0; fi
    rm -rf "runs/$run"
    log "RUN $run"
    if ./venv/bin/python evo_fix.py --homeo distmatch --pop 2560 --conds A,N1,N2,N4,N5 --seed "$s" --gens 600 \
         --probe_every 50 --probe_reps 8 --ptox0 "$pt" --probe_ptox "$pt" --pred_fixed "$pv" \
         --gain_bad 0.25 --tox_curriculum 0 --migrate_every 0 \
         --run "$run" > "runs_$run.log" 2>&1; then
      touch "runs/$run/done.flag"; log "OK $run"
    else
      log "FAIL $run"
    fi
  done
  log "SEED_DONE $s"
done
log "DIST2_DONE"
