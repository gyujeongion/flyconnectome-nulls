#!/bin/bash
# PROTOCOL_BSHAM.md: dose series + boundary-targeted shortcut-free shams, two ecologies, seeds 0-9, 600 generations.
# Runs on the home Mac Studio, one run at a time, everything on the external SSD. The log is mirrored to a small
# file on the internal disk so the launchd notifier (which cannot read external volumes) can see progress.
cd "$(dirname "$0")"
export TMPDIR=/path/to/external-ssd/.tmp
MIRROR="$HOME/.local/share/fly_notify/bsham.log"; mkdir -p "$(dirname "$MIRROR")"
log() { printf "%s %s\n" "$(date "+%F %T")" "$*" | tee -a bsham.log >> "$MIRROR"; }
ENGINE_HASH=$(cat .bsham_engine_hash 2>/dev/null)
[ -n "$ENGINE_HASH" ] && [ "$(shasum -a 256 evo_bsham.py | cut -c1-16)" = "$ENGINE_HASH" ] || { log "ABORT engine hash differs from PROTOCOL_BSHAM"; exit 1; }

log "starting boundary-sham grid"
for s in 0 1 2 3 4 5 6 7 8 9; do
  for pv in 0.0 0.2; do
    run=bsham_P${pv}_T0.0_s${s}
    [ -f "runs/$run/done.flag" ] && { log "SKIP $run"; continue; }
    rm -rf "runs/$run"
    log "RUN $run"
    if ./venv/bin/python evo_bsham.py --homeo global --pop 4608 --conds A,AS1,AS3,AS5,AS10,BS1,BS3,BS5,BS10 --seed "$s" --gens 600 \
         --probe_every 50 --probe_reps 8 --ptox0 0.0 --probe_ptox 0.0 --pred_fixed "$pv" \
         --gain_bad 0.25 --tox_curriculum 0 --migrate_every 0 \
         --run "$run" > "runs_$run.log" 2>&1; then
      touch "runs/$run/done.flag"; log "OK $run"
    else
      log "FAIL $run"
    fi
  done
  log "SEED_DONE $s"
done
log "BSHAM_DONE"
