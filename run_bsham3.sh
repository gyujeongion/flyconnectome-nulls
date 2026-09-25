#!/bin/bash
# Replaces run_bsham2.sh (research log 56). Driver 2 decided "finished" from the last log line being "done", but
# mlx_capped.py prints its memory line after the engine's "done", so every completed capped run was marked FAIL
# and deleted by the retry. Same runs, engine, flags and launcher; differences: "done" is looked for anywhere in
# the log, a completed run is never deleted, and a failed attempt's log is kept as runs_<run>.log.fail<attempt>
# together with its exit code.
cd "$(dirname "$0")"
export TMPDIR=/path/to/external-ssd/.tmp MLX_CACHE_GB=2
MIRROR="$HOME/.local/share/fly_notify/bsham.log"
log() { printf "%s %s\n" "$(date "+%F %T")" "$*" | tee -a bsham.log >> "$MIRROR"; }
ENGINE_HASH=$(cat .bsham_engine_hash 2>/dev/null)
[ -n "$ENGINE_HASH" ] && [ "$(shasum -a 256 evo_bsham.py | cut -c1-16)" = "$ENGINE_HASH" ] || { log "ABORT engine hash differs from PROTOCOL_BSHAM"; exit 1; }
finished() { [ -f "runs/$1/fresh.jsonl" ] && [ "$(wc -l < "runs/$1/fresh.jsonl")" -eq 4 ] && grep -qx "done" "runs_$1.log"; }
peak() { grep -o 'peak [0-9.]* GB' "runs_$1.log" | tail -1; }

if [ -n "${INFLIGHT_PID:-}" ] && [ -n "${INFLIGHT_RUN:-}" ]; then
  log "driver 3: waiting for in-flight $INFLIGHT_RUN (pid $INFLIGHT_PID)"
  while kill -0 "$INFLIGHT_PID" 2>/dev/null; do sleep 60; done
  if finished "$INFLIGHT_RUN"; then touch "runs/$INFLIGHT_RUN/done.flag"; log "OK $INFLIGHT_RUN $(peak "$INFLIGHT_RUN")"
  else mv "runs_$INFLIGHT_RUN.log" "runs_$INFLIGHT_RUN.log.fail_inflight"; log "FAIL $INFLIGHT_RUN (in-flight run did not finish; log kept)"; fi
fi
log "driver 3: memory-capped launcher, retries up to 2, failed logs kept"
for s in 0 1 2 3 4 5 6 7 8 9; do
  for pv in 0.0 0.2; do
    run=bsham_P${pv}_T0.0_s${s}
    [ -f "runs/$run/done.flag" ] && continue
    for attempt in 1 2 3; do
      rm -rf "runs/$run"
      log "RUN $run (attempt $attempt)"
      ./venv/bin/python mlx_capped.py evo_bsham.py --homeo global --pop 4608 --conds A,AS1,AS3,AS5,AS10,BS1,BS3,BS5,BS10 --seed "$s" --gens 600 \
           --probe_every 50 --probe_reps 8 --ptox0 0.0 --probe_ptox 0.0 --pred_fixed "$pv" \
           --gain_bad 0.25 --tox_curriculum 0 --migrate_every 0 \
           --run "$run" > "runs_$run.log" 2>&1
      rc=$?
      if finished "$run"; then
        touch "runs/$run/done.flag"; log "OK $run $(peak "$run")"; break
      fi
      mv "runs_$run.log" "runs_$run.log.fail$attempt"
      log "FAIL $run attempt $attempt: exit $rc $(grep -m1 -oE 'RuntimeError: .{0,90}|Error: .{0,90}' "runs_$run.log.fail$attempt")"
      sleep 120
    done
  done
  log "SEED_DONE $s"
done
log "BSHAM_DONE"
