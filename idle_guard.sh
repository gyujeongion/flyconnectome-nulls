#!/bin/bash
# Keep the evolution runs out of the way of music work on this machine.
#
# A run costs little CPU and RAM (measured: 2.5 GB of 64, two thirds of one core of twenty)
# because the work is on the GPU. That measurement is misleading on its own: the window
# server composites the screen on the same GPU, so a run that looks free by CPU and RAM
# still makes the pointer stutter. Observed on 2026-09-22.
#
# So the default is cautious. The run is suspended whenever a DAW is open OR somebody has
# touched the keyboard or mouse in the last 30 seconds, and only resumes after 5 minutes of
# quiet. Set ALSO_ON_HID=0 to keep running through keyboard activity.
#
# Suspending is SIGSTOP on the run: it keeps its memory, uses no CPU, and resumes exactly
# where it stopped. Nothing is recomputed.
#
#   ./idle_guard.sh &                  # pause while a DAW is running
#   ALSO_ON_HID=1 ./idle_guard.sh &    # pause on any keyboard/mouse activity too
#
PATTERN="${PATTERN:-[e]vo_[a-z]*\.py}"          # the runs to supervise
DAW_PATTERN="${DAW_PATTERN:-a digital audio workstation|Logic Pro|Pro Tools|REAPER|Cubase|FL Studio|Studio One}"
ALSO_ON_HID="${ALSO_ON_HID:-1}"                 # 1 = also pause while the machine is in use (default)
PAUSE_BELOW="${PAUSE_BELOW:-30}"                # (HID mode) idle below this -> pause
RESUME_ABOVE="${RESUME_ABOVE:-300}"             # (HID mode) idle above this -> resume
POLL="${POLL:-5}"

DIR="$(cd "$(dirname "$0")" && pwd)"
LOG="$DIR/idle_guard.log"
STATE="$DIR/idle_guard.state"

idle_seconds() { ioreg -c IOHIDSystem 2>/dev/null | awk '/HIDIdleTime/ {printf "%d", $NF/1000000000; exit}'; }
run_pids()     { pgrep -f "$PATTERN" 2>/dev/null; }
daw_running()  { pgrep -f "$DAW_PATTERN" >/dev/null 2>&1; }
log()          { printf '%s %s\n' "$(date '+%F %T')" "$*" >> "$LOG"; }

state=running
paused_since=0
total_paused=0
started=$(date +%s)

log "GUARD START daw=\"$DAW_PATTERN\" also_on_hid=$ALSO_ON_HID"
trap 'p=$(run_pids); [ -n "$p" ] && kill -CONT $p 2>/dev/null; log "GUARD STOP (resumed any paused run)"; exit 0' INT TERM

while true; do
  i=$(idle_seconds); [ -z "$i" ] && i=999999
  p=$(run_pids)
  reason=""

  # A DAW being open is reason enough on its own; keyboard activity only counts if asked for.
  if daw_running; then
    reason="daw"
  elif [ "$ALSO_ON_HID" = 1 ] && [ "$i" -lt "$PAUSE_BELOW" ]; then
    reason="hid"
  fi

  # Resuming needs the DAW gone, and (in HID mode) the machine idle for long enough that we
  # are not going to flap the moment somebody nudges the mouse.
  can_resume=1
  daw_running && can_resume=0
  [ "$ALSO_ON_HID" = 1 ] && [ "$i" -lt "$RESUME_ABOVE" ] && can_resume=0

  if [ -n "$p" ]; then
    if [ "$state" = running ] && [ -n "$reason" ]; then
      kill -STOP $p 2>/dev/null
      state=paused; paused_since=$(date +%s)
      log "PAUSE  reason=$reason idle=${i}s pids=$(echo $p | tr '\n' ' ')"
    elif [ "$state" = paused ] && [ "$can_resume" = 1 ]; then
      kill -CONT $p 2>/dev/null
      now=$(date +%s); total_paused=$((total_paused + now - paused_since))
      state=running
      log "RESUME idle=${i}s paused_for=$((now - paused_since))s total_paused=${total_paused}s"
    fi
  fi

  now=$(date +%s); cur=$total_paused
  [ "$state" = paused ] && cur=$((cur + now - paused_since))
  printf 'state=%s idle=%ss daw=%s elapsed=%ss paused=%ss compute=%ss\n' \
    "$state" "$i" "$(daw_running && echo yes || echo no)" \
    "$((now - started))" "$cur" "$((now - started - cur))" > "$STATE"

  sleep "$POLL"
done
