#!/bin/bash
# Stops a Gamma-path SCF that is clearly diverging (>= 20 iterations, last estimated accuracy > 1000 Ry)
# so that pipeline.sh falls back to the complex k path at once. Only touches runs of disp-NNNNN.in (not .cplx).
cd /home/koussay/blockc/pbesol
while true; do
  for f in disp-?????.out.part; do
    [ -f "$f" ] || continue
    inp=${f%.out.part}.in
    pgrep -f "pw.x -nk 1 -in $inp\$" >/dev/null || continue
    n=$(grep -c "iteration #" "$f")
    acc=$(grep "estimated scf accuracy" "$f" | tail -n 1 | awk '{print $5}')
    if [ "$n" -ge 20 ] && awk -v a="$acc" 'BEGIN{exit !(a>1000)}'; then
      echo "$(date -u +%FT%TZ) gwatch stop $inp it=$n acc=$acc" >> /home/koussay/blockc/progress.log
      kill $(pgrep -x pw.x)
    fi
  done
  sleep 120
done
