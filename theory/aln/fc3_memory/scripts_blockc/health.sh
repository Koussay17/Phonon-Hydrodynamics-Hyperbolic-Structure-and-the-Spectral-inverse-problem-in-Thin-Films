#!/bin/bash
# One-line health summary of the Block C pipeline (read-only).
cd /home/koussay/blockc
echo "last: $(tail -1 progress.log)"
echo "pw.x=$(pgrep -c pw.x) pbe=$(ls pbe/*.out 2>/dev/null | wc -l)/467 pbesol=$(ls pbesol/*.out 2>/dev/null | wc -l)/467 bad=$(grep -cE 'FAILED|rc=[1-9]' progress.log)"
for f in pbe/*.part pbesol/*.part; do
  [ -f "$f" ] && echo "running $f since $(( ($(date +%s) - $(b=$(stat -c %W "$f"); [ "$b" -gt 0 ] && echo $b || stat -c %Y "$f")) / 60 )) min"
done
