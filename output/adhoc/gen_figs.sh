#!/bin/bash

pylist="
fig_surf.py
fig_scan_p.py
fig_scan_s.py
fig_avgp.py
fig_surf_model2.py
fig_avgp_model2.py
fig_scan_p_model2.py
fig_scan_s_model2.py
fig_minhop.py
"

for s in $pylist; do
    echo $s
    python3 $s
done
