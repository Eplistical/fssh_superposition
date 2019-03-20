#!/bin/bash

Alist="
0.02
0.05
0.10
"

plist=`python -c "for k in range(10,36): print('%d' % k)"`
slist=`python -c "for k in range(21): print('%.2f' % (k * 0.05))"`


rm *.out

for A in $Alist; do
for p in $plist; do
    cat scan_p_ex_raw/e2.A$A.s1.00.px$p.out >> e.A$A.scan_p.out
    cat scan_p_fssh_raw/f2.A$A.s1.00.px$p.out >> f.A$A.scan_p.out
    cat scan_p_fssh_phase_corr_raw/fc2.A$A.s1.00.px$p.out >> fc.A$A.scan_p.out

    cat scan_p_ex_model2_raw/e2.A$A.s1.00.px$p.out >> e.A$A.scan_p.model2.out
    cat scan_p_fssh_model2_raw/f2.A$A.s1.00.px$p.out >> f.A$A.scan_p.model2.out
    cat scan_p_fssh_model2_phase_corr_raw/f2.A$A.s1.00.px$p.out >> fc.A$A.scan_p.model2.out
done
done

for A in $Alist; do
for s in $slist; do
    cat scan_s_ex_raw/e2.A$A.s$s.px30.out >> e.A$A.scan_s.out
    cat scan_s_fssh_raw/f2.A$A.s$s.px30.out >> f.A$A.scan_s.out
    cat scan_s_fssh_phase_corr_raw/fc2.A$A.s$s.px30.out >> fc.A$A.scan_s.out

    cat scan_s_ex_model2_raw/e2.A$A.s$s.px30.out >> e.A$A.scan_s.model2.out
    cat scan_s_fssh_model2_raw/f2.A$A.s$s.px30.out >> f.A$A.scan_s.model2.out
    cat scan_s_fssh_model2_phase_corr_raw/f2.A$A.s$s.px30.out >> fc.A$A.scan_s.model2.out

    cat scan_s_ehrenfest_raw/eh2.A$A.s$s.px30.out >> eh.A$A.scan_s.out
    cat scan_s_ehrenfest_model2_raw/eh2.A$A.s$s.px30.out >> eh.A$A.scan_s.model2.out
done
done

mv *.out adhoc/
