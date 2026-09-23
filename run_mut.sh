#!/bin/bash
# PROTOCOL_MUT.md: mutational-neighbourhood assay on saved genotypes, generation 0/250/500, home ecology
cd "$(dirname "$0")"
for s in 0 1 2 3 4 5 6 7 8 9; do
  for cfg in "0.0 0.0" "0.2 0.0" "0.0 1.0" "0.2 1.0"; do
    pv=${cfg%% *}; pt=${cfg##* }; run=runs/eco_P${pv}_T${pt}_s${s}
    for g in 0 250 500; do
      f=$(printf "%s/mut_g%06d.npz" $run $g); [ -f $f ] && continue
      ./venv/bin/python evo_mut.py --homeo global --seed $s --assay $run --assay_gen $g --assay_set mut \
        --gain_bad 0.25 --pred_fixed $pv --ptox0 $pt >> ${run}_mut.log 2>&1 && echo "OK $run g$g" || echo "FAIL $run g$g"
    done
  done
done
echo MUT_DONE
