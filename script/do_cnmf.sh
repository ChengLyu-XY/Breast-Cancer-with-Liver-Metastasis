#! /bin/bash/
root_dir=/alldata/lc/breast/MP/all_counts #The folder contains count matirix formed as genes × cells of all samples
outdir=/alldata/lc/breast/MP/cnmf_res
ls ${root_dir}|while read id
do
	id="${id::-11}"
	cnmf prepare --output-dir ${outdir} --name ${id} -c ${root_dir}/${id}_counts.txt -k 4 5 6 7 8 9 10 --n-iter 100 --seed 14
	cnmf factorize --output-dir ${outdir} --name ${id} --worker-index 0 --total-workers 1
	cnmf combine --output-dir ${outdir} --name ${id}
	rm -f ${outdir}/${id}/cnmf_tmp/${id}.spectra.k_*.iter_*.df.npz
	cnmf k_selection_plot --output-dir ${outdir} --name ${id}
done