conda activate pyscenic
pyscenic grn --num_workers 16 --output N1_adj.tsv  --method grnboost2 loom/NC_matrix.loom LYR/hs_hgnc_tfs.txt

pyscenic ctx \
N1_adj.tsv \
LYR/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname LYR/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname loom/NC_matrix.loom \
--mode "dask_multiprocessing" \
--output N1_scenic.csv \
--num_workers 16 \
--mask_dropouts

pyscenic aucell \
loom/NC_matrix.loom \
N1_scenic.csv \
--output N1_SCENIC.loom \
--num_workers 16


pyscenic grn --num_workers 16 --output IR2h_adj.tsv  --method grnboost2 loom/IR_2h_matrix.loom LYR/hs_hgnc_tfs.txt
pyscenic ctx IR2h_adj.tsv LYR/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather --annotations_fname LYR/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname loom/IR_2h_matrix.loom --mode "dask_multiprocessing" --output IR2h.csv --num_workers 16 --mask_dropouts
pyscenic aucell loom/IR_2h_matrix.loom IR2h.csv --output IR2h_SCENIC.loom --num_workers 16 

pyscenic grn --num_workers 16 --output int_adj.tsv  --method grnboost2 loom/int_matrix.loom LYR/hs_hgnc_tfs.txt
pyscenic ctx int_adj.tsv LYR/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather --annotations_fname LYR/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname loom/int_matrix.loom --mode "dask_multiprocessing" --output int.csv --num_workers 16 --mask_dropouts
pyscenic aucell loom/int_matrix.loom int.csv --output int_SCENIC.loom --num_workers 16
