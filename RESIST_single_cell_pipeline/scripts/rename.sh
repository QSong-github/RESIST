
cd /blue/qsong1/sen.guo/resist_02/resist_test_share/results


for f in GSE104987_seurat_afterAnno.RDS_* GSE104987_deg.sig_* GSE104987_*; do
  [ -e "$f" ] || continue
  new="${f#GSE104987_seurat_afterAnno.RDS_}"
  new="${new#GSE104987_deg.sig_}"
  new="${new#GSE104987_}"
  [ "$new" != "$f" ] || continue
  if [ -e "$new" ]; then
    echo "SKIP (target exists): $f -> $new"
    continue
  fi
  mv -- "$f" "$new"
done


