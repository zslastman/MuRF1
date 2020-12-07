finalfolder="~/Dropbox/projects/MuRF1/October"
figfolder="~/projects/MuRF1/plots"
mkdir -p $finalfolder


rsync -vr  ~/projects/MuRF1/plots/geneexprplots/* $finalfolder"/exprplots"

rsync -vr  ~/projects/MuRF1/tables/rawcountdata.tsv $finalfolder"/exprplots"
rsync -vr  ~/projects/MuRF1/tables/tpms.tsv $finalfolder"/exprplots"
rsync -vr  ~/projects/MuRF1/tables/allcontrasts.tsv $finalfolder"/exprplots"

find  $finalfolder -type f ! -iname '*DS_Store*' -exec realpath {} \;
