finalfolder="~/Dropbox/projects/MuRF1/October"
figfolder="~/projects/MuRF1/plots"
mkdir -p $finalfolder


rsync -vr  ~/projects/MuRF1/plots/geneexprplots/* $finalfolder"/exprplots"

rsync -vr  ~/projects/MuRF1/tables/rawcountdata.tsv $finalfolder"/exprplots"
rsync -vr  ~/projects/MuRF1/tables/tpms.tsv $finalfolder"/exprplots"
rsync -vr  ~/projects/MuRF1/tables/allcontrasts.tsv $finalfolder"/exprplots"

find  $finalfolder -type f ! -iname '*DS_Store*' -exec realpath {} \;


finalfolder="/Users/dharnet/Dropbox/projects/MuRF1/December"
mkdir -p $finalfolder
figfolder="/fast/work/groups/ag_ohler/dharnet_m//MuRF1/plots"

mkdir -p $finalfolder"/locusplots"
mkdir -p $finalfolder"/dtuplots"

rsync -vr  $figfolder/locusplot*.pdf $finalfolder"/locusplots"
rsync -vr  $figfolder/dtuplots/* $finalfolder"/dtuplots"

rsync -vr '/fast/work/groups/ag_ohler/dharnet_m//MuRF1/src/rnaseq_de_report.nb.html' $finalfolder
rsync -vr '/fast/work/groups/ag_ohler/dharnet_m//MuRF1/tables' $finalfolder




finalfolder="/Users/dharnet/Dropbox/projects/MuRF1/January"
figfolder="/fast/work/groups/ag_ohler/dharnet_m/MuRF1/plots"
mkdir -p $finalfolder
rsync -vr  $figfolder/geneexprplots/* $finalfolder"/exprplots"
rsync -vr  $figfolder/locusplot*.pdf $finalfolder"/locusplots"
rsync -vr  $figfolder/dtuplots/* $finalfolder"/dtuplots"
