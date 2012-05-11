# Author : Sukhdeep Singh
# Usage  : ~/src/shell/foldEnrichment.sh fdr_peaks.bed
# Usage2 : ~/src/shell/foldEnrichment.sh macs_peaks.xls -p

# this script gives the total fold enrichment of whole sample profile wrt to the used control, useful for comparing different controls
if [ $# -ge 2 ]; then
	sed -e '1,24d' $1| sort -k8n $1 | awk '{total=total+$8} END {print total}'
else
	sort -k8n $1 | awk '{total=total+$8} END {print total}'
fi
