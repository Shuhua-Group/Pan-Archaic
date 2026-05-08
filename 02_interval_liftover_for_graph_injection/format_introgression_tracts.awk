BEGIN {
    OFS="\t";
    anc[1]="Denisovan"; anc[2]="Neanderthal"; anc[3]="Mosaic";
}
FNR==NR {
    valid_ids[$1] = 1;
    next;
}
{
    current_id = $10;
    sub(/_/, "#", current_id);
    if (current_id in valid_ids) {
        # Normalize chromosome labels to a single "chr" prefix.
        chrom = $1;
        sub(/^chr/, "", chrom); 
        print "chr"chrom, $2, $3, anc[$5], current_id;
    }
}
