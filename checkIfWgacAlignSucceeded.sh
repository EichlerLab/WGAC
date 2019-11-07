nAlignBothFiles=`ls data/align_both/*/* | wc -l`
szString=`wc -l data/both.parse.defugu.trim.fixed.trim.defrac`
# szString looks like:  121550 data/both.parse.defugu.trim.fixed.trim.defrac
aTokens=( $szString )
n_both_parse_defugu_trim_fixed_trim_degrac_lines=${aTokens[0] }
# don't count the header line

echo $n_both_parse_defugu_trim_fixed_trim_degrac_lines
n_both_parse_defugu_trim_fixed_trim_degrac_lines=$(( $n_both_parse_defugu_trim_fixed_trim_degrac_lines - 1 ))

if [ $nAlignBothFiles -eq $n_both_parse_defugu_trim_fixed_trim_degrac_lines ]
then
    echo "data/align_both/*/*  and data/both.parse.defugu.trim.fixed.trim.defrac are the same size"
    exit 0
else
    echo "data/align_both/*/* size $nAlignBothFiles"
    echo "data/both.parse.defugu.trim.fixed.trim.defrac size $n_both_parse_defugu_trim_fixed_trim_degrac_lines"
    echo "data/align_both/*/* and data/both.parse.defugu.trim.fixed.trim.defrac size $n_both_parse_defugu_trim_fixed_trim_degrac_lines are not the same size so qsub -N global_align must have failed to complete successfully"
    exit 1
fi
