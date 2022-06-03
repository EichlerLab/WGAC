set -vex

# cd data
# cut -f 1,2,3,4,5,6,7,8,23,28 oo.weild10kb.join.all.cull.just.chromosomes > xw.join.all.cull.just.chromosomes

# perl ../xwalign.pl xw.join.all.cull.just.chromosomes > xw.al.just.chromosomes

# perl ../writeGenomeLengHash.pl ../fastalength.log > length_hash
# perl ../writeGenomeLengtab.pl ../fastalength.log > length_tab

# cd ..


echo -e "seqname\tlength" >showseq_fake_chromosome.out

cat fastalength.log | sed 1d | sed '$d' | awk '{print $2"\t"$3}' | ../filterByTokenValue.py --szFileOfLegalValues chromosomes.txt --n0BasedToken 0 | sort -V >> showseq_fake_chromosome.out

# add the fake chromosome
 awk '{x += $2 } END{ print "UNK\t"x }' fakeChromosomeOffsets.txt >>showseq_fake_chromosome.out


# 5kb
# 90%
../parasight751.pl -showseq  showseq_fake_chromosome.out -align fake_chromosome_xw.al \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>5000, -filter2_col=>16, -filter2_min=>0.90, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_5k_90');" -die;

# 95%
../parasight751.pl -showseq  showseq_fake_chromosome.out -align fake_chromosome_xw.al \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>5000, -filter2_col=>16, -filter2_min=>0.95, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_5k_95');" -die;


# 98%
../parasight751.pl -showseq  showseq_fake_chromosome.out -align fake_chromosome_xw.al \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>5000, -filter2_col=>16, -filter2_min=>0.98, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_5k_98');" -die;


# above 3 for 10kb
../parasight751.pl -showseq  showseq_fake_chromosome.out -align fake_chromosome_xw.al \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>10000, -filter2_col=>16, -filter2_min=>0.90, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_10k_90');" -die;
 
../parasight751.pl -showseq  showseq_fake_chromosome.out -align fake_chromosome_xw.al \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>10000, -filter2_col=>16, -filter2_min=>0.95, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_10k_95');" -die;
 
../parasight751.pl -showseq  showseq_fake_chromosome.out -align fake_chromosome_xw.al \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>10000, -filter2_col=>16, -filter2_min=>0.98, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_10k_98');" -die;

# 20kb

../parasight751.pl -showseq  showseq_fake_chromosome.out -align fake_chromosome_xw.al -template ../globalview10k.pst -option '-filterpre2_min=>20000, -filter2_col=>16, -filter2_min=>0.90, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' -precode "&fitlongestline; &print_all(0,'parasight_20k_90');" -die;
 
../parasight751.pl -showseq  showseq_fake_chromosome.out -align fake_chromosome_xw.al -template ../globalview10k.pst -option '-filterpre2_min=>20000, -filter2_col=>16, -filter2_min=>0.95, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' -precode "&fitlongestline; &print_all(0,'parasight_20k_95');" -die;
 
../parasight751.pl -showseq  showseq_fake_chromosome.out -align fake_chromosome_xw.al -template ../globalview10k.pst -option '-filterpre2_min=>20000, -filter2_col=>16, -filter2_min=>0.98, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' -precode "&fitlongestline; &print_all(0,'parasight_20k_98');" -die;


gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=global_view_5k_90.pdf parasight_5k_90.01.01.ps
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=global_view_5k_95.pdf parasight_5k_95.01.01.ps
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=global_view_5k_98.pdf parasight_5k_98.01.01.ps


convert -density 300 -depth 8  -background white -flatten global_view_5k_90.{pdf,png}
convert -density 300 -depth 8  -background white -flatten global_view_5k_95.{pdf,png}
convert -density 300 -depth 8  -background white -flatten global_view_5k_98.{pdf,png}


gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=global_view_10k_90.pdf parasight_10k_90.01.01.ps
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=global_view_10k_95.pdf parasight_10k_95.01.01.ps
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=global_view_10k_98.pdf parasight_10k_98.01.01.ps


convert -density 300 -depth 8  -background white -flatten global_view_10k_90.{pdf,png}
convert -density 300 -depth 8  -background white -flatten global_view_10k_95.{pdf,png}
convert -density 300 -depth 8  -background white -flatten global_view_10k_98.{pdf,png}


gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=global_view_20k_90.pdf parasight_20k_90.01.01.ps
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=global_view_20k_95.pdf parasight_20k_95.01.01.ps
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=global_view_20k_98.pdf parasight_20k_98.01.01.ps



convert -density 300 -depth 8  -background white -flatten global_view_20k_90.{pdf,png}
convert -density 300 -depth 8  -background white -flatten global_view_20k_95.{pdf,png}
convert -density 300 -depth 8  -background white -flatten global_view_20k_98.{pdf,png}


