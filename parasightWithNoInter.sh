mkdir -p globalViewNoInter
cd globalViewNoInter
ln -sf ../globalViewFakeChromosome/fake_chromosome_xw.al
ln -sf ../globalViewFakeChromosome/showseq_fake_chromosome.out

cat fake_chromosome_xw.al | awk '{if ( $1 == $5 ) print }' >xw.al.no_inter

module load perl/5.14.2
set -vex

# 5kb
# 90%
../parasight751.pl -showseq  showseq_fake_chromosome.out -align xw.al.no_inter \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>5000, -filter2_col=>16, -filter2_min=>0.90, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_5k_90');" -die;

# 95%
../parasight751.pl -showseq  showseq_fake_chromosome.out -align xw.al.no_inter \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>5000, -filter2_col=>16, -filter2_min=>0.95, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_5k_95');" -die;


# 98%
../parasight751.pl -showseq  showseq_fake_chromosome.out -align xw.al.no_inter \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>5000, -filter2_col=>16, -filter2_min=>0.98, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_5k_98');" -die;


# above 3 for 10kb
../parasight751.pl -showseq  showseq_fake_chromosome.out -align xw.al.no_inter \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>10000, -filter2_col=>16, -filter2_min=>0.90, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_10k_90');" -die;
 
../parasight751.pl -showseq  showseq_fake_chromosome.out -align xw.al.no_inter \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>10000, -filter2_col=>16, -filter2_min=>0.95, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_10k_95');" -die;
 
../parasight751.pl -showseq  showseq_fake_chromosome.out -align xw.al.no_inter \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>10000, -filter2_col=>16, -filter2_min=>0.98, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_10k_98');" -die;

# 20kb

../parasight751.pl -showseq  showseq_fake_chromosome.out -align xw.al.no_inter -template ../globalview10k.pst -option '-filterpre2_min=>20000, -filter2_col=>16, -filter2_min=>0.90, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' -precode "&fitlongestline; &print_all(0,'parasight_20k_90');" -die;
 
../parasight751.pl -showseq  showseq_fake_chromosome.out -align xw.al.no_inter -template ../globalview10k.pst -option '-filterpre2_min=>20000, -filter2_col=>16, -filter2_min=>0.95, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' -precode "&fitlongestline; &print_all(0,'parasight_20k_95');" -die;
 
../parasight751.pl -showseq  showseq_fake_chromosome.out -align xw.al.no_inter -template ../globalview10k.pst -option '-filterpre2_min=>20000, -filter2_col=>16, -filter2_min=>0.98, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' -precode "&fitlongestline; &print_all(0,'parasight_20k_98');" -die;


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

test -f global_view_5k_90.png -a -f global_view_5k_95.png -a -f global_view_5k_98.png -a -f global_view_10k_90.png -a -f global_view_10k_95.png -a -f global_view_10k_98.png -a -f global_view_20k_90.png -a -f global_view_20k_95.png -a -f global_view_20k_98.png && touch ../parasightWithNoInterDone



cd ..
