set -vex

cd data
cut -f 1,2,3,4,5,6,7,8,23,28 oo.weild10kb.join.all.cull > xw.join.all.cull
perl ../xwalign.pl xw.join.all.cull > xw.al

perl ../writeGenomeLengHash.pl ../fastalength.log > length_hash
perl ../writeGenomeLengtab.pl ../fastalength.log > length_tab

cd ..


echo -e "seqname\tlength" >showseq.out

# changed DG Mar 31, 2022
#sed '/random/d;/chrUn/d' data/length_tab | grep "^chr" | sort -V >> showseq.out

cat fastalength.log | sed 1d | sed '$d' | awk '{print $2"\t"$3}' | ./filterByTokenValue.py --szFileOfLegalValues chromosomes.txt --n0BasedToken 0 | sort -V >> showseq.out



mkdir -p globalView
cd globalView

# 5kb
# 90%
../parasight751.pl -showseq  ../showseq.out -align ../data/xw.al \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>5000, -filter2_col=>16, -filter2_min=>0.90, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_5k_90');" -die;

# 95%
../parasight751.pl -showseq  ../showseq.out -align ../data/xw.al \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>5000, -filter2_col=>16, -filter2_min=>0.95, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_5k_95');" -die;


# 98%
../parasight751.pl -showseq  ../showseq.out -align ../data/xw.al \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>5000, -filter2_col=>16, -filter2_min=>0.98, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_5k_98');" -die;


# above 3 for 10kb
../parasight751.pl -showseq  ../showseq.out -align ../data/xw.al \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>10000, -filter2_col=>16, -filter2_min=>0.90, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_10k_90');" -die;
 
../parasight751.pl -showseq  ../showseq.out -align ../data/xw.al \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>10000, -filter2_col=>16, -filter2_min=>0.95, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_10k_95');" -die;
 
../parasight751.pl -showseq  ../showseq.out -align ../data/xw.al \
    -template ../globalview10k.pst \
    -option '-filterpre2_min=>10000, -filter2_col=>16, -filter2_min=>0.98, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' \
    -precode "&fitlongestline; &print_all(0,'parasight_10k_98');" -die;

# 20kb

../parasight751.pl -showseq  ../showseq.out -align ../data/xw.al -template ../globalview10k.pst -option '-filterpre2_min=>20000, -filter2_col=>16, -filter2_min=>0.90, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' -precode "&fitlongestline; &print_all(0,'parasight_20k_90');" -die;
 
../parasight751.pl -showseq  ../showseq.out -align ../data/xw.al -template ../globalview10k.pst -option '-filterpre2_min=>20000, -filter2_col=>16, -filter2_min=>0.95, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' -precode "&fitlongestline; &print_all(0,'parasight_20k_95');" -die;
 
../parasight751.pl -showseq  ../showseq.out -align ../data/xw.al -template ../globalview10k.pst -option '-filterpre2_min=>20000, -filter2_col=>16, -filter2_min=>0.98, -extra_label_on=>0, -seq_tick_label_fontsize => 28, -seq_label_fontsize => 28, -printer_page_orientation=>0' -precode "&fitlongestline; &print_all(0,'parasight_20k_98');" -die;


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


# added DG, March 31, 2022 since we never use general views
exit 0;


cd ..

# general views

perl ./writeNr5kextra.pl data/stats/mergeHit \
    | sed 's/black/red/g' > wgac5k.extra
 
perl ./writeArrangeUsngLoop.pl showseq.out > arrange.out


mkdir -p generalView
cd generalView

../parasight751_4graph.pl -extra ../wgac5k.extra \
    -template ../template.pst \
    -showseq ../showseq.out -arrangeseq  file:../arrange.out \
    -options "-seq_tick_bp=>10000, -extra_label_on=>0, -sub_on=>0, -printer_page_orientation=>0" \
    -precode "&fitlongestline; &print_all(0,'general-view');" -die;

gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=general_view.pdf *.ps

cd ..

mkdir -p blowups
cd blowups


../fasta_findNs.pl -i ../fastawholesplit -o gap.log
sed -i 's/.fa//' gap.log

head -n 1 gap.log > 5kgap.extra
awk '$3 - $2 > 5000' gap.log >> 5kgap.extra

# ../david_chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull \
#     5kgap.extra 1000 0.90 chr1 1 ../showseq.out
../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull \
    5kgap.extra 1000 0.90 chr1 1 ../showseq.out

# 5 kb
../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull 5kgap.extra 5000 0.90 all 1 ../showseq.out
 
../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull 5kgap.extra 5000 0.94 all 1 ../showseq.out
 
../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull 5kgap.extra 5000 0.98 all 1 ../showseq.out
 
../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull 5kgap.extra 5000 0.99 all 1 ../showseq.out

# 10 kb

../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull 5kgap.extra 10000 0.90 all 1 ../showseq.out
 
../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull 5kgap.extra 10000 0.94 all 1 ../showseq.out
 
../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull 5kgap.extra 10000 0.98 all 1 ../showseq.out
 
../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull 5kgap.extra 10000 0.99 all 1 ../showseq.out

# 20 kb

../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull 5kgap.extra 20000 0.90 all 1 ../showseq.out
 
../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull 5kgap.extra 20000 0.94 all 1 ../showseq.out
 
../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull 5kgap.extra 20000 0.98 all 1 ../showseq.out
 
../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull 5kgap.extra 20000 0.99 all 1 ../showseq.out

# 40 kb

../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull 5kgap.extra 40000 0.90 all 1 ../showseq.out
 
../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull 5kgap.extra 40000 0.94 all 1 ../showseq.out
 
../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull 5kgap.extra 40000 0.98 all 1 ../showseq.out
 
../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull 5kgap.extra 40000 0.99 all 1 ../showseq.out

../chromosome_blowups.pl ../data/oo.weild10kb.join.all.cull 5kgap.extra 100000 0.95 all 1 ../showseq.out



perl ../dirBlowupP2p.pl `pwd`
