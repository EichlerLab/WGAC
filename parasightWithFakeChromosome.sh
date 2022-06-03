mkdir -p globalViewFakeChromosome
cd globalViewFakeChromosome

ln -s ../data/xw.al
ln -s ../fastalength.log
ln -s ../chromosomes.txt
ln -s ../parasight751.pl
ln -s ../globalview10k.pst

./makeFakeChromosome.sh
run_convert_xw.al.py.sh
module load perl/5.14.2 && for_parasight_fake_chromosome.sh

#cd globalViewFakeChromosome

test -f global_view_5k_90.png -a -f global_view_5k_95.png -a -f global_view_5k_98.png -a -f global_view_10k_90.png -a -f global_view_10k_95.png -a -f global_view_10k_98.png -a -f global_view_20k_90.png -a -f global_view_20k_95.png -a -f global_view_20k_98.png && touch ../../parasightWithFakeChromosomeDone

