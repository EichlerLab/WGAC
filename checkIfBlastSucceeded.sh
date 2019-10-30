nFuguSize=`ls fugu2/* | wc -l`
nBlastoutSize=`ls blastout/* | wc -l`

if [ $nFuguSize -eq $nBlastoutSize ]
then
    echo "fugu and blastout are the same size"
    exit 0
else
    echo "fugu and blastout are not the same size"
    exit 1
fi
