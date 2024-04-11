echo -n "Step1: Downloading data.. "
wget -q https://raw.githubusercontent.com/Pas-Kapli/assets/master/aux/Gallotia.COI.fasta

noseqs=$(grep ">" Gallotia.COI.fasta | wc -l)

echo "This alignment contains $noseqs number of sequences"

echo "Step2: Aligning sequences with mafft.."

mafft --quiet --preservecase Gallotia.COI.fasta > Gallotia.COI.aln

echo "Step3: Running RAxML.."
mkdir raxml
cd raxml
raxml-ng --search1 --msa Gallotia.COI.aln --model GTR+G --prefix Gallotia.COI &> /dev/null

echo "Step4: Running Species delimitation with mptp:"
cd ../
mkdir delimit
cd delimit
mptp --minbr_auto ../Gallotia.COI.aln --tree_file ../raxml/Gallotia.COI.raxml.bestTree --output_file minbr > minbr
tmpmin=$(tail -n 1 minbr | grep -oP "[0-9]*\.[0-9]*")
mptp --ml --multi --minbr `echo $tmpmin` --tree_file ../raxml/Gallotia.COI.raxml.bestTree --output_file delimitation &> /dev/null
nosp=$(grep "Number of delimited species" delimitation.txt | cut -f2 -d":" | sed "s/ //g")
echo "Final result: The number of species is $nosp"
