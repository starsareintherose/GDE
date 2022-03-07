#/bin/csh

mkdir bin

#echo "Making blast..."
#cd BLAST
#Install.sh
#cd ..

echo "Making clustal..."
cd CLUSTAL
make
cd ..

echo "Making core GDE editor"
cd CORE
install.csh
cd ..

echo "Making FASTA"
cd FASTA
install.csh
cd ..

echo "Making Harvard Genome Lab functions"
cd HGL_SRC
install.csh
cd ..

echo "Making looptool"
cd LOOPTOOL
make
cd ..

echo "Making PHYLIP"
cd PHYLIP
install.csh
cd ..

echo "Making ReadSeq"
cd READSEQ
install.csh
cd ..

echo "Making other support programs"
cd SUPPORT
make
cd ..

echo "Making Zuker MFOLD"
cd ZUKER
install.csh
cd ..
