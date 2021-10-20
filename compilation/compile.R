saemixDir="/home/eco/work/saemix/saemixextension"

workDir="/home/eco/work/saemix/versions/saemix3.0"

rm -r $workDir/saemix
mkdir $workDir/saemix

cp -rp $saemixDir/R $workDir/saemix
cp -rp $saemixDir/data $workDir/saemix

cp -rp $saemixDir/CHANGES $workDir/saemix/
cp -rp $saemixDir/DESCRIPTION $workDir/saemix/

cd $workDir
R < $saemixDir/compilation/roxygenisePackage.R --no-save  



R CMD build saemix
R CMD check --as-cran --run-donttest saemix_3.0.tar.gz
