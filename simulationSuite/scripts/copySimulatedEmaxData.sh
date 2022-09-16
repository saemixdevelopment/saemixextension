% Copy Hill and Emax simulated data (K=100 first simulated datasets)

fromDir='/home/eco/work/saemix/newLib/zesims/simData/pdhillhigh.rich'
toDir='/home/eco/work/saemix/saemixextension/simulationSuite/cont/hillRichProp/data'

mkdir $toDir
cp $fromDir/data_pdhillhigh?.tab $toDir/
cp $fromDir/data_pdhillhigh??.tab $toDir/
cp $fromDir/data_pdhillhigh100.tab $toDir/
cp $fromDir/param_pdhillhigh?.tab $toDir/
cp $fromDir/param_pdhillhigh??.tab $toDir/
cp $fromDir/param_pdhillhigh100.tab $toDir/


fromDir='/home/eco/work/saemix/newLib/zesims/simData/pdhillhigh.sparse'
toDir='/home/eco/work/saemix/saemixextension/simulationSuite/cont/hillSparseProp/data'

mkdir $toDir
cp $fromDir/data_pdhillhigh?.tab $toDir/
cp $fromDir/data_pdhillhigh??.tab $toDir/
cp $fromDir/data_pdhillhigh100.tab $toDir/
cp $fromDir/param_pdhillhigh?.tab $toDir/
cp $fromDir/param_pdhillhigh??.tab $toDir/
cp $fromDir/param_pdhillhigh100.tab $toDir/


fromDir='/home/eco/work/saemix/newLib/zesims/simData/pdemax.sparse'
toDir='/home/eco/work/saemix/saemixextension/simulationSuite/cont/emaxSparseProp/data'

mkdir $toDir
cp $fromDir/data_pdemax?.tab $toDir/
cp $fromDir/data_pdemax??.tab $toDir/
cp $fromDir/data_pdemax100.tab $toDir/
cp $fromDir/param_pdemax?.tab $toDir/
cp $fromDir/param_pdemax??.tab $toDir/
cp $fromDir/param_pdemax100.tab $toDir/

fromDir='/home/eco/work/saemix/newLib/zesims/simData/pdemax.rich'
toDir='/home/eco/work/saemix/saemixextension/simulationSuite/cont/emaxRichProp/data'

mkdir $toDir
cp $fromDir/data_pdemax?.tab $toDir/
cp $fromDir/data_pdemax??.tab $toDir/
cp $fromDir/data_pdemax100.tab $toDir/
cp $fromDir/param_pdemax?.tab $toDir/
cp $fromDir/param_pdemax??.tab $toDir/
cp $fromDir/param_pdemax100.tab $toDir/
