# First get all names
PathSeg=/home/csudre/NiftySeg_0.9.4/build_comic2/seg-apps
#PathSeg=/Users/csudre/Development/NiftySeg_build_Debug/seg-apps

if [ $# -lt 5 ]
then
echo ""
echo "******************************"
echo "Usage: sh ~/Scripts/LaplaceLobesScripts.sh Path ID PathGIF PathLesion PathT1"
echo "******************************"
echo ""
exit
fi


Path=${1}
ID=${2}
PN=${ID}
PathGIF=${3}
Lesion=${4}
T1=${5}
Model=${6}
Seg=${7}


FrontalArray=(101 102 103 104 105 106 113 114 119 120 121 122 125 126 137 138 141 142 143 144 147 148 151 152 153 154 163 164 165 166 179 180 183 184 187 188 191 192 193 194 205 206)
FrontalArray=(105 106 147 148 137 138 179 180 191 192 143 144 163 164 165 166 205 206 125 126 141 142 187 188 153 154 183 184 151 152 193 194 119 120 113 114 121 122)


ParietalArray=(107 108  149 150 169 170 175 176 177 178 195 196 199 200)
ParietalArray=(169 170 175 176 195 196 199 200 107 108 177 178 149 150)

OccipitalArray=(109 110 115 116 129 130 135 136 145 146 157 158 161 162 197 198)
OccipitalArray=(115 116 109 110 135 136 161 162 197 198 129 130 145 146 157 158)

TemporalArray=(117 118 123 124 133 134 155 156 181 182 185 186 201 202 203 204 207 208)
TemporalArray=(117 118 123 124 171 172 133 134 155 156 201 202 203 204 181 182 185 186 207 208)

BGArray=(16 24 31 32 33 37 38 56 57 58 59 60 61 76 77)
ITArray=(35 36 39 40 41 42 43 44 72 73 74)
VentrArray=(5 16 12 47 50 51 52 53)
ITArray=(35 36 39 40 41 42 43 44 72 73 74)
NameScript=${Path}/${ID}_scriptLayers.sh
echo "" > ${NameScript}

#if ((${#Model}<2))
#then
#read -a ValueModelText <<< $(ls ${Path}/T1*BiASM*${OptCross}.txt)
#echo "${Path}/T1*BiASM*${OptCross}*"
#read -a ValueModelNii <<<$(ls ${Path}/T1*BiASM*${OptCross}.nii.gz )
#else
#ValueModelText=(${Path}/${Model})
#ValueModelNii=(${Path}/${Seg})
#fi

read -a ValueLesion <<< $(ls ${Lesion}*)
# Get segmentation from GIF

stringAddition="${Path}/${ID}_ITConstruction_0.nii.gz "
for ((k=0;k<${#ITArray[@]};k++))
do
Value=${ITArray[k]}
echo "${PathSeg}/seg_maths ${PathGIF}/*Parcellation.nii.gz -equal $Value -bin ${Path}/${ID}_ITConstruction_${k}.nii.gz " >> ${NameScript}
stringAddition="${stringAddition} -add ${Path}/${ID}_ITConstruction_${k}.nii.gz "
done
echo "${PathSeg}/seg_maths ${stringAddition} -bin -odt char ${Path}/${ID}_Infratentorial.nii.gz " >> ${NameScript}
echo "rm ${Path}/*ITConstruction* " >> ${NameScript}


if [ ! -f ${Path}/Laplace*WMDGMLes*_4_${ID}_* ]
then
echo "${PathSeg}/seg_maths ${PathGIF}/*Segmentation* -tpmax ${Path}/${ID}_Seg1.nii.gz" >> ${NameScript}
echo "${PathSeg}/seg_maths ${Path}/${ID}_Seg1.nii.gz -thr 0.5 -uthr 1.5 ${Path}/${ID}_CSF_bin.nii.gz" >> ${NameScript}
echo "${PathSeg}/seg_maths ${Path}/${ID}_Seg1.nii.gz -thr 1.5 -uthr 2.5 ${Path}/${ID}_CGM_bin.nii.gz" >> ${NameScript}
echo "${PathSeg}/seg_maths ${Path}/${ID}_Seg1.nii.gz -thr 2.5 -uthr 3.5 ${Path}/${ID}_WM_bin.nii.gz" >> ${NameScript}
echo "${PathSeg}/seg_maths ${Path}/${ID}_Seg1.nii.gz -thr 3.5 -uthr 4.5 ${Path}/${ID}_DGM_bin.nii.gz" >> ${NameScript}
echo "${PathSeg}/seg_maths ${Path}/${ID}_Seg1.nii.gz -thr 4.5 -uthr 5.5 ${Path}/${ID}_Brainstem_bin.nii.gz" >> ${NameScript}
echo "${PathSeg}/seg_maths ${Path}/${ID}_DGM_bin.nii.gz -add ${Path}/${ID}_WM_bin.nii.gz ${Path}/${ID}_WMDGM_bin.nii.gz" >> ${NameScript}
echo "${PathSeg}/seg_maths ${Path}/${ID}_WMDGM_bin.nii.gz -smo 1 ${Path}/${ID}_WMDGM_ext.nii.gz" >> ${NameScript}
echo "${PathSeg}/seg_maths ${Path}/${ID}_CGM_bin.nii.gz -smo 1 ${Path}/${ID}_CGM_ext.nii.gz" >> ${NameScript}
# Creation of Ventr/Sim segmentation
echo "${PathSeg}/seg_maths ${PathGIF}/*Parcellation** -thr 49.5 -uthr 53.5 -bin ${Path}/${ID}_Ventr2.nii.gz" >> ${NameScript}
echo "${PathSeg}/seg_maths ${PathGIF}/*Parcellation* -thr 4.5 -uthr 16.5 -bin ${Path}/${ID}_Ventr3.nii.gz" >> ${NameScript}
echo "${PathSeg}/seg_maths ${Path}/${ID}_Ventr2.nii.gz -add ${Path}/${ID}_Ventr3.nii.gz -sub ${ValueLesion[0]} -thr 0.5 ${Path}/${ID}_VentrTot.nii.gz" >> ${NameScript}
# Creation of Parenchymal filled
echo "${PathSeg}/seg_maths ${Path}/${ID}_WMDGM_bin.nii.gz -add ${ValueLesion[0]} -add ${Path}/${ID}_Brainstem_bin.nii.gz -add ${Path}/${ID}_VentrTot.nii.gz -bin  -lconcomp ${Path}/${ID}_WMDGMLes.nii.gz " >> ${NameScript}
echo "${PathSeg}/Seg_Analysis -LapIO ${Path}/${ID}_VentrTot.nii.gz ${Path}/${ID}_WMDGMLes.nii.gz  -numbLaminae 1  4  -nameLap WMDGMLes" >> ${NameScript}
echo "cp ${Path}/LaplaceLayers*${PN}* ${Path}/Layers_${PN}.nii.gz" >> ${NameScript}
fi
echo "rm ${Path}/${ID}*bin.nii.gz ${Path}/${ID}*ext.nii.gz ${Path}/${ID}*Ventr*.nii.gz" >> ${NameScript}

PathScratch=/scratch0/${PN}_Generic${Opt}Lobes
echo "mkdir /scratch0/${PN}_Generic${Opt}Lobes ">> ${NameScript}
ChangePathString="-inChangePath ${PathScratch}/"
echo "${PathSeg}/seg_maths ${T1} -odt float ${PathScratch}/T1_${PN}.nii.gz " >> ${NameScript}
echo "cp ${Path}/Lap* ${PathScratch}/." >> ${NameScript}
echo "cp ${Path}/Brainmask/* ${PathScratch}/." >> ${NameScript}
#echo "cp ${Path}/*${Opt}.* ${PathScratch}/." >> ${NameScript}
#echo "cp ${Path}/*${Opt}_3.* ${PathScratch}/." >> ${NameScript}
echo "cp ${PathGIF}/* ${PathScratch}/." >> ${NameScript}
echo "cp ${PathGIF}/*Artefacts* ${PathScratch}/." >> ${NameScript}
echo "cp ${PathScratch}/*NeuroMorph*Parcellation* ${PathScratch}/${PN}_Parcellation.nii.gz" >> ${NameScript}

echo "${PathSeg}/seg_maths ${PathScratch}/*TIV.nii.gz -thr 0.5 -bin -odt char ${PathScratch}/${PN}_GIF_B1.nii.gz" >> ${NameScript}
read -a ListParc <<< $(ls ${PathGIF}/*Parcellation*)
echo "cp ${ListParc[0]} ${PathScratch}/${PN}_Parcellation.nii.gz" >> ${NameScript}
read -a ListTIV <<< $(ls ${PathGIF}/*TIV*)
echo ${ListTIV[*]}
echo ${ListParc[*]}

echo "cp ${ListTIV[0]} ${PathScratch}/${PN}_TIV.nii.gz" >> ${NameScript}
echo "${PathSeg}/seg_maths ${ListTIV[0]} -bin -odt char ${PathScratch}/${PN}_GIF_B1.nii.gz" >> ${NameScript}
 

stringAddition1="${PathScratch}/${PN}_FrontalConstruction_0.nii.gz "
stringAddition2="${PathScratch}/${PN}_FrontalConstruction_1.nii.gz "
for ((i=0;i<${#FrontalArray[@]};i++))
do
Value=${FrontalArray[i]}
echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_Parcellation.nii.gz -equal $Value  -bin ${PathScratch}/${PN}_FrontalConstruction_${i}.nii.gz " >> ${NameScript}
stringAddition1="${stringAddition1} -add ${PathScratch}/${PN}_FrontalConstruction_${i}.nii.gz "
((i++))
Value=${FrontalArray[i]}
echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_Parcellation.nii.gz -equal $Value -bin ${PathScratch}/${PN}_FrontalConstruction_${i}.nii.gz " >> ${NameScript}
stringAddition2="${stringAddition2} -add ${PathScratch}/${PN}_FrontalConstruction_${i}.nii.gz "
#((i++))
done
echo "${PathSeg}/seg_maths ${stringAddition2} -bin -odt char ${PathScratch}/${PN}_FrontalLobeL.nii.gz " >> ${NameScript}
echo "${PathSeg}/seg_maths ${stringAddition1} -bin -odt char ${PathScratch}/${PN}_FrontalLobeR.nii.gz " >> ${NameScript}
echo "rm ${PathScratch}/*FrontalConstruction* " >> ${NameScript}
echo "cp ${PathScratch}/*FrontalLobe*.nii.gz ${PathGIF}/." >> ${NameScript}

stringAddition1="${PathScratch}/${PN}_ParietalConstruction_0.nii.gz "
stringAddition2="${PathScratch}/${PN}_ParietalConstruction_1.nii.gz "
for ((i=0;i<${#ParietalArray[@]};i++))
do
Value=${ParietalArray[i]}
echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_Parcellation.nii.gz -equal $Value  -bin ${PathScratch}/${PN}_ParietalConstruction_${i}.nii.gz " >> ${NameScript}
stringAddition1="${stringAddition1} -add ${PathScratch}/${PN}_ParietalConstruction_${i}.nii.gz "
((i++))
Value=${ParietalArray[i]}
echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_Parcellation.nii.gz -equal $Value  -bin ${PathScratch}/${PN}_ParietalConstruction_${i}.nii.gz " >> ${NameScript}
stringAddition2="${stringAddition2} -add ${PathScratch}/${PN}_ParietalConstruction_${i}.nii.gz "

done
echo "${PathSeg}/seg_maths ${stringAddition2} -bin -odt char ${PathScratch}/${PN}_ParietalLobeL.nii.gz " >> ${NameScript}
echo "${PathSeg}/seg_maths ${stringAddition1} -bin -odt char ${PathScratch}/${PN}_ParietalLobeR.nii.gz " >> ${NameScript}
echo "rm ${PathScratch}/*ParietalConstruction* " >> ${NameScript}
echo "cp ${PathScratch}/*ParietalLobe*.nii.gz ${PathGIF}/." >> ${NameScript}

stringAddition1="${PathScratch}/${PN}_OccipitalConstruction_0.nii.gz "
stringAddition2="${PathScratch}/${PN}_OccipitalConstruction_1.nii.gz "
for ((i=0;i<${#OccipitalArray[@]};i++))
do
Value=${OccipitalArray[i]}
echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_Parcellation.nii.gz -equal $Value  -bin ${PathScratch}/${PN}_OccipitalConstruction_${i}.nii.gz " >> ${NameScript}
stringAddition1="${stringAddition1} -add ${PathScratch}/${PN}_OccipitalConstruction_${i}.nii.gz "

((i++))
Value=${OccipitalArray[i]}
echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_Parcellation.nii.gz -equal $Value  -bin ${PathScratch}/${PN}_OccipitalConstruction_${i}.nii.gz " >> ${NameScript}
stringAddition2="${stringAddition2} -add ${PathScratch}/${PN}_OccipitalConstruction_${i}.nii.gz "

done
echo "${PathSeg}/seg_maths ${stringAddition2} -bin -odt char ${PathScratch}/${PN}_OccipitalLobeL.nii.gz " >> ${NameScript}
echo "${PathSeg}/seg_maths ${stringAddition1} -bin -odt char ${PathScratch}/${PN}_OccipitalLobeR.nii.gz " >> ${NameScript}
echo "rm ${PathScratch}/*OccipitalConstruction* " >> ${NameScript}
echo "cp ${PathScratch}/*OccipitalLobe*.nii.gz ${PathGIF}/." >> ${NameScript}

stringAddition1="${PathScratch}/${PN}_TemporalConstruction_0.nii.gz "
stringAddition2="${PathScratch}/${PN}_TemporalConstruction_1.nii.gz "
for ((i=0;i<${#TemporalArray[@]};i++))
do
Value=${TemporalArray[i]}
echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_Parcellation.nii.gz -equal $Value  -bin ${PathScratch}/${PN}_TemporalConstruction_${i}.nii.gz " >> ${NameScript}
stringAddition1="${stringAddition1} -add ${PathScratch}/${PN}_TemporalConstruction_${i}.nii.gz "

((i++))
Value=${TemporalArray[i]}
echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_Parcellation.nii.gz -equal $Value  -bin ${PathScratch}/${PN}_TemporalConstruction_${i}.nii.gz " >> ${NameScript}
stringAddition2="${stringAddition2} -add ${PathScratch}/${PN}_TemporalConstruction_${i}.nii.gz "
done
echo "${PathSeg}/seg_maths ${stringAddition2} -bin -odt char ${PathScratch}/${PN}_TemporalLobeL.nii.gz " >> ${NameScript}
echo "${PathSeg}/seg_maths ${stringAddition1} -bin -odt char ${PathScratch}/${PN}_TemporalLobeR.nii.gz " >> ${NameScript}
echo "rm ${PathScratch}/*TemporalConstruction* " >> ${NameScript}
echo "cp ${PathScratch}/*TemporalLobe*.nii.gz ${PathGIF}/." >> ${NameScript}
stringAddition="${PathScratch}/${PN}_BGConstruction_0.nii.gz "
for ((i=0;i<${#BGArray[@]};i++))
do
Value=${BGArray[i]}
echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_Parcellation.nii.gz -equal $Value -bin -ero 1 -lconcomp -dil 1 ${PathScratch}/${PN}_BGConstruction_${i}.nii.gz " >> ${NameScript}
stringAddition="${stringAddition} -add ${PathScratch}/${PN}_BGConstruction_${i}.nii.gz "
done
echo "${PathSeg}/seg_maths ${stringAddition} -bin -dil 4 -fill -ero 4 -odt char ${PathScratch}/${PN}_BasalGanglia.nii.gz " >> ${NameScript}
echo "rm ${PathScratch}/*BGConstruction* " >> ${NameScript}

echo "cp ${PathScratch}/*BasalGanglia.nii.gz ${PathGIF}/." >> ${NameScript}

stringAddition="${PathScratch}/${PN}_ITConstruction_0.nii.gz "
for ((i=0;i<${#ITArray[@]};i++))
do
Value=${ITArray[i]}

echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_Parcellation.nii.gz -equal $Value  -bin ${PathScratch}/${PN}_ITConstruction_${i}.nii.gz " >> ${NameScript}
stringAddition="${stringAddition} -add ${PathScratch}/${PN}_ITConstruction_${i}.nii.gz "
done
echo "${PathSeg}/seg_maths ${stringAddition} -bin -odt char ${PathScratch}/${PN}_Infratentorial.nii.gz " >> ${NameScript}
echo "rm ${PathScratch}/*ITConstruction* " >> ${NameScript}
echo "cp ${PathScratch}/*Infratentorial.nii.gz ${PathGIF}/." >> ${NameScript}


stringAddition="${PathScratch}/${PN}_VentrConstruction_0.nii.gz "
for ((i=0;i<${#VentrArray[@]};i++))
do
Value=${VentrArray[i]}

echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_Parcellation.nii.gz -equal $Value  -bin ${PathScratch}/${PN}_VentrConstruction_${i}.nii.gz " >> ${NameScript}
stringAddition="${stringAddition} -add ${PathScratch}/${PN}_VentrConstruction_${i}.nii.gz "
done
echo "${PathSeg}/seg_maths ${stringAddition} -bin -odt char ${PathScratch}/${PN}_Ventricles.nii.gz " >> ${NameScript}
echo "rm ${PathScratch}/*VentrConstruction* " >> ${NameScript}
echo "cp ${PathScratch}/*Ventricle*.nii.gz ${PathGIF}/." >> ${NameScript}

LaplaceSolFile=${Path}/LaplaceSol_${ID}_VentrTot.nii.gz
LesionArray=$(ls ${Path}/LesionCorr* )
echo "${PathSeg}/seg_maths ${PathScratch}/T1_${PN}.nii.gz -odt float ${PathScratch}/T1_${PN}.nii.gz" >> ${NameScript}
echo "${PathSeg}/Seg_Analysis -LaplaceNormSol ${LaplaceSolFile} -inLobes 10 ${PathScratch}/${PN}_FrontalLobeL.nii.gz ${PathScratch}/${PN}_FrontalLobeR.nii.gz ${PathScratch}/${PN}_ParietalLobeL.nii.gz ${PathScratch}/${PN}_ParietalLobeR.nii.gz ${PathScratch}/${PN}_OccipitalLobeL.nii.gz ${PathScratch}/${PN}_OccipitalLobeR.nii.gz ${PathScratch}/${PN}_TemporalLobeL.nii.gz ${PathScratch}/${PN}_TemporalLobeR.nii.gz ${PathScratch}/${PN}_BasalGanglia.nii.gz ${PathScratch}/${PN}_Infratentorial.nii.gz -mask ${PathScratch}/${PN}_GIF_B1.nii.gz -inLesCorr ${PathScratch}/T1_${PN}.nii.gz -inVentricleSeg ${PathScratch}/${PN}_Ventricles.nii.gz -ParcellationIn ${PathScratch}/${PN}_Parcellation.nii.gz" >> ${NameScript}
echo "${PathSeg}"
echo "cp ${PathScratch}/DistanceChoice_T1_${PN}.nii.gz ${Path}/Lobes_${PN}.nii.gz" >> ${NameScript}

echo "rm -r ${PathScratch} " >> ${NameScript}

if [[ $2 == "sh" ]]
then
sh ${NameScript}
else
qsub -S /bin/bash  -l h_rt=1:0:0 -l h_vmem=3.9G -l tmem=3.9G -l tscratch=20G -N  BaMoSRec${PN} ${NameScript}
fi
