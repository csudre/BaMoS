#!/bin/bash
# set -x

#### WARNING !!!! THIS SCRIPT ASSUMES THAT THE FLAIR IMAGE WILL BE REGISTERED TO THE T1 IMAGE. ALL RESULTS SHOULD BE VISUALISED IN THE T1 SPACE. ALL SCANS ARE PUT IN THE RAS ORIENTATION BY DEFAULT

ArrayModalities=("T1" "T2" "FLAIR" "PD" "SWI")
ListCorr=(102 103 187 188 47 48 49)
AW=0
JC=0
MRE=250
OptSP=0
OptCL=2
OptOL=0
OptWMI=0 # If >0, indicates that a higher level of sensitivity to lesion should be considered from voxels coming from WMI
TVS=1
arrayMod=("T1" "FLAIR" )
arrayModNumber=(1 3)
ArtefactArray=(101 102 139 140 187 188 167 168) 
ArtefactArray=(5 12 32 33 39 40 50 51 62 63 64 65 72 73 74 48 49 101 102 117 118 139 140 171 172 187 188 167 168 202 203 181 182 117 118 123 124 133 134 155 156 181 182 185 186 201 202 203 204 207 208)
ListCorr=(5 12 32 33 39 40 50 51 62 63 64 65 72 73 74 48 49 101 102 117 118 139 140 171 172 187 188 167 168 202 203 181 182 117 118 123 124 133 134 155 156 181 182 185 186 201 202 203 204 207 208)
#ArtefactArray=(101 102 139 140 187 188 167 168 39 40 117 118 123 124 125 126 133 134 137 138 141 142 147 148 155 156 181 182 185 186 187 188 201 202 203 204 207 208)

RuleFileName=/mounts/auto/xnat/pipelines/BiASM/GenericRule_CSF.txt

NameGMatrix=/mounts/auto/xnat/pipelines/BiASM/GMatrix4_Low3.txt

OW=0.01
#OW=0.01
#OW=0.10
#OW=0.50
#OW=0.001
#OW=0.99
JC=1
Opt=${OptCross}
OptVL=1
OptSP=1
OptTA=1
# Flags for levels of correction
flag_NoCerrebellum=0
flag_furtherCorr=1
flag_corrDGM=1
 
# PathReg=/share/apps/cmic/niftyreg-10.4.15/bin
# PathSeg=/home/csudre/NiftySeg_0.9.4/build_comic2/seg-apps
# PathFSL=$FSLDIR/bin
# PathICBM=/cluster/project0/SegBiASM/ICBM_Priors
PathReg=/Users/csudre/Development/niftyreg_install//bin/
PathSeg=/Users/csudre/Development/NiftySeg_build_Debug//seg-apps/
PathFSL=$FSLDIR/bin
PathICBM=/Users/csudre/Documents/Priors/ICBM_Priors


RuleFileName=/cluster/project0/SegBiASM/DataToTryBaMoS/GenericRule_CSF.txt
NameGMatrix=/cluster/project0/SegBiASM/DataToTryBaMoS/GMatrix4_Low3.txt
NameGMatrix=/Users/csudre/Dropbox/Scripts/GMatrix4_Low3.txt
RuleFileName=/Users/csudre/Dropbox/Scripts/GenericRule_CSF.txt

if [ $# -lt 5 ] || [ $# -gt 6 ]
then
	echo ""
	echo "*******************************************************************************"
	echo "Usage: $0 ID ImageFLAIR ImageT1 GIF_results_path PathResults JUMP_START"
	echo "*******************************************************************************"
	echo ""
	exit
fi

ID=$1
ImageFLAIR=$2
ImageT1=$3
GIF_results_path=$4
PathResults=$5
JUMP_START=$6
ModalitiesTot=T1FLAIR
echo ${PathResults}

if [ ! -d $PathResults ]
then
    mkdir $PathResults
fi

cd $PathResults

Space=1 #(FLAIR space)

ChangePathString="-inChangePath ${PathResults}/"
StringPriorsIO=" "

if [ -z "$JUMP_START" ] || [ "$JUMP_START" -eq 0 ]; then
    echo "Reorientation to standard space of T1 image"
    ${PathFSL}/fslreorient2std ${ImageT1} ${PathResults}/T1_${ID}.nii.gz
    ${PathSeg}/seg_maths ${PathResults}/T1_${ID}.nii.gz -odt float ${PathResults}/T1_${ID}.nii.gz

    echo "Alignment of FLAIR image to T1"
    ${PathReg}/reg_aladin -ref ${PathResults}/T1_${ID}.nii.gz -flo ${ImageFLAIR} -aff ${PathResults}/Aff_FLAIRtoT1.txt -res ${PathResults}/FLAIR_${ID}.nii.gz -rigOnly
    ${PathSeg}/seg_maths ${PathResults}/FLAIR_${ID}.nii.gz -odt float ${PathResults}/FLAIR_${ID}.nii.gz

    # echo "Alignment of T2 image to FLAIR"
    # ${PathReg}/reg_aladin -ref ${PathResults}/FLAIR_${ID}.nii.gz -flo ${ImageT2} -aff ${PathResults}/Aff_T2toFLAIR.txt -res ${PathResults}/T2_${ID}.nii.gz
    # ${PathSeg}/seg_maths ${PathResults}/T2_${ID}.nii.gz -odt float ${PathResults}/T2_${ID}.nii.gz


    echo "Reorientation of T1 TIV mask and binarisation"
    ${PathFSL}/fslreorient2std ${GIF_results_path}/*TIV* ${PathResults}/Mask_${ID}.nii.gz
    ${PathSeg}/seg_maths ${PathResults}/Mask_${ID}.nii.gz -bin -odt char ${PathResults}/GIF_${ID}_B1.nii.gz

    echo "Reorientation of GIF priors and creation of GIF priors array"
    ${PathFSL}/fslreorient2std ${GIF_results_path}/*prior* ${PathResults}/GIF_prior_${ID}.nii.gz

#  If GIF not performed in FLAIR space, need to resample parcellation to FLAIR space
    echo "Reorientation of GIF parcellation"
    ${PathFSL}/fslreorient2std ${GIF_results_path}/*Parcellation* ${PathResults}/GIF_Parcellation_${ID}.nii.gz 

echo "Reorientation of GIF segmentation"
    ${PathFSL}/fslreorient2std ${GIF_results_path}/*Segmentation* ${PathResults}/GIF_Segmentation_${ID}.nii.gz 

    PriorsArray=("Out" "CSF" "CGM" "WMI" "DGM" "Brainstem")
    array_Priors=()

    for ((p=0;p<6;p++))
    do
        ${PathSeg}/seg_maths ${PathResults}/GIF_prior* -tp ${p} ${PathResults}/GIF_${PriorsArray[p]}_${ID}.nii.gz
    done

    
    

    echo "Creation of artefact map based on T1 parcellation and registration to FLAIR"

    if ((${#ArtefactArray[@]}>0))
    then
        stringAddition="${PathResults}/${ID}_ArtConstruction_0.nii.gz "
        for ((i=0;i<${#ArtefactArray[@]};i++))
        do
            Value=${ArtefactArray[i]}
            ValueMin=`echo "$Value - 0.5" | bc -l`
            ValueMax=`echo "$Value + 0.5" | bc -l`
            ${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz -thr $ValueMin -uthr $ValueMax -bin ${PathResults}/${ID}_ArtConstruction_${i}.nii.gz
            stringAddition="${stringAddition} -add ${PathResults}/${ID}_ArtConstruction_${i}.nii.gz "
        done
        ${PathSeg}/seg_maths ${stringAddition} -bin ${PathResults}/${ID}_Artefacts.nii.gz
    fi

    rm ${PathResults}/*_ArtConstruction_*

    echo "Registration of ICBM template and creation of ICBM atlases"
    ${PathReg}/reg_aladin -ref ${PathResults}/T1_${ID}.nii.gz -flo ${PathICBM}/ICBM_Template.nii.gz -aff ${PathResults}/${ID}_AffTransf.txt
    ${PathReg}/reg_f3d -ref ${PathResults}/T1_${ID}.nii.gz -flo ${PathICBM}/ICBM_Template.nii.gz -aff ${PathResults}/${ID}_AffTransf.txt -cpp ${PathResults}/${ID}_cpp.nii.gz

    for p in CGM DGM ECSF ICSF Out WM
    do
        ${PathReg}/reg_resample -ref ${PathResults}/T1_${ID}.nii.gz -flo ${PathICBM}/ICBM_${p}.nii.gz -cpp ${PathResults}/${ID}_cpp.nii.gz -res ${PathResults}/ICBM_${p}_${ID}.nii.gz
    done

    ${PathSeg}/seg_maths ${PathResults}/ICBM_DGM_${ID}.nii.gz -add ${PathResults}/ICBM_CGM_${ID}.nii.gz ${PathResults}/ICBM_GM_${ID}.nii.gz

    ${PathSeg}/seg_maths ${PathResults}/ICBM_ICSF_${ID}.nii.gz -add ${PathResults}/ICBM_ECSF_${ID}.nii.gz ${PathResults}/ICBM_CSFs_${ID}.nii.gz

    ${PathSeg}/seg_maths ${PathResults}/ICBM_DGM_${ID}.nii.gz -bin -mul ${PathResults}/GIF_DGM_${ID}.nii.gz ${PathResults}/GIF_DGM_${ID}_bis.nii.gz
    
${PathSeg}/seg_maths ${PathResults}/GIF_CGM_${ID}.nii.gz -add ${PathResults}/GIF_DGM_${ID}_bis.nii.gz ${PathResults}/GIF_GM_${ID}.nii.gz
    ${PathSeg}/seg_maths ${PathResults}/GIF_DGM_${ID}.nii.gz -sub ${PathResults}/GIF_DGM_${ID}_bis.nii.gz -add ${PathResults}/GIF_WMI_${ID}.nii.gz ${PathResults}/GIF_WMI_${ID}_bis.nii.gz

    ${PathSeg}/seg_maths ${PathResults}/GIF_WMI_${ID}_bis.nii.gz -add ${PathResults}/GIF_Brainstem_${ID}.nii.gz ${PathResults}/GIF_WM_${ID}.nii.gz

    

    arrayImage=()
    arrayModNumber=()
    echo "Modalities to put are ${arrayMod[0]}"
    for ((m=0;m<${#arrayMod[@]};m++))
    do
        for ((pos=0;pos<${#ArrayModalities[@]};pos++))
        do
            TmpModa=${ArrayModalities[pos]}
            TmptestModa=${arrayMod[m]}
            Subtracted="${TmptestModa/$TmpModa}"
    #echo ${#Subtracted}
            if ((${#Subtracted}<${#TmptestModa}))
            then
    # echo "Testing between $TmptestModa $TmpModa"
                FinModa=$((pos+1))
                arrayModNumber=(${arrayModNumber[*]} $FinModa)
            fi
        done
        arrayImage=( "${arrayImage[*]}" "${PathResults}/${arrayMod[m]}_${ID}.nii.gz" )
    done

array_Priors=("${PathResults}/GIF_GM_${ID}.nii.gz" "${PathResults}/GIF_WM_${ID}.nii.gz" "${PathResults}/GIF_CSF_${ID}.nii.gz" "${PathResults}/GIF_Out_${ID}.nii.gz")
    echo "Segmentation SegBiASM"
    #${PathSeg}/Seg_BiASM -in 2 ${arrayImage[0]} ${arrayImage[1]} -priors 4 ${array_Priors[*]} -mask ${PathResults}/GIF_${ID}_B1.nii.gz -out 2 ${PathResults}/${arrayMod[0]}${arrayMod[1]}_BiASM_${ID}_${Opt}.nii.gz ${PathResults}/${arrayMod[0]}${arrayMod[1]}_BiASMG_${ID}.nii.gz -txt_out ${PathResults}/${arrayMod[0]}${arrayMod[1]}_BiASM_${ID}_${Opt}.txt -bc_order 3 -CovPriors 8 -BFP 1 -maxRunEM ${MRE} -AtlasSmoothing 1 1 -AtlasWeight 1 ${AW}  -SMOrder 0 -KernelSize 3 -PriorsKept 5 -VLkappa 1.5 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -deleteUW 1 -outliersM 3 -outliersW ${OW} -init_splitUnif 0 -splitAccept 0 -unifTot 1 -MRF 1 -GMRF ${NameGMatrix} -juxtaCorr ${JC} -progMod 0 -priorDGM ${PathResults}/ICBM_DGM_${ID}.nii.gz -TypicalityAtlas ${OptTA} ${StringPriorsIO}

echo"
    ${PathSeg}/Seg_BiASM -VLkappa 3 -in 2 ${arrayImage[*]} -priors 4 ${array_Priors[*]} -mask ${PathResults}/GIF_${ID}_B1.nii.gz -out 2 ${PathResults}/${ModalitiesTot}_BiASM_${ID}_${Opt}.nii.gz ${PathResults}/${ModalitiesTot}_BiASMG_${ID}.nii.gz -txt_out ${PathResults}/${ModalitiesTot}_BiASM_${ID}_${Opt}.txt -bc_order 3 -CovPriors 8 -BFP 1 -maxRunEM ${MRE} -AtlasSmoothing 1 1 -AtlasWeight 1 ${AW}  -SMOrder 0 -KernelSize 3 -PriorsKept 5 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -deleteUW 1 -outliersM 3 -outliersW ${OW} -init_splitUnif 0 -splitAccept 0 -unifTot 1 -MRF 1 -GMRF ${NameGMatrix} -juxtaCorr ${JC} -progMod 0 -priorDGM ${PathResults}/ICBM_DGM_${ID}.nii.gz -TypicalityAtlas ${OptTA} "



${PathSeg}/Seg_BiASM -VLkappa 3 -in 2 ${arrayImage[*]} -priors 4 ${array_Priors[*]} -mask ${PathResults}/GIF_${ID}_B1.nii.gz -out 2 ${PathResults}/${ModalitiesTot}_BiASM_${ID}_${Opt}.nii.gz ${PathResults}/${ModalitiesTot}_BiASMG_${ID}.nii.gz -txt_out ${PathResults}/${ModalitiesTot}_BiASM_${ID}_${Opt}.txt -bc_order 3 -CovPriors 8 -BFP 1 -maxRunEM ${MRE} -AtlasSmoothing 1 1 -AtlasWeight 1 ${AW}  -SMOrder 0 -KernelSize 3 -PriorsKept 5 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -deleteUW 1 -outliersM 3 -outliersW ${OW} -init_splitUnif 0 -splitAccept 0 -unifTot 1 -MRF 1 -GMRF ${NameGMatrix} -juxtaCorr ${JC} -progMod 0 -priorDGM ${PathResults}/ICBM_DGM_${ID}.nii.gz -TypicalityAtlas ${OptTA}

rm ${PathResults}/BG* ${PathResults}/MRF*

    echo "Lesion segmentation"
    ${PathSeg}/Seg_Analysis  -inTxt2 ${PathResults}/${ModalitiesTot}_BiASM_${ID}_${Opt}.txt ${PathResults}/${ModalitiesTot}_BiASM_${ID}_${Opt}.nii.gz -mask ${PathResults}/GIF_${ID}_B1.nii.gz -Package 1 -SegType 1 -WeightedSeg 3 3 1 -connect -correct -inModa 2 1 3 -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${PathResults}/ICBM_ICSF_${ID}.nii.gz -inPriorsDGM ${PathResults}/ICBM_DGM_${ID}.nii.gz -inPriorsCGM ${PathResults}/ICBM_CGM_${ID}.nii.gz -inPriorsECSF ${PathResults}/ICBM_ECSF_${ID}.nii.gz -TO 1 -ParcellationIn ${PathResults}/GIF_Parcellation_${ID}.nii.gz -typeVentrSeg ${TVS} -Simple -Secondary 60 -juxtaCorr ${JC}  -SP ${OptSP} -LevelCorrection ${OptCL} ${ChangePathString} -LesWMI ${OptWMI} -Neigh 18

fi
echo "cp ${PathResults}/LesionCorrected*WS3WT3WC1* ${PathResults}/PrimaryLesions_${ID}.nii.gz"
 echo  " cp ${PathResults}/SecondaryCorrected*WS3WT3WC1* ${PathResults}/SecondaryLesions_${ID}.nii.gz "
    cp ${PathResults}/LesionCorrected*WS3WT3WC1* ${PathResults}/PrimaryLesions_${ID}.nii.gz
    cp ${PathResults}/SecondaryCorrected*WS3WT3WC1* ${PathResults}/SecondaryLesions_${ID}.nii.gz

    rm ${PathResults}/DataR* ${PathResults}/LesionT* ${PathResults}/Summ*

    echo "Correction for septum pellucidum"
    ${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz -thr 65.5 -uthr 67.5 -bin ${PathResults}/VentricleLining.nii.gz 
    ${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz  -equal 52 -bin -euc -uthr 0 -abs ${PathResults}/Temp1.nii.gz
    ${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz  -equal 53 -bin -euc -uthr 0 -abs -add ${PathResults}/Temp1.nii.gz -uthr 5 -bin ${PathResults}/PotentialInterVentr.nii.gz

    ${PathSeg}/seg_maths ${PathResults}/VentricleLining.nii.gz -dil 1 ${PathResults}/ExpandedVentrLin.nii.gz
    ${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz  -thr 49.5 -uthr 53.5 -bin -sub ${PathResults}/ExpandedVentrLin.nii.gz -thr 0 -bin ${PathResults}/PotentialCP.nii.gz

    #${PathSeg}/seg_maths ${GIF_results_path}/*Parcellation*  -thr 122.5 -uthr 124.5 -bin ${PathResults}/Temp4.nii.gz

    ${PathSeg}/seg_maths ${PathResults}/GIF_Segmentation_${ID}.nii.gz -tp 4 -thr 0.5 -bin -sub ${PathResults}/VentricleLining* -thr 0 ${PathResults}/SegDGM.nii.gz
    ${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz  -equal 87 -mul -1 -add ${PathResults}/PotentialInterVentr.nii.gz -sub ${PathResults}/SegDGM.nii.gz -thr 0 ${PathResults}/PotentialSP3.nii.gz
    ${PathSeg}/seg_maths  ${PathResults}/PotentialSP3.nii.gz -add ${PathResults}/PotentialCP.nii.gz ${PathResults}/PotentialSPCP.nii.gz
    ${PathSeg}/seg_maths ${PathResults}/ICBM_ICSF* -thr 0.3 -bin -mul ${PathResults}/PotentialSPCP.nii.gz -add ${PathResults}/PotentialSP3.nii.gz -thr 0.2  -bin  ${PathResults}/PotentialSPCP2.nii.gz

    ${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz -thr 83.5 -uthr 84.5 -bin -dil 5 ${PathResults}/WMDil.nii.gz
    ${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz -thr 91.5 -uthr 92.5 -bin -dil 5 ${PathResults}/WMDil2.nii.gz
    ${PathSeg}/seg_maths ${PathResults}/PotentialSPCP2.nii.gz -sub ${PathResults}/WMDil.nii.gz -sub ${PathResults}/WMDil2.nii.gz -thr 0 ${PathResults}/PotentialSPCP3.nii.gz 

    ${PathSeg}/seg_maths ${PathResults}/PrimaryLesions_* -sub ${PathResults}/PotentialSPCP3.nii.gz -thr 0 ${PathResults}/PrimarySPCP_${ID}.nii.gz
    ${PathSeg}/seg_maths ${PathResults}/SecondaryLesions_* -sub ${PathResults}/PotentialSPCP3.nii.gz -thr 0 ${PathResults}/SecondarySPCP_${ID}.nii.gz


    echo "Refined correction for third ventricle"
    ${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz  -sub ${PathResults}/GIF_Parcellation_${ID}.nii.gz ${PathResults}/Artefacts_${ID}.nii.gz 

for ((i=0;i<${#ListCorr[@]};i++))
do
${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz  -equal ${ListCorr[i]} -bin -add ${PathResults}/Artefacts_${ID}.nii.gz 
done
${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz  -thr 65.5 -uthr 67.5 ${PathResults}/VentrLin.nii.gz 
${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz  -equal 87 -add ${PathResults}/VentrLin.nii.gz ${PathResults}/VentrLin.nii.gz 
${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz  -equal 47 -euc  -abs -uthr 5 -bin -mul ${PathResults}/VentrLin.nii.gz ${PathResults}/VentrLinSP.nii.gz 



# ${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz  -equal 52 -bin -euc  -abs  ${PathResults}//Ventr1.nii.gz 

# ${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz  -equal 53 -bin -euc -abs -uthr 5 -sub  ${PathResults}/Ventr1.nii.gz -thr -5 -uthr 5 -bin -mul  ${PathResults}/GIF_Parcellation_${ID}.nii.gz  -equal 52 -bin -euc -abs -uthr 5 -mul ${PathResults}/VentrLin.nii.gz  -add ${PathResults}/VentrLinSP.nii.gz ${PathResults}/VentrLinSP.nii.gz
# ${PathSeg}/seg_maths ${PathResults}/VentrLinSP.nii.gz -add ${PathResults}/Artefacts_${ID}.nii.gz -mul -1 -add ${PathResults}/MergedLesion_${ID}_SPCP.nii.gz -thr 0 ${PathResults}/MergedLesion_${ID}_corr.nii.gz 

${PathSeg}/seg_maths ${PathResults}/VentrLinSP.nii.gz -add ${PathResults}/Artefacts_${ID}.nii.gz -mul -1 -add ${PathResults}/SecondarySPCP_${ID}.nii.gz -thr 0 ${PathResults}/SecondarySPCP_${ID}_corr.nii.gz 

${PathSeg}/seg_maths ${PathResults}/VentrLinSP.nii.gz -add ${PathResults}/Artefacts_${ID}.nii.gz -mul -1 -add ${PathResults}/PrimarySPCP_${ID}.nii.gz -thr 0 ${PathResults}/PrimarySPCP_${ID}_corr.nii.gz 

${PathSeg}/seg_maths ${PathResults}/PrimarySPCP_${ID}_corr.nii.gz -merge 1 4 ${PathResults}/SecondarySPCP_${ID}_corr.nii.gz -tmax ${PathResults}/MergedLesionSPCP_${ID}_corr.nii.gz

echo "Optional second level of correction"
if ((flag_furtherCorr==1))
then
    ${PathSeg}/Seg_Analysis -inArtefact ${PathResults}/${ID}_Artefacts.nii.gz -inTxt2 ${PathResults}/${ModalitiesTot}_BiASM_${ID}_${Opt}.txt ${PathResults}/${ModalitiesTot}_BiASM_${ID}_${Opt}.nii.gz -correct -connect -inLesCorr ${PathResults}/MergedLesionSPCP_${ID}_corr.nii.gz  -mask ${PathResults}/GIF_${ID}_B1.nii.gz -inPriorsICSF ${PathResults}/ICBM_ICSF_${ID}.nii.gz -inPriorsDGM ${PathResults}/ICBM_DGM_${ID}.nii.gz -inPriorsCGM ${PathResults}/ICBM_CGM_${ID}.nii.gz -inPriorsECSF ${PathResults}/ICBM_ECSF_${ID}.nii.gz -outConnect 1 -outWM 1 -ParcellationIn ${PathResults}/GIF_Parcellation_${ID}.nii.gz  -inChangePath ${PathResults}/ -inModa 2 1 3 -WeightedSeg 3 3 1 -LesWMI ${OptWMI} -Neigh 6

fi

if ((flag_corrDGM==1))
then
    Array=(24   31  32  56  57  58  59  60  61  76  77  37  38)

# Get segmentation of DGM from GIF

if [ ! -f ${PathResults}/${ID}_DGM.nii.gz ]
then
if ((${#Array[@]}>0))
then
stringAddition="${PathResults}/${ID}_DGMConstruction_0.nii.gz "
for ((k=0;k<${#Array[@]};k++))
do
Value=${Array[k]}
${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz -equal ${Value} -bin ${PathResults}/${ID}_DGMConstruction_${k}.nii.gz 
stringAddition="${stringAddition} -add ${PathResults}/${ID}_DGMConstruction_${k}.nii.gz "
done
${PathSeg}/seg_maths ${stringAddition} -bin -odt char ${PathResults}/${ID}_DGM.nii.gz 
rm ${PathResults}/*DGMConstruction* 
fi
fi
if [ ! -f ${PathResults}/WMminIn* ]
then
${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz -thr 81.5 -uthr 83.5 ${PathResults}/Insula_${ID}.nii.gz 
${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz -thr 89.5 -uthr 91.5 -add ${PathResults}/Insula_${ID}.nii.gz -bin  ${PathResults}/Insula_${ID}.nii.gz 
${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz -thr 95.5 -uthr 97.5 -add ${PathResults}/Insula_${ID}.nii.gz -bin  ${PathResults}/Insula_${ID}.nii.gz 
${PathSeg}/seg_maths ${PathResults}/GIF_Parcellation_${ID}.nii.gz -thr 79 -uthr 98 -bin -sub ${PathResults}/Insula_${ID}.nii.gz ${PathResults}/WMminIn_${ID}.nii.gz
fi
${PathSeg}/seg_maths ${PathResults}/*Parcellation* -thr 23 -uthr 45 -bin ${PathResults}/InfraDGM_${ID}.nii.gz
${PathSeg}/seg_maths ${PathResults}/${ID}_DGM* -add ${PathResults}/Infra* -add ${PathResults}/Insula* -bin -dil 1 -sub ${PathResults}/WMminIn_${ID}.nii.gz -thr 0 ${PathResults}/Correction_${ID}.nii.gz
read -a Mahal <<< $(ls ${PathResults}/LesionMahal*)
read -a Lesion <<< $(ls ${PathResults/Corr*Mer*})
${PathSeg}/seg_maths ${Mahal[0]} -tp 2 -uthr 4 -bin -mul ${Lesion[0]} -mul ${PathResults}/Correction_${ID}.nii.gz -mul -1 -add ${Lesion[0]} -thr 0 ${PathResults}/Lesion_${ID}_corr.nii.gz

${PathSeg}/Seg_Analysis -LesWMI ${OptWMI}  -inLesCorr ${PathResults}/Lesion_${ID}_corr.nii.gz -inTxt2 ${PathResults}/${ModalitiesTot}_BiASM_${ID}_${Opt}.txt ${PathResults}/${ModalitiesTot}_BiASM_${ID}_${Opt}.nii.gz -mask ${PathResults}/GIF_${ID}_B1.nii.gz -Package 1 -SegType 1 -WeightedSeg 3 3 1 -connect -correct -inModa 2 1 3 -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${PathResults}/ICBM_ICSF_${ID}.nii.gz -inPriorsDGM ${PathResults}/ICBM_DGM_${ID}.nii.gz -inPriorsCGM ${PathResults}/ICBM_CGM_${ID}.nii.gz -inPriorsECSF ${PathResults}/ICBM_ECSF_${ID}.nii.gz -TO 1 -juxtaCorr 1 -SP ${OptSP} -LevelCorrection ${OptCL} -inArtefact ${PathResults}/${ID}_Artefacts.nii.gz ${ChangePathString} -ParcellationIn ${PathResults}/GIF_Parcellation_${ID}.nii.gz -typeVentrSeg 1 -outWM 1 -outConnect 1 -Neigh 6

fi

rm ${PathResults}/DataR* ${PathResults}/LesionT* ${PathResults}/Summ*



#${PathSeg}/Seg_Analysis -inToAnalyse ${PathResults}/CorrectMergedLesion_${ID}_corr.nii.gz -minSize 3 -Neigh 6

