
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
OptWMI=2 # If >0, indicates that a higher level of sensitivity to lesion should be considered from voxels coming from WMI
TVS=1
arrayMod=("T1" "FLAIR" )
arrayModNumber=(1 3)
# ArtefactArray=(101 102 139 140 187 188 167 168) 
# ArtefactArray=(5 12 32 33 39 40 50 51 62 63 64 65 72 73 74 48 49 101 102 117 118 139 140 171 172 187 188 167 168 202 203 181 182 117 118 123 124 133 134 155 156 181 182 185 186 201 202 203 204 207 208)
# ListCorr=(5 12 32 33 39 40 50 51 62 63 64 65 72 73 74 48 49 101 102 117 118 139 140 171 172 187 188 167 168 202 203 181 182 117 118 123 124 133 134 155 156 181 182 185 186 201 202 203 204 207 208)
# #ArtefactArray=(101 102 139 140 187 188 167 168 39 40 117 118 123 124 125 126 133 134 137 138 141 142 147 148 155 156 181 182 185 186 187 188 201 202 203 204 207 208)


ListCorr=(101 102 103 104 105 106 187 188 47 48 49 32 33 173 174)
ArtefactArray=(5 12 24 31 39 40 50 51 62 63 64 65 72 73 74 48
                    49 101 102 105 106 139 140 187 188 167 168 39 40 117 118 123 124 125 126 133 134 137 138 141 142 147 148 155 156 181 182 185 186 187 188 201 202 203 204 207 208)

RuleFileName=/cluster/project0/SegBiASM/DataToTryBaMoS/GenericRule_CSF.txt
NameGMatrix=/cluster/project0/SegBiASM/DataToTryBaMoS/GMatrix4_Low3.txt

# RuleFileName=/mounts/auto/xnat/pipelines/BiASM/GenericRule_CSF.txt

# NameGMatrix=/mounts/auto/xnat/pipelines/BiASM/GMatrix4_Low3.txt
Mem=7.9
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
OptTxt=Test
Opt=TA
# Flags for levels of correction
flag_NoCerrebellum=0
flag_furtherCorr=1
flag_corrDGM=1
 
PathReg=/share/apps/cmic/niftyreg-10.4.15/bin
PathReg=/share/apps/cmic/niftyreg_v1.5.43/bin
PathSeg=/home/csudre/NiftySeg_0.9.4/build_comic2/seg-apps
PathFSL=/share/apps/fsl-5.0.8/bin
PathICBM=/cluster/project0/SegBiASM/ICBM_Priors
# PathReg=/Users/csudre/Development/niftyreg_install//bin/
# PathSeg=/Users/csudre/Development/NiftySeg_build_Debug//seg-apps/
# PathFSL=$FSLDIR/bin
# PathICBM=/Users/csudre/Documents/Priors/ICBM_Priors


# RuleFileName=/cluster/project0/SegBiASM/DataToTryBaMoS/GenericRule_CSF.txt
# NameGMatrix=/cluster/project0/SegBiASM/DataToTryBaMoS/GMatrix4_Low3.txt
# NameGMatrix=/Users/csudre/Dropbox/Scripts/GMatrix4_Low3.txt
# RuleFileName=/Users/csudre/Dropbox/Scripts/GenericRule_CSF.txt

if [ $# -lt 5 ] || [ $# -gt 10 ]
then
	echo ""
	echo "*******************************************************************************"
	echo "Usage: $0 ID ImageFLAIR ImageT1 GIF_results_path PathResults JUMP_START Opt Mem Space"
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
Opt=$7
Mem=$8
Space=$9
Verbose=$10
PN=${ID}
ModalitiesTot=T1FLAIR
echo ${PathScratch}
PathGIF=${GIF_results_path}

if [ ! -d $PathResults ]
then
    mkdir $PathResults
fi

NameScript=${PathResults}/ScriptBaMoS_${PN}_${Opt}.sh
echo \#\!/bin/sh > ${NameScript}
echo "#$ -S /bin/bash" > ${NameScript}
PathScratch=/scratch0/csudre/${ID}_BaMoSCross_${Opt}_${JOB_ID}
echo "mkdir -p ${PathScratch} ">> ${NameScript}
ChangePathString="-inChangePath ${PathScratch}/"

echo "cp ${PathResults}/* ${PathScratch}/." >> ${NameScript}


StringPriorsIO=" "

echo "echo "Assuming T1 and GIF outputs are already in FLAIR Space" " >> ${NameScript}
    echo "cp ${ImageT1} ${PathScratch}/T1_${ID}.nii.gz" >> ${NameScript}
    #echo "${PathFSL}/fslreorient2std ${ImageT1} ${PathScratch}/T1_${ID}.nii.gz" >> ${NameScript}
   
    echo "cp ${GIF_results_path}/*TIV* ${PathScratch}/Mask_${ID}.nii.gz" >> ${NameScript}
    #echo "${PathFSL}/fslreorient2std ${GIF_results_path}/*TIV* ${PathScratch}/Mask_${ID}.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/Mask_${ID}.nii.gz -bin -odt char ${PathScratch}/GIF_${ID}_B1.nii.gz" >> ${NameScript}

    
    echo "cp ${GIF_results_path}/*prior* ${PathScratch}/GIF_prior_${ID}.nii.gz" >> ${NameScript}
    #echo "${PathFSL}/fslreorient2std ${GIF_results_path}/*prior* ${PathScratch}/GIF_prior_${ID}.nii.gz" >> ${NameScript}

echo "cp ${ImageFLAIR} ${PathScratch}/FLAIR_${ID}.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/FLAIR_${ID}.nii.gz -odt float ${PathScratch}/FLAIR_${ID}.nii.gz" >> ${NameScript}




for p in Segmentation prior TIV Parcellation
do
echo "cp ${PathGIF}/*${p}* ${PathScratch}/GIF_${ID}_${p}.nii.gz" >> ${NameScript}
echo "cp ${PathGIF}/*${p}* ${PathScratch}/GIF_${p}_${ID}.nii.gz" >> ${NameScript}
done



echo "${PathSeg}/seg_maths ${PathScratch}/GIF_${ID}_TIV.nii.gz -odt char ${PathScratch}/GIF_${ID}_B1.nii.gz" >>${NameScript}
echo "${PathSeg}/seg_maths ${PathScratch}/GIF_${ID}_TIV.nii.gz -odt char ${PathScratch}/GIF_B1_${ID}.nii.gz" >> ${NameScript}

#  If GIF not performed in FLAIR space, need to resample parcellation to FLAIR space
#    echo "echo "Reorientation of GIF parcellation"" >> ${NameScript}
 #    echo "cp ${GIF_results_path}/*Parcellation* ${PathScratch}/GIF_Parcellation_${ID}.nii.gz " >> ${NameScript}
    #echo "${PathFSL}/fslreorient2std ${GIF_results_path}/*Parcellation* ${PathScratch}/GIF_Parcellation_${ID}.nii.gz " >> ${NameScript}
echo "echo "Reorientation of GIF segmentation"" >> ${NameScript}
 #echo "cp ${GIF_results_path}/*Segmentation* ${PathScratch}/GIF_Segmentation_${ID}.nii.gz " >> ${NameScript}
  #  echo "${PathFSL}/fslreorient2std ${GIF_results_path}/*Segmentation* ${PathScratch}/GIF_Segmentation_${ID}.nii.gz " >> ${NameScript}

    echo "cp ${PathScratch}/* ${PathResults}/." >> ${NameScript}

    PriorsArray=("Out" "CSF" "CGM" "WMI" "DGM" "Brainstem")
    array_Priors=()

    for ((p=0;p<6;p++))
    do
       echo " ${PathSeg}/seg_maths ${PathScratch}/GIF_prior* -tp ${p} ${PathScratch}/GIF_${PriorsArray[p]}_${ID}.nii.gz" >> ${NameScript}
    done

    
    

    echo "echo "Creation of artefact map based on T1 parcellation and registration to FLAIR"" >> ${NameScript}

    if ((${#ArtefactArray[@]}>0))
    then
        stringAddition="${PathScratch}/${ID}_ArtConstruction_0.nii.gz "
        for ((i=0;i<${#ArtefactArray[@]};i++))
        do
            Value=${ArtefactArray[i]}
            ValueMin=`echo "$Value - 0.5" | bc -l`
            ValueMax=`echo "$Value + 0.5" | bc -l`
            echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz -thr $ValueMin -uthr $ValueMax -bin ${PathScratch}/${ID}_ArtConstruction_${i}.nii.gz" >> ${NameScript}
            stringAddition="${stringAddition} -add ${PathScratch}/${ID}_ArtConstruction_${i}.nii.gz "
        done
        echo "${PathSeg}/seg_maths ${stringAddition} -bin ${PathScratch}/${ID}_Artefacts.nii.gz" >> ${NameScript}
    fi

    echo "rm ${PathScratch}/*_ArtConstruction_*" >> ${NameScript}

if [ "$JUMP_START" -eq 0 ]
then
    echo "echo "Registration of ICBM template and creation of ICBM atlases"" >> ${NameScript}
    echo "${PathReg}/reg_aladin -ref ${PathScratch}/T1_${ID}.nii.gz -flo ${PathICBM}/ICBM_Template.nii.gz -aff ${PathScratch}/${ID}_AffTransf.txt" >> ${NameScript}
    echo "${PathReg}/reg_f3d -ref ${PathScratch}/T1_${ID}.nii.gz -flo ${PathICBM}/ICBM_Template.nii.gz -aff ${PathScratch}/${ID}_AffTransf.txt -cpp ${PathScratch}/${ID}_cpp.nii.gz" >> ${NameScript}

    echo "cp ${PathScratch}/*cpp* ${PathResults}/." >> ${NameScript} 

    for p in CGM DGM ECSF ICSF Out WM
    do
        echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${ID}.nii.gz -flo ${PathICBM}/ICBM_${p}.nii.gz -cpp ${PathScratch}/${ID}_cpp.nii.gz -res ${PathScratch}/ICBM_${p}_${ID}.nii.gz" >> ${NameScript}
    done

    echo "${PathSeg}/seg_maths ${PathScratch}/ICBM_DGM_${ID}.nii.gz -add ${PathScratch}/ICBM_CGM_${ID}.nii.gz ${PathScratch}/ICBM_GM_${ID}.nii.gz" >> ${NameScript}

    echo "${PathSeg}/seg_maths ${PathScratch}/ICBM_ICSF_${ID}.nii.gz -add ${PathScratch}/ICBM_ECSF_${ID}.nii.gz ${PathScratch}/ICBM_CSFs_${ID}.nii.gz" >> ${NameScript}

    echo "${PathSeg}/seg_maths ${PathScratch}/ICBM_DGM_${ID}.nii.gz -bin -mul ${PathScratch}/GIF_DGM_${ID}.nii.gz ${PathScratch}/GIF_DGM_${ID}_bis.nii.gz" >> ${NameScript}
    
echo "${PathSeg}/seg_maths ${PathScratch}/GIF_CGM_${ID}.nii.gz -add ${PathScratch}/GIF_DGM_${ID}_bis.nii.gz ${PathScratch}/GIF_GM_${ID}.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_DGM_${ID}.nii.gz -sub ${PathScratch}/GIF_DGM_${ID}_bis.nii.gz -add ${PathScratch}/GIF_WMI_${ID}.nii.gz ${PathScratch}/GIF_WMI_${ID}_bis.nii.gz" >> ${NameScript}

    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_WMI_${ID}_bis.nii.gz -add ${PathScratch}/GIF_Brainstem_${ID}.nii.gz ${PathScratch}/GIF_WM_${ID}.nii.gz" >> ${NameScript}

    

    arrayImage=()
    arrayModNumber=()
    echo "echo "Modalities to put are ${arrayMod[0]}"" >> ${NameScript}
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
        arrayImage=( "${arrayImage[*]}" "${PathScratch}/${arrayMod[m]}_${ID}.nii.gz" )
    done

array_Priors=("${PathScratch}/GIF_GM_${ID}.nii.gz" "${PathScratch}/GIF_WM_${ID}.nii.gz" "${PathScratch}/GIF_CSF_${ID}.nii.gz" "${PathScratch}/GIF_Out_${ID}.nii.gz")
    echo "echo "Segmentation SegBiASM"" >> ${NameScript}
    #${PathSeg}/Seg_BiASM -in 2 ${arrayImage[0]} ${arrayImage[1]} -priors 4 ${array_Priors[*]} -mask ${PathScratch}/GIF_${ID}_B1.nii.gz -out 2 ${PathScratch}/${arrayMod[0]}${arrayMod[1]}_BiASM_${ID}_${Opt}.nii.gz ${PathScratch}/${arrayMod[0]}${arrayMod[1]}_BiASMG_${ID}.nii.gz -txt_out ${PathScratch}/${arrayMod[0]}${arrayMod[1]}_BiASM_${ID}_${Opt}.txt -bc_order 3 -CovPriors 8 -BFP 1 -maxRunEM ${MRE} -AtlasSmoothing 1 1 -AtlasWeight 1 ${AW}  -SMOrder 0 -KernelSize 3 -PriorsKept 5 -VLkappa 1.5 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -deleteUW 1 -outliersM 3 -outliersW ${OW} -init_splitUnif 0 -splitAccept 0 -unifTot 1 -MRF 1 -GMRF ${NameGMatrix} -juxtaCorr ${JC} -progMod 0 -priorDGM ${PathScratch}/ICBM_DGM_${ID}.nii.gz -TypicalityAtlas ${OptTA} ${StringPriorsIO}




echo "${PathSeg}/Seg_BiASM -VLkappa 3 -in 2 ${arrayImage[*]} -priors 4 ${array_Priors[*]} -mask ${PathScratch}/GIF_${ID}_B1.nii.gz -out 2 ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${Opt}.nii.gz ${PathScratch}/${ModalitiesTot}_BiASMG_${ID}.nii.gz -txt_out ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${Opt}.txt -bc_order 3 -CovPriors 8 -BFP 1 -maxRunEM ${MRE} -AtlasSmoothing 1 1 -AtlasWeight 1 ${AW}  -SMOrder 0 -KernelSize 3 -PriorsKept 5 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -deleteUW 1 -outliersM 3 -outliersW ${OW} -init_splitUnif 0 -splitAccept 0 -unifTot 1 -MRF 1 -GMRF ${NameGMatrix} -juxtaCorr ${JC} -progMod 0 -priorDGM ${PathScratch}/ICBM_DGM_${ID}.nii.gz -TypicalityAtlas ${OptTA}" >> ${NameScript}

echo "rm ${PathScratch}/BG* ${PathScratch}/MRF*" >> ${NameScript}
echo "cp ${PathScratch}/T1FLAIR* ${PathResults}/. " >> ${NameScript}
echo "cp ${PathScratch}/Data*T1FLAIR* ${PathResults}/. " >> ${NameScript}

    echo "echo "Lesion segmentation"" >> ${NameScript}
    echo "${PathSeg}/Seg_Analysis  -inTxt2 ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${Opt}.txt ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${Opt}.nii.gz -mask ${PathScratch}/GIF_${ID}_B1.nii.gz -Package 1 -SegType 1 -WeightedSeg 3 3 1 -connect -correct -inModa 2 1 3 -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${PathScratch}/ICBM_ICSF_${ID}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${ID}.nii.gz -inPriorsCGM ${PathScratch}/ICBM_CGM_${ID}.nii.gz -inPriorsECSF ${PathScratch}/ICBM_ECSF_${ID}.nii.gz -TO 1 -ParcellationIn ${PathScratch}/GIF_Parcellation_${ID}.nii.gz -typeVentrSeg ${TVS} -Simple -Secondary 60 -juxtaCorr ${JC}  -SP ${OptSP} -LevelCorrection ${OptCL} ${ChangePathString} -LesWMI ${OptWMI} -Neigh 18" >> ${NameScript}

fi

if [ "$JUMP_START" -eq 1 ]
then

    for p in CGM DGM ECSF ICSF Out WM
    do
        echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${ID}.nii.gz -flo ${PathICBM}/ICBM_${p}.nii.gz -cpp ${PathScratch}/${ID}_cpp.nii.gz -res ${PathScratch}/ICBM_${p}_${ID}.nii.gz" >> ${NameScript}
    done

        PriorsArray=("Out" "CSF" "CGM" "WMI" "DGM" "Brainstem")
    array_Priors=()

    for ((p=0;p<6;p++))
    do
       echo " ${PathSeg}/seg_maths ${PathScratch}/GIF_prior* -tp ${p} ${PathScratch}/GIF_${PriorsArray[p]}_${ID}.nii.gz" >> ${NameScript}
    done

    echo "${PathSeg}/seg_maths ${PathScratch}/ICBM_DGM_${ID}.nii.gz -add ${PathScratch}/ICBM_CGM_${ID}.nii.gz ${PathScratch}/ICBM_GM_${ID}.nii.gz" >> ${NameScript}

    echo "${PathSeg}/seg_maths ${PathScratch}/ICBM_ICSF_${ID}.nii.gz -add ${PathScratch}/ICBM_ECSF_${ID}.nii.gz ${PathScratch}/ICBM_CSFs_${ID}.nii.gz" >> ${NameScript}

    echo "${PathSeg}/seg_maths ${PathScratch}/ICBM_DGM_${ID}.nii.gz -bin -mul ${PathScratch}/GIF_DGM_${ID}.nii.gz ${PathScratch}/GIF_DGM_${ID}_bis.nii.gz" >> ${NameScript}
    
echo "${PathSeg}/seg_maths ${PathScratch}/GIF_CGM_${ID}.nii.gz -add ${PathScratch}/GIF_DGM_${ID}_bis.nii.gz ${PathScratch}/GIF_GM_${ID}.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_DGM_${ID}.nii.gz -sub ${PathScratch}/GIF_DGM_${ID}_bis.nii.gz -add ${PathScratch}/GIF_WMI_${ID}.nii.gz ${PathScratch}/GIF_WMI_${ID}_bis.nii.gz" >> ${NameScript}

    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_WMI_${ID}_bis.nii.gz -add ${PathScratch}/GIF_Brainstem_${ID}.nii.gz ${PathScratch}/GIF_WM_${ID}.nii.gz" >> ${NameScript}
fi

 echo "${PathSeg}/Seg_Analysis  -inTxt2 ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${Opt}.txt ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${Opt}.nii.gz -mask ${PathScratch}/GIF_${ID}_B1.nii.gz -Package 1 -SegType 1 -WeightedSeg 3 3 1 -connect -correct -inModa 2 1 3 -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${PathScratch}/ICBM_ICSF_${ID}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${ID}.nii.gz -inPriorsCGM ${PathScratch}/ICBM_CGM_${ID}.nii.gz -inPriorsECSF ${PathScratch}/ICBM_ECSF_${ID}.nii.gz -TO 1 -ParcellationIn ${PathScratch}/GIF_Parcellation_${ID}.nii.gz -typeVentrSeg ${TVS} -Simple -Secondary 60 -juxtaCorr ${JC}  -SP ${OptSP} -LevelCorrection ${OptCL} ${ChangePathString} -LesWMI ${OptWMI} -Neigh 18" >> ${NameScript}
#echo "echo "cp ${PathScratch}/LesionCorrected*WS3WT3WC1* ${PathScratch}/PrimaryLesions_${ID}.nii.gz"" >> ${NameScript}
 #echo "echo  " cp ${PathScratch}/SecondaryCorrected*WS3WT3WC1* ${PathScratch}/SecondaryLesions_${ID}.nii.gz "" >> ${NameScript}
    echo "cp ${PathScratch}/LesionCorrected*WS3WT3WC1* ${PathScratch}/PrimaryLesions_${ID}.nii.gz" >> ${NameScript}
    echo "cp ${PathScratch}/SecondaryCorrected*WS3WT3WC1* ${PathScratch}/SecondaryLesions_${ID}.nii.gz" >> ${NameScript}

    echo "cp ${PathScratch}/Primary* ${PathResults}/. " >> ${NameScript}
    echo "cp ${PathScratch}/Secondary* ${PathResults}/. " >> ${NameScript}

    echo "rm ${PathScratch}/DataR* ${PathScratch}/LesionT* ${PathScratch}/Summ*" >> ${NameScript}

    echo "echo "Correction for septum pellucidum"" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz -thr 65.5 -uthr 67.5 -bin ${PathScratch}/VentricleLining.nii.gz " >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -equal 52 -bin -euc -uthr 0 -abs ${PathScratch}/Temp1.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -equal 53 -bin -euc -uthr 0 -abs -add ${PathScratch}/Temp1.nii.gz -uthr 5 -bin ${PathScratch}/PotentialInterVentr.nii.gz" >> ${NameScript}

    echo "${PathSeg}/seg_maths ${PathScratch}/VentricleLining.nii.gz -dil 1 ${PathScratch}/ExpandedVentrLin.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -thr 49.5 -uthr 53.5 -bin -sub ${PathScratch}/ExpandedVentrLin.nii.gz -thr 0 -bin ${PathScratch}/PotentialCP.nii.gz" >> ${NameScript}

    #${PathSeg}/seg_maths ${GIF_results_path}/*Parcellation*  -thr 122.5 -uthr 124.5 -bin ${PathScratch}/Temp4.nii.gz

    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Segmentation_${ID}.nii.gz -tp 4 -thr 0.5 -bin -sub ${PathScratch}/VentricleLining* -thr 0 ${PathScratch}/SegDGM.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -equal 87 -mul -1 -add ${PathScratch}/PotentialInterVentr.nii.gz -sub ${PathScratch}/SegDGM.nii.gz -thr 0 ${PathScratch}/PotentialSP3.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths  ${PathScratch}/PotentialSP3.nii.gz -add ${PathScratch}/PotentialCP.nii.gz ${PathScratch}/PotentialSPCP.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/ICBM_ICSF* -thr 0.3 -bin -mul ${PathScratch}/PotentialSPCP.nii.gz -add ${PathScratch}/PotentialSP3.nii.gz -thr 0.2  -bin  ${PathScratch}/PotentialSPCP2.nii.gz" >> ${NameScript}

    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz -thr 83.5 -uthr 84.5 -bin -dil 5 ${PathScratch}/WMDil.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz -thr 91.5 -uthr 92.5 -bin -dil 5 ${PathScratch}/WMDil2.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/PotentialSPCP2.nii.gz -sub ${PathScratch}/WMDil.nii.gz -sub ${PathScratch}/WMDil2.nii.gz -thr 0 ${PathScratch}/PotentialSPCP3.nii.gz " >> ${NameScript}

    echo "${PathSeg}/seg_maths ${PathScratch}/PrimaryLesions_* -sub ${PathScratch}/PotentialSPCP3.nii.gz -thr 0 ${PathScratch}/PrimarySPCP_${ID}_${OptTxt}.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/SecondaryLesions_* -sub ${PathScratch}/PotentialSPCP3.nii.gz -thr 0 ${PathScratch}/SecondarySPCP_${ID}_${OptTxt}.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/PrimarySPCP*${OptTxt}* -merge 1 4 ${PathScratch}/SecondarySPCP_*${OptTxt}* -tmax ${PathScratch}/MergedLesion_${ID}_SPCP_${OptTxt}.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathResults}/MergedLesion_${ID}_SPCP_${OptTxt}.nii.gz -sub ${PathResults}/${ID}_Artefacts.nii.gz -thr 0 ${PathResults}/MergedLesion_${ID}_SPCP_${OptTxt}.nii.gz" >> ${NameScript}


    echo "echo "Refined correction for third ventricle"" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -sub ${PathScratch}/GIF_Parcellation_${ID}.nii.gz ${PathScratch}/Artefacts_${ID}.nii.gz " >> ${NameScript}

for ((i=0;i<${#ListCorr[@]};i++))
do
echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -equal ${ListCorr[i]} -bin -add ${PathScratch}/Artefacts_${ID}.nii.gz " >> ${NameScript}
done
echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -thr 65.5 -uthr 67.5 ${PathScratch}/VentrLin.nii.gz " >> ${NameScript}
echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -equal 87 -add ${PathScratch}/VentrLin.nii.gz ${PathScratch}/VentrLin.nii.gz " >> ${NameScript}
echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -equal 47 -euc  -abs -uthr 5 -bin -mul ${PathScratch}/VentrLin.nii.gz ${PathScratch}/VentrLinSP.nii.gz " >> ${NameScript}

echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -equal 52 -bin -euc  -abs  ${PathScratch}//Ventr1.nii.gz " >> ${NameScript}

echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -equal 53 -bin -euc -abs -uthr 5 -sub  ${PathScratch}/Ventr1.nii.gz -thr -5 -uthr 5 -bin -mul  ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -equal 52 -bin -euc -abs -uthr 5 -mul ${PathScratch}/VentrLin.nii.gz  -add ${PathScratch}/VentrLinSP.nii.gz ${PathScratch}/VentrLinSP.nii.gz" >> ${NameScript}
echo "${PathSeg}/seg_maths ${PathScratch}/VentrLinSP.nii.gz -add ${PathScratch}/Artefacts_${ID}.nii.gz -mul -1 -add ${PathScratch}/MergedLesion_${ID}_SPCP_${OptTxt}.nii.gz -thr 0 ${PathScratch}/MergedLesion_${ID}_${OptTxt}_corr.nii.gz " >> ${NameScript}

echo "${PathSeg}/Seg_Analysis -LesWMI ${OptWMI} ${OptInfarcts} -inLesCorr ${PathScratch}/MergedLesion_${ID}_${OptTxt}_corr.nii.gz  -inTxt2 ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${Opt}.txt ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${Opt}.nii.gz -mask ${PathScratch}/GIF_${ID}_B1.nii.gz -Package 1 -SegType 1 -WeightedSeg 3 3 1 -connect -correct -inModa ${#arrayModNumber[@]} ${arrayModNumber[*]} -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${PathScratch}/ICBM_ICSF_${ID}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${ID}.nii.gz -inPriorsCGM ${PathScratch}/ICBM_CGM_${ID}.nii.gz -inPriorsECSF ${PathScratch}/ICBM_ECSF_${ID}.nii.gz -TO 1 -juxtaCorr 1 -SP ${OptSP} -LevelCorrection ${OptCL} -inArtefact ${PathScratch}/${ID}_Artefacts.nii.gz ${ChangePathString} -ParcellationIn ${PathScratch}/GIF_Parcellation_${ID}.nii.gz -typeVentrSeg 1 -outWM 1 -outConnect 1 -Neigh 6" >> ${NameScript}

echo "rm ${PathScratch}/LesionWeigh* ${PathScratch}/Binary* ${PathScratch}/WMDil* ${PathScratch}/WMCard* ${PathScratch}/LesionInit* ${PathScratch}/DataR* ${PathScratch}/DataT* ${PathScratch}/Summ* ${PathScratch}/LesSegHard* ${PathScratch}/Check* ${PathScratch}/BinaryNIV*" >> ${NameScript}



# ${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -equal 52 -bin -euc  -abs  ${PathScratch}//Ventr1.nii.gz 

# ${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -equal 53 -bin -euc -abs -uthr 5 -sub  ${PathScratch}/Ventr1.nii.gz -thr -5 -uthr 5 -bin -mul  ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -equal 52 -bin -euc -abs -uthr 5 -mul ${PathScratch}/VentrLin.nii.gz  -add ${PathScratch}/VentrLinSP.nii.gz ${PathScratch}/VentrLinSP.nii.gz
# ${PathSeg}/seg_maths ${PathScratch}/VentrLinSP.nii.gz -add ${PathScratch}/Artefacts_${ID}.nii.gz -mul -1 -add ${PathScratch}/MergedLesion_${ID}_SPCP.nii.gz -thr 0 ${PathScratch}/MergedLesion_${ID}_corr.nii.gz 

# echo "${PathSeg}/seg_maths ${PathScratch}/VentrLinSP.nii.gz -add ${PathScratch}/Artefacts_${ID}.nii.gz -mul -1 -add ${PathScratch}/SecondarySPCP_${ID}.nii.gz -thr 0 ${PathScratch}/SecondarySPCP_${ID}_corr.nii.gz " >> ${NameScript}

# echo "${PathSeg}/seg_maths ${PathScratch}/VentrLinSP.nii.gz -add ${PathScratch}/Artefacts_${ID}.nii.gz -mul -1 -add ${PathScratch}/PrimarySPCP_${ID}.nii.gz -thr 0 ${PathScratch}/PrimarySPCP_${ID}_corr.nii.gz " >> ${NameScript}

# echo "${PathSeg}/seg_maths ${PathScratch}/PrimarySPCP_${ID}_corr.nii.gz -merge 1 4 ${PathScratch}/SecondarySPCP_${ID}_corr.nii.gz -tmax ${PathScratch}/MergedLesionSPCP_${ID}_corr.nii.gz" >> ${NameScript}

# echo "echo "Optional second level of correction"" >> ${NameScript}
# if ((flag_furtherCorr==1))
# then
#     echo "${PathSeg}/Seg_Analysis -inArtefact ${PathScratch}/Artefacts_${ID}.nii.gz -inTxt2 ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${Opt}.txt ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${Opt}.nii.gz -correct -connect -inLesCorr ${PathScratch}/MergedLesionSPCP_${ID}_corr.nii.gz  -mask ${PathScratch}/GIF_${ID}_B1.nii.gz -inPriorsICSF ${PathScratch}/ICBM_ICSF_${ID}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${ID}.nii.gz -inPriorsCGM ${PathScratch}/ICBM_CGM_${ID}.nii.gz -inPriorsECSF ${PathScratch}/ICBM_ECSF_${ID}.nii.gz -outConnect 1 -outWM 1 -ParcellationIn ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -inChangePath ${PathScratch}/ -inModa 2 1 3 -WeightedSeg 3 3 1 -LesWMI ${OptWMI} -Neigh 6" >> ${NameScript}

# fi

if ((flag_corrDGM==1))
then
    Array=(24   31  32  56  57  58  59  60  61  76  77  37  38)

# Get segmentation of DGM from GIF

if [ ! -f ${PathScratch}/${ID}_DGM.nii.gz ]
then
if ((${#Array[@]}>0))
then
stringAddition="${PathScratch}/${ID}_DGMConstruction_0.nii.gz "
for ((k=0;k<${#Array[@]};k++))
do
Value=${Array[k]}
echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz -equal ${Value} -bin ${PathScratch}/${ID}_DGMConstruction_${k}.nii.gz " >> ${NameScript}
stringAddition="${stringAddition} -add ${PathScratch}/${ID}_DGMConstruction_${k}.nii.gz "
done
echo "${PathSeg}/seg_maths ${stringAddition} -bin -odt char ${PathScratch}/${ID}_DGM.nii.gz " >> ${NameScript}
echo "rm ${PathScratch}/*DGMConstruction* " >> ${NameScript}
fi
fi
if [ ! -f ${PathScratch}/WMminIn* ]
then
echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz -thr 81.5 -uthr 83.5 ${PathScratch}/Insula_${ID}.nii.gz " >> ${NameScript}
echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz -thr 89.5 -uthr 91.5 -add ${PathScratch}/Insula_${ID}.nii.gz -bin  ${PathScratch}/Insula_${ID}.nii.gz " >> ${NameScript}
echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz -thr 95.5 -uthr 97.5 -add ${PathScratch}/Insula_${ID}.nii.gz -bin  ${PathScratch}/Insula_${ID}.nii.gz " >> ${NameScript}
echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz -thr 79 -uthr 98 -bin -sub ${PathScratch}/Insula_${ID}.nii.gz ${PathScratch}/WMminIn_${ID}.nii.gz" >> ${NameScript}
fi
echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz -thr 23 -uthr 45 -bin ${PathScratch}/InfraDGM_${ID}.nii.gz" >> ${NameScript}
echo "${PathSeg}/seg_maths ${PathScratch}/${ID}_DGM* -add ${PathScratch}/Infra* -add ${PathScratch}/Insula* -bin -dil 1 -sub ${PathScratch}/WMminIn_${ID}.nii.gz -thr 0 ${PathScratch}/Correction_${ID}.nii.gz" >> ${NameScript}

echo "${PathSeg}/seg_maths ${PathScratch}/LesionMahal* -tp 2 -uthr 4 -bin -mul ${PathScratch}/Corr*Mer* -mul ${PathScratch}/Correction_${ID}.nii.gz -mul -1 -add ${PathScratch}/Corr*Mer* -thr 0 ${PathScratch}/Lesion_${ID}_corr.nii.gz" >> ${NameScript}

echo "${PathSeg}/Seg_Analysis -LesWMI ${OptWMI}  -inLesCorr ${PathScratch}/Lesion_${ID}_corr.nii.gz -inTxt2 ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${Opt}.txt ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${Opt}.nii.gz -mask ${PathScratch}/GIF_${ID}_B1.nii.gz -Package 1 -SegType 1 -WeightedSeg 3 3 1 -connect -correct -inModa 2 1 3 -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${PathScratch}/ICBM_ICSF_${ID}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${ID}.nii.gz -inPriorsCGM ${PathScratch}/ICBM_CGM_${ID}.nii.gz -inPriorsECSF ${PathScratch}/ICBM_ECSF_${ID}.nii.gz -TO 1 -juxtaCorr 1 -SP ${OptSP} -LevelCorrection ${OptCL} -inArtefact ${PathScratch}/Artefacts_${ID}.nii.gz ${ChangePathString} -ParcellationIn ${PathScratch}/GIF_Parcellation_${ID}.nii.gz -typeVentrSeg 1 -outWM 1 -outConnect 1 -Neigh 6" >> ${NameScript}

fi

echo "rm ${PathScratch}/LesionWeigh* ${PathScratch}/Binary* ${PathScratch}/WMDil* ${PathScratch}/WMCard* ${PathScratch}/ICBM* ${PathScratch}/LesionInit* ${PathScratch}/DataR* ${PathScratch}/DataT* ${PathScratch}/Summ* ${PathScratch}/LesSegHard* ${PathScratch}/Check* ${PathScratch}/BinaryNIV*" >> ${NameScript}
echo "cp ${PathScratch}/*Co* ${PathResults}/." >> ${NameScript}
echo "cp ${PathScratch}/*Co*.txt ${PathResults}/." >> ${NameScript}
echo "cp ${PathScratch}/LesionMahal* ${PathResults}/.">>${NameScript}
echo "cp ${PathScratch}/Txt* ${PathResults}/." >> ${NameScript}
echo "cp ${PathScratch}/Out* ${PathResults}/." >> ${NameScript}
echo "cp ${PathScratch}/Autho* ${PathResults}/." >> ${NameScript}
echo "cp ${PathScratch}/*Infar* ${PathResults}/." >> ${NameScript}

echo "cp ${PathScratch}/*Artefacts.nii.gz ${PathResults}/." >> ${NameScript}
function finish {
     rm -rf ${PathScratch}
}
echo "trap finish EXIT ERR " >> ${NameScript}
echo "rm -rf ${PathScratch}" >> ${NameScript}

if [ -z "$JUMP_START" ] || [ "$JUMP_START" -eq 0 ]
then
qsub -l h_rt=20:0:0 -l tmem=${Mem}G -l h_vmem=${Mem}G -N BaMoS${ID} ${NameScript}
else
qsub -l h_rt=2:0:0 -l tmem=${Mem}G -l h_vmem=${Mem}G -N BaMoSLes${ID} ${NameScript}	
fi

function finish {
    rm -rf ${PathScratch}
}

echo "trap finish EXIT ERR" >> ${NameScript}


#${PathSeg}/Seg_Analysis -inToAnalyse ${PathScratch}/CorrectMergedLesion_${ID}_corr.nii.gz -minSize 3 -Neigh 6

