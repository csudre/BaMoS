#!/bin/sh
# set -x
 


if [ $# -gt 12 ]
then
	echo ""
	echo "*******************************************************************************"
	echo "One argument is expected to run this script:"
	echo "- File with contains the altas creation parameters"
	echo "example: $0 param_groupwise_niftyreg.sh "
	echo "*******************************************************************************"
	echo ""
	exit
fi

ArtefactArray=(101 102 139 140 187 188 167 168)
ArtefactArray=(5 12 32 33 39 40 50 51 62 63 64 65 72 73 74 48 49 101 102 117 118 139 140 171 172 187 188 167 168)

ListCorr=(101 102 103 104 105 106 187 188 47 48 49 32 33 173 174)
ArtefactArray=(5 12 24 31 39 40 50 51 62 63 64 65 72 73 74 48
                    49 101 102 105 106 139 140 187 188 167 168 39 40 117 118 123 124 125 126 133 134 137 138 141 142 147 148 155 156 181 182 185 186 187 188 201 202 203 204 207 208)

RuleFileName=./GenericRule_CSF.txt
NameGMatrix=./GMatrix4_Low3.txt
ScriptGW="./GroupwiseReg_BaMoSLong_UB3_orderFit1.sh"

PossibleModalities=(T1 FLAIR T2)
CodeModality=(1 3 2)
OptWMI=2
TVS=1
DO_LES=1
OW=0.01
AW=0
#OW=0.01
#OW=0.10
#OW=0.50
#OW=0.001
#OW=0.99
MRE=250
JC=1
Opt=${OptCross}
OptVL=1
OptSP=1
OptTA=1
OptTxt=Test
Opt=TA
OptCL=2
TVS=1
# Flags for levels of correction
flag_NoCerrebellum=0
flag_furtherCorr=1
flag_corrDGM=1
 

PathReg=PathToNiftyReg/bin
PathSeg=PathToSeg/seg-apps
PathICBM=PathToICBM_Priors



OptArt=1
#############################################################################
# read the input parameters
FinalMod=(T1 FLAIR)

TAString=""
OptTP=""
TAOpt=""

if [ $# -lt 5 ] || [ $# -gt 11 ]
then
    echo ""
    echo "*******************************************************************************"
    echo "Usage: $0 ID ImageFLAIR ImageT1 GIF_results_path TPList PathResults DO_LES JUMP_START OptTP Opt OptGIF"
    echo "*******************************************************************************"
    echo ""
    exit
fi



ArgFLAIR=${2%\\)}
ArgFLAIR=${ArgFLAIR#\\(}
ImageFLAIR=(${ArgFLAIR})
echo $ArgFLAIR


ArgT1=${3%\\)}
echo HalfArgT1 $ArgT1
ArgT1=${ArgT1#\\(}
echo NewArgT1 $ArgT1
ImageT1=(${ArgT1})
echo $ImageT1

ArgGIF=${4%\\)}
ArgGIF=${ArgGIF#\\(}
echo ArgGIF $ArgGIF
GIF_results_path=(${ArgGIF})
echo $GIF_results_path

ArgTP=${5%\\)}
ArgTP=${ArgTP#\\(}
TP=(${ArgTP})
echo $ArgTP


ID=$1
#ImageFLAIR=($2)
#ImageT1=($3)

#GIF_results_path=($4)
#TP=($5)
PathResults=$6
DO_LES=$7
JUMP_START=$8
OptTP=$9
Opt=$10
OptGIF=${11}
ModalitiesTot=T1FLAIR
echo ${PathResults}
echo GIF is ${GIF_results_path[@]}
echo ImageFLAIR is ${ImageFLAIR[@]}
if [ ! -d $PathResults ]
then
    mkdir $PathResults
fi


Space=1 #(FLAIR space)

ChangePathString="-inChangePath ${PathResults}/"
StringPriorsIO=" "

if ((OptTA==1))
then
TAString="-TypicalityAtlas 1"
TAOpt="TA1"
fi

PN=${ID}
echo ${OptTP}
NameScript=${PathResults}/LongWMH_${PN}${OptTP}.sh
echo \#\!/bin/sh > ${NameScript}
NameScriptGroupwise=${PathResults}/GroupwiseScript_${PN}${OptTP}.sh

# Be careful here with Simu : not possible to run for different Tests
PathScratch=/scratch0/csudre/${PN}_Long${OptTP}_${JOB_ID}
echo $PathScratch
ChangePathString="-inChangePath ${PathScratch}/"
echo "mkdir -p ${PathScratch}" >> ${NameScript}
echo "cp -r ${PathResults}/* ${PathScratch}/." >> ${NameScript}
echo ${TP[1]} ${ImageFLAIR[1]} ${ImageT1[1]} ${DO_LES} ${OptGIF}


if ((OptGIF==1))
then
ArrayGIF=(Segmentation Parcellation prior TIV)
else
ArrayGIF=(desc-seg label prior TIV-binary)
fi

for ((tp=0;tp<${#TP[@]};tp++))
do
    echo "cp ${ImageT1[tp]} ${PathScratch}/T1_${ID}_${TP[tp]}.nii.gz">> ${NameScript}
    echo "cp ${ImageFLAIR[tp]} ${PathScratch}/FLAIR_${ID}_${TP[tp]}_init.nii.gz">> ${NameScript}
    ArrayNew=(Segmentation Parcellation prior TIV)
    for p in 0 1 2 3 
    do
        echo "cp ${GIF_results_path[tp]}/*${ArrayGIF[p]}* ${PathScratch}/${ID}_${TP[tp]}_${ArrayNew[p]}.nii.gz">> ${NameScript}
        echo "cp ${GIF_results_path[tp]}/*${ArrayGIF[p]}* ${PathScratch}/GIF_${ID}_${TP[tp]}_${ArrayNew[p]}.nii.gz">> ${NameScript}
        echo "cp ${GIF_results_path[tp]}/*${ArrayGIF[p]}* ${PathScratch}/GIF_${ArrayNew[p]}_${ID}_${TP[tp]}.nii.gz">> ${NameScript}
    done

    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_${PN}_*${TP[tp]}_TIV.nii.gz  -odt char ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz" >> ${NameScript}
done



#Create folders needed for performing task, copy the appropriate data, do the registration and stacking
#mkdir ${PathID}/Study

for ((tp=0;tp<${#TP[@]};tp++))
do



# Create the ArrayModalityCode
    ArrayModaCode=()
    for ((mod=0;mod<${#FinalMod[@]};mod++))
    do
        for ((i=0;i<${#PossibleModalities[@]};i++))
        do
            if [[ ${FinalMod[mod]} == ${PossibleModalities[i]} ]]
            then
                ArrayModaCode=(${ArrayModaCode[*]} ${CodeModality[i]})
            fi
        done
    done
# Copy the initial preprocessed data
#for (( mod=0;mod<${#ArrayCorrespondance[@]}/2;mod++))
#do
#cp ${PathID}/*${PN}_${TP[tp]}*${ArrayCorrespondance[2*mod]}* ${PathID}/${PN}_${TP[tp]}/${ArrayCorrespondance[2*mod+1]}_${PN}_${TP[tp]}_init.nii.gz
#done

#reg_aladin on the FLAIR image and array to make stacking either afterwards
    array_ToStack=()
    NameTogether=""

    for ((mod=0;mod<${#FinalMod[@]};mod++))
    do

        echo "${PathReg}/reg_aladin -ref ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -flo ${PathScratch}/${FinalMod[mod]}_${PN}_${TP[tp]}_init.nii.gz -res ${PathScratch}/${FinalMod[mod]}_${PN}_${TP[tp]}.nii.gz -rigOnly" >> ${NameScript}
 	echo "cp ${PathScratch}/${FinalMod[mod]}_${PN}_${TP[tp]}.nii.gz ${PathResults}/." >> ${NameScript} 
        array_ToStack=(${array_ToStack[*]} "${PathScratch}/${FinalMod[mod]}_${PN}_${TP[tp]}.nii.gz")
        NameTogether=${NameTogether}${FinalMod[mod]}
    done

    SignModa=""
    if ((${#FinalMod[@]}>2))
    then
        SignModa="_${#FinalMod[@]}"
    fi

#Perform the stacking
    let "NM  = ${#array_ToStack[@]}-1 "
    echo "${PathSeg}/seg_maths ${array_ToStack[0]} -merge ${NM} 4 ${array_ToStack[@]:1} ${PathScratch}/ST_${PN}_${TP[tp]}.nii.gz" >> ${NameScript}
		

#All what concerns the T1 only (brain mask + priors registration)
#Beginning with Brainmask

#echo "${PathSeg}/seg_maths ${PathID}/${PN}_${TP[tp]}/T1_${PN}_${TP[tp]}.nii.gz -thr 0 -bin -dil 1 -fill -ero 1 -odt char ${PathID}/${PN}_${TP[tp]}/Brainmask/${PN}_${TP[tp]}_B1.nii.gz" >> ${NameScript}

# Taking care of priors registration
    if [ ! -f ${PathResults}/${PN}_${TP[tp]}_AffTransf.txt ]
    then
		echo "${PathReg}/reg_aladin -ref ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -flo ${PathICBM}/ICBM_Template.nii.gz -aff ${PathScratch}/${PN}_${TP[tp]}_AffTransf.txt" >> ${NameScript}
        echo "cp ${PathScratch}/${PN}_${TP[tp]}_AffTransf.txt ${PathResults}/${PN}_${TP[tp]}_AffTransf.txt " >> ${NameScript}
    fi
    if [ ! -f ${PathResults}/${PN}_${TP[tp]}_cpp.nii.gz ]
    then
		echo "${PathReg}/reg_f3d -ref ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -flo ${PathICBM}/ICBM_Template.nii.gz -aff ${PathScratch}/${PN}_${TP[tp]}_AffTransf.txt -cpp ${PathScratch}/${PN}_${TP[tp]}_cpp.nii.gz" >> ${NameScript}
        echo "cp ${PathScratch}/${PN}_${TP[tp]}_cpp.nii.gz ${PathResults}/${PN}_${TP[tp]}_cpp.nii.gz " >> ${NameScript}
    fi
    for p in CGM DGM ECSF ICSF Out WM
    do
        if [ ! -f ${PathResults}/ICBM_${p}_${PN}_${TP[tp]}.nii.gz ]
        then
            echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -flo ${PathICBM}/ICBM_${p}.nii.gz -cpp ${PathScratch}/${PN}_${TP[tp]}_cpp.nii.gz -res ${PathScratch}/ICBM_${p}_${PN}_${TP[tp]}.nii.gz" >> ${NameScript}
        fi
    done

# summing to obtain CSFs and GM
    if [ ! -f ${PathResults}/*_CSFs_${PN}_${TP[tp]}.nii.gz ]
    then
		echo " ${PathSeg}/seg_maths ${PathScratch}/ICBM_ICSF_${PN}_${TP[tp]}.nii.gz -add ${PathScratch}/ICBM_ECSF_${PN}_${TP[tp]}.nii.gz ${PathScratch}/ICBM_CSFs_${PN}_${TP[tp]}.nii.gz" >> ${NameScript}
    fi
    if [ ! -f ${PathResults}/*_GM_${PN}_${TP[tp]}.nii.gz ]
    then
        echo "${PathSeg}/seg_maths ${PathScratch}/ICBM_CGM_${PN}_${TP[tp]}.nii.gz -add ${PathScratch}/ICBM_DGM_${PN}_${TP[tp]}.nii.gz ${PathScratch}/ICBM_GM_${PN}_${TP[tp]}.nii.gz" >> ${NameScript}
    fi

    array_PriorsICBM=("${PathScratch}/ICBM_GM_${PN}_${TP[tp]}.nii.gz" "${PathScratch}/ICBM_WM_${PN}_${TP[tp]}.nii.gz" "${PathScratch}/ICBM_CSFs_${PN}_${TP[tp]}.nii.gz" "${PathScratch}/ICBM_Out_${PN}_${TP[tp]}.nii.gz")

    array_PriorsGIF=("${PathScratch}/GIF_GM_${PN}_${TP[tp]}.nii.gz" "${PathScratch}/GIF_WM_${PN}_${TP[tp]}.nii.gz" "${PathScratch}/GIF_CSF_${PN}_${TP[tp]}.nii.gz" "${PathScratch}/GIF_Out_${PN}_${TP[tp]}.nii.gz")

    
        PriorsArray=("Out" "CSF" "CGM" "WMI" "DGM" "Brainstem")
        for ((p=0;p<6;p++))
        do
            echo "${PathSeg}/seg_maths ${PathScratch}/GIF_${PN}_${TP[tp]}_prior.nii.gz -tp ${p} ${PathScratch}/GIF_${PriorsArray[p]}_${PN}_${TP[tp]}.nii.gz" >> ${NameScript}
        done
        echo "${PathSeg}/seg_maths ${PathScratch}/ICBM_DGM_${PN}_${TP[tp]}.nii.gz -thr 0.3 -bin -mul ${PathScratch}/GIF_DGM_${PN}_${TP[tp]}.nii.gz ${PathScratch}/GIF_DGM_${PN}_${TP[tp]}_corr.nii.gz" >> ${NameScript}
	echo " ${PathSeg}/seg_maths ${PathScratch}/GIF_DGM_${PN}_${TP[tp]}.nii.gz -sub ${PathScratch}/GIF_DGM_${PN}_${TP[tp]}_corr.nii.gz -thr 0 -add ${PathScratch}/GIF_WMI_${PN}_${TP[tp]}.nii.gz ${PathScratch}/GIF_WMI_${PN}_${TP[tp]}.nii.gz" >> ${NameScript}

        echo "${PathSeg}/seg_maths ${PathScratch}/GIF_CGM_${PN}_${TP[tp]}.nii.gz -add ${PathScratch}/GIF_DGM_${PN}_${TP[tp]}_corr.nii.gz ${PathScratch}/GIF_GM_${PN}_${TP[tp]}.nii.gz" >> ${NameScript}
        echo "${PathSeg}/seg_maths ${PathScratch}/GIF_WMI_${PN}_${TP[tp]}.nii.gz -add ${PathScratch}/GIF_Brainstem_${PN}_${TP[tp]}.nii.gz ${PathScratch}/GIF_WM_${PN}_${TP[tp]}.nii.gz" >> ${NameScript}
        if  [ -f ${GIF_results_path[tp]}/*prior* ]
        then
            array_PriorsToUse=(${array_PriorsGIF[*]})
        else
            array_PriorsToUse=(${array_PriorsICBM[*]})
        fi
    


#if [ ! -f ${PathID}/PN_GW/Priors/${PN}_GW${OptTP}_Parcellation.nii.gz ]
#  then
    echo "cp ${PathID}/${PN}_${TP[tp]}/${PathGIF}/*Parcellation* ${PathScratch}/${PN}_${TP[tp]}_Parcellation.nii.gz ">>${NameScript}
# fi
#101 102 139 140 187 188 169 170
    if [ ! -f ${PathResults}/${PN}_${TP[tp]}_Artefacts.nii.gz ]
    then
        if ((${#ArtefactArray[@]}>0))
        then
            stringAddition="${PathScratch}/${PN}_${TP[tp]}_ArtConstruction_${i}.nii.gz "
            for ((i=0;i<${#ArtefactArray[@]};i++))
            do
                Value=${ArtefactArray[i]}
                ValueMin=`echo "$Value - 0.5" | bc -l`
                ValueMax=`echo "$Value + 0.5" | bc -l`
                echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_${TP[tp]}_Parcellation.nii.gz -equal  $Value -bin ${PathScratch}/${PN}_${TP[tp]}_ArtConstruction_${i}.nii.gz " >> ${NameScript}
                stringAddition="${stringAddition} -add ${PathScratch}/${PN}_${TP[tp]}_ArtConstruction_${i}.nii.gz "
            done
            echo "${PathSeg}/seg_maths ${stringAddition} -bin ${PathScratch}/${PN}_${TP[tp]}_Artefacts.nii.gz " >> ${NameScript}
            echo "cp ${PathScratch}/${PN}_${TP[tp]}_Artefacts.nii.gz ${PathResults}/." >> ${NameScript}
            echo "rm ${PathScratch}/*ArtConstr* " >> ${NameScript}
        fi
    fi






    ProgMod=""
    if ((${#array_ToStack[@]}==2))
    then
        ProgMod="-progMod 1 2 "
    else
        ProgMod="-progMod 2 2 1"
    fi
    # Perform Seg_BiASM EMfO on each
    if [ ! -f ${PathResults}/${NameTogether}_EMfO_${PN}_${TP[tp]}*.nii.gz ]
    then
        echo "${PathSeg}/seg_maths ${PathScratch}/*${PN}_${TP[tp]}*TIV.nii.gz -bin -odt char ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz" >> ${NameScript}
        echo "cp ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz ${PathScratch}/${PN}_${TP[tp]}_GIF_B1.nii.gz" >> ${NameScript}
		echo "${PathSeg}/Seg_BiASM -in ${#array_ToStack[@]} ${array_ToStack[*]} -priors 4 ${array_PriorsToUse[*]} -mask ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz -out 2 ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}.nii.gz ${PathScratch}/${NameTogether}_EMfOf_${PN}_${TP[tp]}.nii.gz -txt_out ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}.txt -bc_order 3 -CovPriors 8 -BFP 1 -maxRunEM 1 ${ProgMod} -AtlasSmoothing 1 1 -AtlasWeight 2 0 1  -SMOrder 0 -KernelSize 3 -PriorsKept 8 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -deleteUW 1 -outliersM 3 -outliersW 0.01 -init_splitUnif 0 -splitAccept 0 -unifTot 1 -MRF 1 -GMRF ${NameGMatrix} -juxtaCorr 1 -meanPriors 0 -priorDGM ${PathScratch}/ICBM_DGM_${PN}_${TP[tp]}.nii.gz -BiASM 0 ${TAString}" >> ${NameScript}
        echo "rm ${PathScratch}/BG* " >> ${NameScript}
        echo "rm ${PathScratch}/MRF* " >> ${NameScript}
        echo "cp ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}* ${PathResults}/." >> ${NameScript}

	fi
# Seg_Analysis for IO

    if [ ! -f ${PathResults}/I_${NameTogether}_EMfO_${PN}_${TP[tp]}${SignModa}.nii.gz ]
    then
		echo "${PathSeg}/Seg_Analysis -inTxt2 ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}${SignModa}.txt ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}${SignModa}.nii.gz -IO 1 -mask ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz ${ChangePathString}" >> ${NameScript}
    fi
	echo "cp ${PathScratch}/*${NameTogether}*EMfO*${PN}_${TP[tp]}* ${PathResults}/." >> ${NameScript}
done


if  [ ! -f ${PathResults}/T1*GW* ] 
then

# Create the file for the group wise
    NameFileGW=${PathResults}/GroupwiseParam_${PN}.sh
    echo "" > ${NameFileGW}
    echo "#!/bin/sh" >> ${NameFileGW}
    echo "export PathSeg=$PathSeg" >> ${NameFileGW}
    echo "export PathReg=$PathReg" >> ${NameFileGW}
    IMG_INPUT_TOPUT=()
    IMG_INPUT_INLIERS_TOPUT=()
    MASK_INPUT_TOPUT=()
    for ((tp=0;tp<${#TP[@]};tp++))
    do
        IMG_INPUT_TOPUT=(${IMG_INPUT_TOPUT[*]} "${PathScratch}/DataCorrected_${NameTogether}_EMfO_${PN}_${TP[tp]}.nii.gz")
        IMG_INPUT_INLIERS_TOPUT=(${IMG_INPUT_INLIERS_TOPUT[*]} "${PathScratch}/I_${NameTogether}_EMfO_${PN}_${TP[tp]}${SignModa}.nii.gz")
        MASK_INPUT_TOPUT=(${MASK_INPUT_TOPUT[*]} "${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz")
    done
    printf "%s\n" "${testa[@]}" > file.txt
    echo "export IMG_INPUT=( "${IMG_INPUT_TOPUT[@]}")" >> ${NameFileGW}
# echo "export IMG_INPUT_MASK=("${MASK_INPUT_TOPUT[@]}")" >> ${NameFileGW}
    echo " export IMG_INPUT_INLIERS=( "${IMG_INPUT_INLIERS_TOPUT[@]}")" >> ${NameFileGW}
# echo "export IMG_INPUT=(`ls ${PathID}/${PN}_*/Results/DataCorrected*${NameTogether}_EMfO_${PN}_*.nii.gz`)" >> ${NameFileGW}
#echo " export IMG_INPUT_INLIERS=(`ls ${PathID}/${PN}_*/Results/I_*${NameTogether}_EMfO_${PN}_*.nii.gz`)" >> ${NameFileGW}

#IMG_INPUT=(`ls ${PathID}/${PN}_*/Results/DataCorrected*${NameTogether}_EMfO_${PN}_*.nii.gz`)
#IMG_INPUT=${IMG_INPUT_TOPUT[*]}
    echo "export TEMPLATE=${IMG_INPUT_TOPUT[0]} " >> ${NameFileGW}
#echo "export TEMPLATE_MASK=${MASK_INPUT_TOPUT[0]}" >> ${NameFileGW}
    echo "export RES_FOLDER="${PathScratch}"" >> ${NameFileGW}
    echo 'export NRR_args="-ln 3 -lp 3 -maxit 100 -vel"' >> ${NameFileGW}
    echo 'export AFF_IT_NUM=3' >> ${NameFileGW}
    echo 'export NRR_IT_NUM=2' >> ${NameFileGW}

#echo 'export QSUB_CMD="qsub -l h_rt=05:00:00 -l tmem=3.5G -l h_vmem=3.5G -l vf=3.5G -l s_stack=10240  -j y -S /bin/csh -b y -cwd -V"'>> ${NameFileGW}

    echo " sh ${ScriptGW} ${NameFileGW} ${NameScriptGroupwise}" >> ${NameScript}
    echo "sh ${NameScriptGroupwise}" >> ${NameScript}
fi
#Once done, get the priors on the new mean model 

#Destack the result
    array_DCGW=()
    NameFinalAverage=${PathScratch}/nrr_2/average_nonrigid_it_2.nii.gz
    for ((i=0;i<${#array_ToStack[@]};i++))
    do
        echo "${PathSeg}/seg_maths ${NameFinalAverage} -tp ${i} ${PathScratch}/${FinalMod[i]}_${PN}_GW${OptTP}.nii.gz" >> ${NameScript}
            array_DCGW=(${array_DCGW[*]} "${PathScratch}/${FinalMod[i]}_${PN}_GW${OptTP}.nii.gz" )
    done

    echo "cp -r ${PathScratch}/nrr_2 ${PathResults}/." >> ${NameScript}

# First do a reg_resample on all time points and all atlases using the final cpp transformation
    echo "${#TP[@]} time points "





    for p in CGM DGM GM WM Out ICSF ECSF CSFs
    do
        if [ ! -f ${PathResults}/ICBM_${p}_${PN}_GW${OptTP}.nii.gz ]
        then
            array_priorspre=()
            for ((tp=0;tp<${#TP[@]};tp++))
            do
                
                    echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_GW${OptTP}.nii.gz -flo  ${PathScratch}/ICBM_${p}_${PN}_${TP[tp]}.nii.gz -cpp  ${PathScratch}/nrr_2/nrr_cpp_DataCorrected_T1FLAIR_EMfO_${PN}_${TP[tp]}_it2.nii.gz -res ${PathScratch}/ICBM_${p}_${PN}_${TP[tp]}_preGW.nii.gz" >> ${NameScript}
                
                array_priorspre=(${array_priorspre[*]} "${PathScratch}/ICBM_${p}_${PN}_${TP[tp]}_preGW.nii.gz")
            done

    # Mean over the transformed time points
            if [ ! -f ${PathResults}/ICBM_${p}_${PN}_GW${OptTP}.nii.gz ]
            then
                echo "${PathSeg}/Seg_Analysis -meanImages ${#TP[@]} ${array_priorspre[*]} -outMean ${PathScratch}/ICBM_${p}_${PN}_GW${OptTP}.nii.gz" >> ${NameScript}
                echo "rm ${PathID}/${PN}_GW${OptTP}/Priors/*preGW*" >> ${NameScript}
            fi
        fi
    done


    if [ ! -f ${PathResults}/MeanT1_${PN}.nii.gz ]
    then
        array_meanpre=();
        for ((tp=0;tp<${#TP[@]};tp++))
        do
            if [ ! -f ${PathResults}/T1_${PN}_${TP[tp]}_preGW.nii.gz ]
            then
                echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_GW${OptTP}.nii.gz -flo  ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -cpp  ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}_it2.nii.gz -res ${PathScratch}/T1_${PN}_${TP[tp]}_preGW.nii.gz" >> ${NameScript}
            fi
            array_meanpre=(${array_meanpre[*]} "${PathScratch}/T1_${PN}_${TP[tp]}_preGW.nii.gz")
        done

        if [ -f ${PathResults}/MeanT1_${PN}.nii.gz ]
        then
            echo "rm ${PathScratch}/*/*log*" >> ${NameScript}
            echo  "rm ${PathScratch}/*/*inliers* " >> ${NameScript}
            echo "rm ${PathScratch}/*/*res* " >> ${NameScript}
            echo "rm ${PathScratch}/*/mask* " >> ${NameScript}
        fi

# Mean over the transformed time points

        echo "${PathSeg}/Seg_Analysis -meanImages ${#TP[@]} ${array_meanpre[*]} -outMean ${PathScratch}/MeanT1_${PN}.nii.gz" >> ${NameScript}
        echo "cp ${PathScratch}/Mean* ${PathResults}/." >> ${NameScript}

        echo "rm -r ${PathScratch}/*preGW* " >> ${NameScript}
    fi





    array_artefacts=()
    for ((tp=0;tp<${#TP[@]};tp++))
    do
        if [ ! -f ${PathScratch}/${PN}_${TP[tp]}_preGW_Artefacts.nii.gz ]
        then
            echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_GW${OptTP}.nii.gz -flo  ${PathScratch}/${PN}_${TP[tp]}_Artefacts.nii.gz -cpp  ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}_it2.nii.gz -res ${PathScratch}/${PN}_${TP[tp]}_preGW_Artefacts.nii.gz" >> ${NameScript}
        fi
        array_artefacts=(${array_artefacts[*]} "${PathScratch}/${PN}_${TP[tp]}_preGW_Artefacts.nii.gz")
    done



    if [ ! -f ${PathResults}/${PN}_GW${OptTP}_Artefacts.nii.gz ]
    then
        echo "${PathSeg}/Seg_Analysis -meanImages ${#TP[@]} ${array_artefacts[*]} -outMean ${PathScratch}/${PN}_GW${OptTP}_Artefacts.nii.gz" >> ${NameScript}
        echo "rm ${PathResults}/*preGW*" >> ${NameScript}
    fi


    for p in CGM DGM GM WM Out CSF
    do
        if [ ! -f ${PathResults}/GIF_${p}_${PN}_GW${OptTP}.nii.gz ]
        then
            array_priorspre=()
            for ((tp=0;tp<${#TP[@]};tp++))
            do
                if [ ! -f ${PathScratch}/ICBM_${p}_${PN}_${TP[tp]}_preGW.nii.gz ]
                then
                    echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_GW${OptTP}.nii.gz -flo  ${PathScratch}/GIF_${p}_${PN}_${TP[tp]}.nii.gz -cpp  ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}_it2.nii.gz -res ${PathScratch}/GIF_${p}_${PN}_${TP[tp]}_preGW.nii.gz" >> ${NameScript}
                fi
                array_priorspre=(${array_priorspre[*]} "${PathScratch}/GIF_${p}_${PN}_${TP[tp]}_preGW.nii.gz")
            done

# Mean over the transformed time points

            echo "${PathSeg}/Seg_Analysis -meanImages ${#TP[@]} ${array_priorspre[*]} -outMean ${PathScratch}/GIF_${p}_${PN}_GW${OptTP}.nii.gz" >> ${NameScript}
        fi
    done

    array_brainmask=()
    for ((tp=0;tp<${#TP[@]};tp++))
    do
        if [ ! -f ${PathResults}/GIF_B1_${PN}_${TP[tp]}_preGW.nii.gz ]
        then
            echo "cp ${PathScratch}/${PN}_${TP[tp]}_GIF_TIV.nii.gz ${PathScratch}/${PN}_${TP[tp]}_GIF_B1.nii.gz" >> ${NameScript}
            echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_GW${OptTP}.nii.gz -flo  ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz -cpp  ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}_it2.nii.gz -res ${PathScratch}/GIF_B1_${PN}_${TP[tp]}_preGW.nii.gz" >> ${NameScript}
        fi
        array_brainmask=(${array_brainmask[*]} "${PathScratch}/GIF_B1_${PN}_${TP[tp]}_preGW.nii.gz")
    done

# Mean over the transformed time points

    echo "${PathSeg}/Seg_Analysis -meanImages ${#TP[@]} ${array_brainmask[*]} -outMean ${PathScratch}/GIF_B1_${PN}_GW${OptTP}.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_B1_${PN}_GW${OptTP}.nii.gz -thr 0.5 -bin -odt char ${PathScratch}/${PN}_GW${OptTP}_B1.nii.gz" >> ${NameScript}
    echo "cp ${PathScratch}/GIF_B1_${PN}_GW${OptTP}.nii.gz ${PathResults}/${PN}_GW${OptTP}_GIF_B1.nii.gz" >> ${NameScript}
    echo "cp ${PathScratch}/GIF_B1_${PN}_GW${OptTP}.nii.gz ${PathScratch}/${PN}_GW${OptTP}_GIF_B1.nii.gz" >> ${NameScript}
    array_brainmask=()








# Model selection
 echo "We have to do MS seg"
if (($DO_LES==1))
then
array_DCGW=()
NameFinalAverage=${PathScratch}/nrr_2/average_nonrigid_it_2.nii.gz
for ((i=0;i<${#array_ToStack[@]};i++))
do
    echo "${PathSeg}/seg_maths ${NameFinalAverage} -tp ${i} ${PathScratch}/${FinalMod[i]}_${PN}_GW${OptTP}.nii.gz" >> ${NameScript}
    array_DCGW=(${array_DCGW[*]} "${PathScratch}/${FinalMod[i]}_${PN}_GW${OptTP}.nii.gz" )
done

array_PriorsICBM=("${PathScratch}/ICBM_GM_${PN}_GW${OptTP}.nii.gz" "${PathScratch}/ICBM_WM_${PN}_GW${OptTP}.nii.gz" "${PathScratch}/ICBM_CSFs_${PN}_GW${OptTP}.nii.gz" "${PathScratch}/ICBM_Out_${PN}_GW${OptTP}.nii.gz")

array_PriorsGIF=("${PathScratch}/GIF_GM_${PN}_GW${OptTP}.nii.gz" "${PathScratch}/GIF_WM_${PN}_GW${OptTP}.nii.gz" "${PathScratch}/GIF_CSF_${PN}_GW${OptTP}.nii.gz" "${PathScratch}/GIF_Out_${PN}_GW${OptTP}.nii.gz")

NameBrainmaskGW=""

if [ -f ${PathResults}/Mean*prior* ]
then
    PriorsArray=("Out" "CSF" "CGM" "WMI" "DGM" "Brainstem")
    echo "${PathSeg}/seg_maths ${PathScratch}/Mean*TIV.nii.gz -thr 0.5 -bin -odt char ${PathScratch}/${PN}_GW${OptTP}_GIF_B1.nii.gz" >> ${NameScript}
    for ((p=0;p<6;p++))
    do
        echo "${PathSeg}/seg_maths ${PathScratch}/Mean*prior* -tp ${p} ${PathScratch}/GIF_${PriorsArray[p]}_${PN}_GW${OptTP}.nii.gz" >> ${NameScript}
    done
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_CGM_${PN}_GW${OptTP}.nii.gz -add ${PathScratch}/GIF_DGM_${PN}_GW${OptTP}.nii.gz ${PathScratch}/GIF_GM_${PN}_GW${OptTP}.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_WMI_${PN}_GW${OptTP}.nii.gz -add ${PathScratch}/GIF_Brainstem_${PN}_GW${OptTP}.nii.gz ${PathScratch}/GIF_WM_${PN}_GW${OptTP}.nii.gz" >> ${NameScript}
    array_PriorsGIF=("${PathScratch}/GIF_GM_${PN}_GW${OptTP}.nii.gz" "${PathScratch}/GIF_WM_${PN}_GW${OptTP}.nii.gz" "${PathScratch}/GIF_CSF_${PN}_GW${OptTP}.nii.gz" "${PathScratch}/GIF_Out_${PN}_GW${OptTP}.nii.gz")
    array_PriorsToUse=(${array_PriorsGIF[*]})
    NameBrainmaskGW=${PathScratch}/${PN}_GW${OptTP}_GIF_B1.nii.gz
fi
if [ -f ${PathResults}/GIF_GM_${PN}_${TP[0]}.nii.gz ]
then
    array_PriorsToUse=(${array_PriorsGIF[*]})
    NameBrainmaskGW=${PathScratch}/${PN}_GW${OptTP}_GIF_B1.nii.gz
elif [ -f ${PathResults}/GIF_GM_${PN}_GW${OptTP}.nii.gz ] || [ -f ${GIF_results_path[0]}/*prior* ]
then
    array_PriorsToUse=(${array_PriorsGIF[*]})
    NameBrainmaskGW=${PathScratch}/${PN}_GW${OptTP}_GIF_B1.nii.gz
else
    array_PriorsToUse=(${array_PriorsICBM[*]})
    NameBrainmaskGW=${PathScratch}/${PN}_GW${OptTP}_B1.nii.gz
fi
    # else
    #     array_PriorsToUse=(${array_PriorsICBM[*]})
    # fi

NameBrainmaskGW=${PathScratch}/${PN}_GW${OptTP}_GIF_B1.nii.gz

if [ ! -f ${PathResults}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.txt ]
then
        # if ((Do2==0))
        # then
     echo "Do2=0"
    echo "${PathSeg}/Seg_BiASM -in ${#array_DCGW[@]} ${array_DCGW[*]} -priors 4 ${array_PriorsToUse[*]} -mask ${NameBrainmaskGW} -out 2 ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz ${PathScratch}/${NameTogether}_BiASMG_${PN}_GW${OptTP}.nii.gz -txt_out ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.txt -bc_order 0 -inDC_flag 1 -CovPriors 8 -BFP 1 -maxRunEM ${MRE} ${ProgMod} -AtlasSmoothing 1 1 -AtlasWeight 2 ${AW} 1  -SMOrder 0 -KernelSize 3 -PriorsKept 8 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -deleteUW 1 -outliersM 3 -outliersW ${OW} -init_splitUnif 0 -splitAccept 0 -unifTot 1 -MRF 1 -GMRF ${NameGMatrix} -juxtaCorr 1 -priorDGM ${PathScratch}/ICBM_DGM_${PN}_GW${OptTP}.nii.gz ${TAString}" >> ${NameScript}

        # else
        #     if [ ! -f ${PathID}/${PN}_GW${OptTP}/Results/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.txt ]
        #     then
        #         echo "Segmentation 2 to do"
        #         echo "${PathSeg}/Seg_BiASM -in 2 ${array_DCGW[0]} ${array_DCGW[1]} -priors 4 ${array_PriorsToUse[*]} -mask ${NameBrainmaskGW} -out 2 ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz ${PathScratch}/${NameTogether}_BiASMG_${PN}_GW${OptTP}.nii.gz -txt_out ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.txt -bc_order 0 -inDC_flag 1 -CovPriors 8 -BFP 1 -maxRunEM ${MRE} -AtlasSmoothing 1 1 -AtlasWeight 2 ${AW} 1  -SMOrder 0 -KernelSize 3 -PriorsKept 8 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -deleteUW 1 -outliersM 3 -outliersW ${OW} -init_splitUnif 0 -splitAccept 0 -unifTot 1 -MRF 1 -GMRF ${NameGMatrix} -juxtaCorr 1 -progMod 0 -priorDGM ${PathScratch}/ICBM_DGM_${PN}_GW${OptTP}.nii.gz ${TAString} " >> ${NameScript}
        #     fi

        # fi
        #fi
fi
#fi

echo "${PathSeg}/Seg_Analysis -inTxt2 ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.txt ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz -mask ${NameBrainmaskGW} -Package 1 -SegType 1 -WeightedSeg 3 3 1 -connect -correct -inModa 2 ${ArrayModaCode[0]} ${ArrayModaCode[1]}  -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${PathScratch}/ICBM_ICSF_${PN}_GW${OptTP}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${PN}_GW${OptTP}.nii.gz -TO 1 -juxtaCorr 1 ${ChangePathString}" >> ${NameScript}

# EMf using summarised Corrected for atlases at each time point
#read -a NameSummarisedCorr <<< $(ls ${PathID}/${PN}_GW${OptTP}/Results/SummarisedCorr*${NameTogether}_BiASM_${PN}_GW${OptTP}.nii.gz)
NameSummarisedCorr=${PathScratch}/SummarisedCorrected_WS3WT3WC1JC1ST1CL2CIV1_TO${NameTogether}_BiASM_${PN}_GW${OptTP}.nii.gz


# Obtention of final segmentation


if [ ! -f ${PathResults}/SummarisedCorr*${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}* ] & [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/${NameTogether}_BiASM_${PN}_${TP[tp]}_${SignModa}*HM${OptTP}${TAOpt}.nii.gz ]
then
    echo "${PathSeg}/Seg_Analysis -inTxt2 ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.txt ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz -mask ${PathScratch}/${PN}_GW${OptTP}_B1.nii.gz -Package 1 -SegType 1 -WeightedSeg 3 3 1 -connect -correct -inModa 2 ${ArrayModaCode[0]} ${ArrayModaCode[1]}  -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${PathScratch}/ICBM_ICSF_${PN}_GW${OptTP}.nii.gz -inPriorsECSF ${PathScratch}/ICBM_ECSF_${PN}_GW${OptTP}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${PN}_GW${OptTP}.nii.gz -inPriorsCGM ${PathScratch}/ICBM_CGM_${PN}_GW${OptTP}.nii.gz -TO 1 -juxtaCorr 1 ${ChangePathString}" >> ${NameScript}
fi

if [ ! -f ${PathResults}/IO_${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.nii.gz ]
then
    echo "${PathSeg}/Seg_Analysis -inTxt2 ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.txt ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.nii.gz -IO 1 -mask ${NameBrainmaskGW} ${ChangePathString}" >> ${NameScript}
fi


NameDataCorrectedGW=${PathScratch}/DataCorrected_${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz

# EMf using summarised Corrected for atlases at each time point
#read -a NameSummarisedCorr <<< $(ls ${PathID}/${PN}_GW${OptTP}/Results/SummarisedCorr*${NameTogether}_BiASM_${PN}_GW${OptTP}${SignModa}.nii.gz)
NameSummarisedCorr=${PathScratch}/SummarisedCorrected_WS3WT3WC1JC1SP1ST1CL2CIV1_TO${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.nii.gz
echo "Sum corr to use is ${NameSummarisedCorr}"
array_Brainmask=()
for ((tp=0;tp<${#TP[@]};tp++))
do

#Destack the DataCorrected result
#NameDC=()
#read -a NameDC <<< $(ls ${PathID}/${PN}_${TP[tp]}/Results/DataCorrected*EMfO*)
    NameDC=${PathScratch}/DataCorrected_${NameTogether}_EMfO_${PN}_${TP[tp]}.nii.gz
    array_DC=()

    for ((i=0;i<${#array_ToStack[@]};i++))
    do
        echo "${PathSeg}/seg_maths ${NameDC} -tp ${i} ${PathScratch}/${FinalMod[i]}_${PN}_${TP[tp]}_DC.nii.gz" >> ${NameScript}
        array_DC=(${array_DC[*]} "${PathScratch}/${FinalMod[i]}_${PN}_${TP[tp]}_DC.nii.gz" )
    done
	

# resample SummarisedCorr into initial space
    echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -flo ${NameSummarisedCorr} -cpp ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}_it2*backward*.nii.gz -res  ${PathScratch}/GC_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz " >> ${NameScript}
    echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -flo ${PathScratch}/IO_${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.nii.gz -cpp ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}_it2*backward*.nii.gz -res  ${PathScratch}/IO_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz -inter 1" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/IO_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz -tp 0 ${PathScratch}/I_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz ">> ${NameScript}
    echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -flo ${PathScratch}/${PN}_GW${OptTP}_GIF_B1.nii.gz -cpp ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}_it2*backward*.nii.gz -res  ${PathScratch}/${PN}_${TP[tp]}_GIF.nii.gz " >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz -bin -fill -odt char ${PathScratch}/${PN}_${TP[tp]}_GIF_B1.nii.gz " >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz -bin -odt char ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz"  >> ${NameScript}  
array_Brainmask=("${array_Brainmask[*]}" ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz )
            # else
            #     array_Brainmask=("${array_Brainmask[*]}" ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz )
            # fi


# Destack it into priors
    array_Priors=()
    for i in 0 1 2 3
    do
        echo "${PathSeg}/seg_maths ${PathScratch}/GC_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz -tp ${i} ${PathScratch}/GC_${FinalPriors[i]}_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz " >> ${NameScript}
        array_Priors=(${array_Priors[*]} "${PathScratch}/GC_${FinalPriors[i]}_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz ")
                    #echo ${PathScratch}
    done

    echo ${SignModa} "Having to do the EMf with new atlases and progressive inclusion but no BF ?"

# # EMf using these atlases with progressive and no BF
#                 if [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.nii.gz ] & [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/${NameTogether}_BiASM_${PN}_${TP[tp]}_${TAOpt}${SignModa}*HM${OptTP}${TAOpt}.nii.gz ]
# 				then 
#                     echo "Yes we have to !" ${array_Brainmask[tp]} ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz 
#                     echo "${PathSeg}/Seg_BiASM -in ${#array_DC[@]} ${array_DC[*]} -priors 4 ${array_Priors[*]} -mask ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz -out 2 ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}.nii.gz ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}.nii.gz -txt_out ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}.txt -bc_order 0 -CovPriors 8 -BFP 1 -maxRunEM 1 ${ProgMod} -AtlasSmoothing 1 1 -AtlasWeight 2 ${AW} 1  -SMOrder 0 -KernelSize 3 -PriorsKept 8 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -deleteUW 1 -outliersM 3 -outliersW 0.01 -init_splitUnif 0 -splitAccept 0 -unifTot 1 -MRF 1 -GMRF ${NameGMatrix} -juxtaCorr 1 -BiASM 0 -inDC_flag 1 -priorDGM ${PathScratch}/ICBM_DGM_${PN}_${TP[tp]}.nii.gz -BiASM 0 -maxRunEM 1 -juxtaCorr 1 ${TAString}" >> ${NameScript}
#                     echo "rm ${PathScratch}/BG*" >> ${NameScript}
#                     echo "rm ${PathScratch}/GC_*" >> ${NameScript}
#                 fi

# # Get IO atlases from that one
#                 if [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/IO_${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${SignModa}.nii.gz ] & [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/${NameTogether}_BiASM_${PN}_${TP[tp]}_${TAOpt}${SignModa}*HM${OptTP}${TAOpt}.nii.gz ]
#                 then
#                     echo "${PathSeg}/Seg_Analysis -inTxt2 ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.txt ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.nii.gz -IO 1 -mask ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz ${ChangePathString}" >> ${NameScript}
#                 fi
    echo "${PathReg}/reg_resample -ref ${PathScratch}/DataCorrected_${NameTogether}_EMfO_${PN}_${TP[tp]}.nii.gz -flo ${PathScratch}/DataCorrected_${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz -cpp ${PathScratch}/nrr_2/*cpp*${NameTogether}*${PN}_${TP[tp]}_it2*backward*.nii.gz -res  ${PathScratch}/DataCorrected_GWto${TP[tp]}_${NameTogether}_${PN}.nii.gz " >> ${NameScript}
           
# echo "${PathSeg}/Seg_BiASM -averageGC 0.5 ${PathScratch}/GC_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz  ${PathScratch}/GC_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz -averageIO 0.5 ${PathScratch}/IO_${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.nii.gz ${PathScratch}/IO_${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.nii.gz -inDC ${PathScratch}/DataCorrected_HM_GWto${TP[tp]}_${NameTogether}_${PN}_0.nii.gz -txt_in ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.txt -meanPriors 1 -txt_out ${PathScratch}/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.txt -out 2 ${PathScratch}/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.nii.gz ${PathScratch}/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.nii.gz -meanPriors 1 -juxtaCorr 1 -mask ${NameBrainmask} -MRF 1 -GMRF ${NameGMatrix} -bc_order 0 -CovPriors 8 -priorDGM ${PathScratch}/ICBM_DGM_${PN}_${TP[tp]}.nii.gz -BiASM 0 -maxRunEM 1 -outliersM 3 -PriorsKept 5 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -AtlasSmoothing 1 0.3" >> ${NameScript}
NameBrainmask=${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz
if [ ! -f ${PathResults}/${NameTogether}*BiASM*HM${OptTP}${TAOpt}*.nii.gz ]
then
    NameBrainmask=${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz
    echo "${PathSeg}/Seg_Analysis -mask ${NameBrainmask} -maskMatch ${PathScratch}/I_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz -orderFit 1 -matchRef 1 ${PathScratch}/DataCorrected_GWto${TP[tp]}_${NameTogether}_${PN}.nii.gz -matchFloat 1 ${PathScratch}/DataCorrected_${NameTogether}_EMfO_${PN}_${TP[tp]}.nii.gz -outFit ${PathScratch}/DataCorrected_HM_GWto${TP[tp]}_${NameTogether}_${PN}.nii.gz" >> ${NameScript}
    echo "${PathSeg}/Seg_BiASM -averageGC 0.5 ${PathScratch}/GC_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz  ${PathScratch}/GC_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz -averageIO 0.5 ${PathScratch}/IO_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz ${PathScratch}/IO_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz -inDC ${PathScratch}/DataCorrected_HM_GWto${TP[tp]}_${NameTogether}_${PN}_0.nii.gz -txt_in ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.txt -meanPriors 1 -txt_out ${PathScratch}/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}${SignModa}.txt -out 2 ${PathScratch}/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}${SignModa}.nii.gz ${PathScratch}/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}${SignModa}.nii.gz -meanPriors 1 -juxtaCorr 1 -mask ${NameBrainmask} -MRF 1 -GMRF ${NameGMatrix} -bc_order 0 -CovPriors 8 -priorDGM ${PathScratch}/ICBM_DGM_${PN}_${TP[tp]}.nii.gz -BiASM 0 -maxRunEM 1 -outliersM 3 -PriorsKept 5 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -AtlasSmoothing 1 0.3" >> ${NameScript}
    fi        

            # Segmentation and correction // Cross sectional version

    Opt=HM${OptTP}${TAOpt}${SignModa}

    if ((${#ArtefactArray[@]}>0))
    then
        stringAddition="${PathScratch}/${PN}_${TP[tp]}_ArtConstruction_0.nii.gz "
        for ((i=0;i<${#ArtefactArray[@]};i++))
        do
            Value=${ArtefactArray[i]}
            ValueMin=`echo "$Value - 0.5" | bc -l`
            ValueMax=`echo "$Value + 0.5" | bc -l`
            echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${PN}_${TP[tp]}.nii.gz -equal $Value -bin ${PathScratch}/${PN}_${TP[tp]}_ArtConstruction_${i}.nii.gz " >> ${NameScript}
            stringAddition="${stringAddition} -add ${PathScratch}/${PN}_${TP[tp]}_ArtConstruction_${i}.nii.gz "
        done
        echo "${PathSeg}/seg_maths ${stringAddition} -bin ${PathScratch}/${PN}_${TP[tp]}_Artefacts.nii.gz" >> ${NameScript}
    fi
            #${ID}_${TP[tp]}

    echo "echo "Lesion segmentation"" >> ${NameScript}
if [ ! -f ${PathResults}/PrimaryLesion*${ID}_${TP[tp]}* ]
then
    echo "${PathSeg}/Seg_Analysis  -inTxt2 ${PathScratch}/${ModalitiesTot}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}${SignModa}.txt ${PathScratch}/${ModalitiesTot}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}${SignModa}.nii.gz -mask ${PathScratch}/GIF_${ID}_${TP[tp]}_B1.nii.gz -Package 1 -SegType 1 -WeightedSeg 3 3 1 -connect -correct -inModa 2 1 3 -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${PathScratch}/ICBM_ICSF_${ID}_${TP[tp]}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${ID}_${TP[tp]}.nii.gz -inPriorsCGM ${PathScratch}/ICBM_CGM_${ID}_${TP[tp]}.nii.gz -inPriorsECSF ${PathScratch}/ICBM_ECSF_${ID}_${TP[tp]}.nii.gz -TO 1 -ParcellationIn ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz -typeVentrSeg ${TVS} -Simple -Secondary 60 -juxtaCorr ${JC}  -SP ${OptSP} -LevelCorrection ${OptCL} ${ChangePathString} -LesWMI ${OptWMI} -Neigh 18" >> ${NameScript}
fi

echo "echo "cp ${PathScratch}/LesionCorrected*WS3WT3WC1*${ID}_${TP[tp]}* ${PathScratch}/PrimaryLesions_${ID}_${TP[tp]}.nii.gz"" >> ${NameScript}
 echo "echo  " cp ${PathScratch}/SecondaryCorrected*WS3WT3WC1*${ID}_${TP[tp]}* ${PathScratch}/SecondaryLesions_${ID}_${TP[tp]}.nii.gz "" >> ${NameScript}
    echo "cp ${PathScratch}/LesionCorrected*WS3WT3WC1*${ID}_${TP[tp]}* ${PathScratch}/PrimaryLesions_${ID}_${TP[tp]}.nii.gz" >> ${NameScript}
    echo "cp ${PathScratch}/SecondaryCorrected*WS3WT3WC1*${ID}_${TP[tp]}* ${PathScratch}/SecondaryLesions_${ID}_${TP[tp]}.nii.gz" >> ${NameScript}

    echo "cp ${PathScratch}/Primary* ${PathResults}/. " >> ${NameScript}
    echo "cp ${PathScratch}/Secondary* ${PathResults}/. " >> ${NameScript}

    echo "rm ${PathScratch}/DataR* ${PathScratch}/LesionT* ${PathScratch}/Summ*${ID}_${TP[tp]}*" >> ${NameScript}

    echo "echo "Correction for septum pellucidum"" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz -thr 65.5 -uthr 67.5 -bin ${PathScratch}/VentricleLining.nii.gz " >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz  -equal 52 -bin -euc -uthr 0 -abs ${PathScratch}/Temp1.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz  -equal 53 -bin -euc -uthr 0 -abs -add ${PathScratch}/Temp1.nii.gz -uthr 5 -bin ${PathScratch}/PotentialInterVentr.nii.gz" >> ${NameScript}

    echo "${PathSeg}/seg_maths ${PathScratch}/VentricleLining.nii.gz -dil 1 ${PathScratch}/ExpandedVentrLin.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz  -thr 49.5 -uthr 53.5 -bin -sub ${PathScratch}/ExpandedVentrLin.nii.gz -thr 0 -bin ${PathScratch}/PotentialCP.nii.gz" >> ${NameScript}

    #${PathSeg}/seg_maths ${GIF_results_path}/*Parcellation*  -thr 122.5 -uthr 124.5 -bin ${PathScratch}/Temp4.nii.gz

    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Segmentation_${ID}_${TP[tp]}.nii.gz -tp 4 -thr 0.5 -bin -sub ${PathScratch}/VentricleLining* -thr 0 ${PathScratch}/SegDGM.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz  -equal 87 -mul -1 -add ${PathScratch}/PotentialInterVentr.nii.gz -sub ${PathScratch}/SegDGM.nii.gz -thr 0 ${PathScratch}/PotentialSP3.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths  ${PathScratch}/PotentialSP3.nii.gz -add ${PathScratch}/PotentialCP.nii.gz ${PathScratch}/PotentialSPCP.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/ICBM_ICSF_${ID}_${TP[tp]}.nii.gz -thr 0.3 -bin -mul ${PathScratch}/PotentialSPCP.nii.gz -add ${PathScratch}/PotentialSP3.nii.gz -thr 0.2  -bin  ${PathScratch}/PotentialSPCP2.nii.gz" >> ${NameScript}

    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz -thr 83.5 -uthr 84.5 -bin -dil 5 ${PathScratch}/WMDil.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz -thr 91.5 -uthr 92.5 -bin -dil 5 ${PathScratch}/WMDil2.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/PotentialSPCP2.nii.gz -sub ${PathScratch}/WMDil.nii.gz -sub ${PathScratch}/WMDil2.nii.gz -thr 0 ${PathScratch}/PotentialSPCP3.nii.gz " >> ${NameScript}

    echo "${PathSeg}/seg_maths ${PathScratch}/PrimaryLesions_*${ID}_${TP[tp]}* -sub ${PathScratch}/PotentialSPCP3.nii.gz -thr 0 ${PathScratch}/PrimarySPCP_${ID}_${TP[tp]}.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/SecondaryLesions_*${ID}_${TP[tp]}* -sub ${PathScratch}/PotentialSPCP3.nii.gz -thr 0 ${PathScratch}/SecondarySPCP_${ID}_${TP[tp]}.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/PrimarySPCP*${ID}_${TP[tp]}* -merge 1 4 ${PathScratch}/SecondarySPCP*${ID}_${TP[tp]}* -tmax ${PathScratch}/MergedLesion_${ID}_${TP[tp]}_${OptTP}.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/MergedLesion_${ID}_${TP[tp]}_${OptTP}.nii.gz -sub ${PathScratch}/${ID}_${TP[tp]}_Artefacts.nii.gz -thr 0 ${PathScratch}/MergedLesion_${ID}_${TP[tp]}_SPCP_${OptTP}.nii.gz" >> ${NameScript}


    echo "echo "Refined correction for third ventricle"" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz  -sub ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz ${PathScratch}/Artefacts_${ID}_${TP[tp]}.nii.gz " >> ${NameScript}

    for ((i=0;i<${#ListCorr[@]};i++))
    do
        echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz  -equal ${ListCorr[i]} -bin -add ${PathScratch}/Artefacts_${ID}_${TP[tp]}.nii.gz " >> ${NameScript}
    done
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz  -thr 65.5 -uthr 67.5 ${PathScratch}/VentrLin.nii.gz " >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz  -equal 87 -add ${PathScratch}/VentrLin.nii.gz ${PathScratch}/VentrLin.nii.gz " >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz  -equal 47 -euc  -abs -uthr 5 -bin -mul ${PathScratch}/VentrLin.nii.gz ${PathScratch}/VentrLinSP.nii.gz " >> ${NameScript}

    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz  -equal 52 -bin -euc  -abs  ${PathScratch}//Ventr1.nii.gz " >> ${NameScript}

    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz  -equal 53 -bin -euc -abs -uthr 5 -sub  ${PathScratch}/Ventr1.nii.gz -thr -5 -uthr 5 -bin -mul  ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz  -equal 52 -bin -euc -abs -uthr 5 -mul ${PathScratch}/VentrLin.nii.gz  -add ${PathScratch}/VentrLinSP.nii.gz ${PathScratch}/VentrLinSP.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/VentrLinSP.nii.gz -add ${PathScratch}/Artefacts_${ID}_${TP[tp]}.nii.gz -mul -1 -add ${PathScratch}/MergedLesion_${ID}_${TP[tp]}_SPCP_${OptTP}.nii.gz -thr 0 ${PathScratch}/MergedLesion_${ID}_${TP[tp]}_${OptTP}_corr.nii.gz " >> ${NameScript}

    echo "${PathSeg}/Seg_Analysis -LesWMI ${OptWMI} ${OptInfarcts} -inLesCorr ${PathScratch}/MergedLesion_${ID}_${TP[tp]}_${OptTP}_corr.nii.gz  -inTxt2 ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${TP[tp]}_${Opt}.txt ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${TP[tp]}_${Opt}.nii.gz -mask ${PathScratch}/GIF_${ID}_${TP[tp]}_B1.nii.gz -Package 1 -SegType 1 -WeightedSeg 3 3 1 -connect -correct -inModa ${#ArrayModaCode[@]} ${ArrayModaCode[*]} -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${PathScratch}/ICBM_ICSF_${ID}_${TP[tp]}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${ID}_${TP[tp]}.nii.gz -inPriorsCGM ${PathScratch}/ICBM_CGM_${ID}_${TP[tp]}.nii.gz -inPriorsECSF ${PathScratch}/ICBM_ECSF_${ID}_${TP[tp]}.nii.gz -TO 1 -juxtaCorr 1 -SP ${OptSP} -LevelCorrection ${OptCL} -inArtefact ${PathScratch}/${ID}_${TP[tp]}_Artefacts.nii.gz ${ChangePathString} -ParcellationIn ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz -typeVentrSeg 1 -outWM 1 -outConnect 1 -Neigh 6" >> ${NameScript}

    echo "rm ${PathScratch}/LesionWeigh* ${PathScratch}/Binary* ${PathScratch}/WMDil* ${PathScratch}/WMCard*  ${PathScratch}/LesionInit* ${PathScratch}/DataR* ${PathScratch}/DataT*  ${PathScratch}/LesSegHard* ${PathScratch}/Check* ${PathScratch}/BinaryNIV*" >> ${NameScript}



    # ${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -equal 52 -bin -euc  -abs  ${PathScratch}//Ventr1.nii.gz 

    # ${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -equal 53 -bin -euc -abs -uthr 5 -sub  ${PathScratch}/Ventr1.nii.gz -thr -5 -uthr 5 -bin -mul  ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -equal 52 -bin -euc -abs -uthr 5 -mul ${PathScratch}/VentrLin.nii.gz  -add ${PathScratch}/VentrLinSP.nii.gz ${PathScratch}/VentrLinSP.nii.gz
    # ${PathSeg}/seg_maths ${PathScratch}/VentrLinSP.nii.gz -add ${PathScratch}/Artefacts_${ID}.nii.gz -mul -1 -add ${PathScratch}/MergedLesion_${ID}_SPCP.nii.gz -thr 0 ${PathScratch}/MergedLesion_${ID}_corr.nii.gz 

    # echo "${PathSeg}/seg_maths ${PathScratch}/VentrLinSP.nii.gz -add ${PathScratch}/Artefacts_${ID}.nii.gz -mul -1 -add ${PathScratch}/SecondarySPCP_${ID}.nii.gz -thr 0 ${PathScratch}/SecondarySPCP_${ID}_corr.nii.gz " >> ${NameScript}

    # echo "${PathSeg}/seg_maths ${PathScratch}/VentrLinSP.nii.gz -add ${PathScratch}/Artefacts_${ID}.nii.gz -mul -1 -add ${PathScratch}/PrimarySPCP_${ID}.nii.gz -thr 0 ${PathScratch}/PrimarySPCP_${ID}_corr.nii.gz " >> ${NameScript}

    # echo "${PathSeg}/seg_maths ${PathScratch}/PrimarySPCP_${ID}_corr.nii.gz -merge 1 4 ${PathScratch}/SecondarySPCP_${ID}_corr.nii.gz -tmax ${PathScratch}/MergedLesionSPCP_${ID}_corr.nii.gz" >> ${NameScript}

    # echo "echo "Optional second level of correction"" >> ${NameScript}
    # if ((flag_furtherCorr==1))
    # then
    #     echo "${PathSeg}/Seg_Analysis -inArtefact ${PathScratch}/Artefacts_${ID}.nii.gz -inTxt2 ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${Opt}.txt ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${Opt}.nii.gz -correct -connect -inLesCorr ${PathScratch}/MergedLesionSPCP_${ID}_corr.nii.gz  -mask ${PathScratch}/GIF_${ID}_${TP[tp]}_B1.nii.gz -inPriorsICSF ${PathScratch}/ICBM_ICSF_${ID}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${ID}.nii.gz -inPriorsCGM ${PathScratch}/ICBM_CGM_${ID}.nii.gz -inPriorsECSF ${PathScratch}/ICBM_ECSF_${ID}.nii.gz -outConnect 1 -outWM 1 -ParcellationIn ${PathScratch}/GIF_Parcellation_${ID}.nii.gz  -inChangePath ${PathScratch}/ -inModa 2 1 3 -WeightedSeg 3 3 1 -LesWMI ${OptWMI} -Neigh 6" >> ${NameScript}

    # fi

    if ((flag_corrDGM==1))
    then
        Array=(24   31  32  56  57  58  59  60  61  76  77  37  38)

    # Get segmentation of DGM from GIF

        if [ ! -f ${PathScratch}/${ID}_DGM.nii.gz ]
        then
            if ((${#Array[@]}>0))
            then
                stringAddition="${PathScratch}/${ID}_${TP[tp]}_DGMConstruction_0.nii.gz "
                for ((k=0;k<${#Array[@]};k++))
                do
                    Value=${Array[k]}
                    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz -equal ${Value} -bin ${PathScratch}/${ID}_${TP[tp]}_DGMConstruction_${k}.nii.gz " >> ${NameScript}
                    stringAddition="${stringAddition} -add ${PathScratch}/${ID}_${TP[tp]}_DGMConstruction_${k}.nii.gz "
                done
                echo "${PathSeg}/seg_maths ${stringAddition} -bin -odt char ${PathScratch}/${ID}_${TP[tp]}_DGM.nii.gz " >> ${NameScript}
                echo "rm ${PathScratch}/*DGMConstruction* " >> ${NameScript}
            fi
        fi
        if [ ! -f ${PathScratch}/WMminIn* ]
        then
            echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz -thr 81.5 -uthr 83.5 ${PathScratch}/Insula_${ID}_${TP[tp]}.nii.gz " >> ${NameScript}
            echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz -thr 89.5 -uthr 91.5 -add ${PathScratch}/Insula_${ID}_${TP[tp]}.nii.gz -bin  ${PathScratch}/Insula_${ID}_${TP[tp]}.nii.gz " >> ${NameScript}
            echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz -thr 95.5 -uthr 97.5 -add ${PathScratch}/Insula_${ID}_${TP[tp]}.nii.gz -bin  ${PathScratch}/Insula_${ID}_${TP[tp]}.nii.gz " >> ${NameScript}
            echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz -thr 79 -uthr 98 -bin -sub ${PathScratch}/Insula_${ID}_${TP[tp]}.nii.gz ${PathScratch}/WMminIn_${ID}_${TP[tp]}.nii.gz" >> ${NameScript}
        fi
        echo "${PathSeg}/seg_maths ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz -thr 23 -uthr 45 -bin ${PathScratch}/InfraDGM_${ID}_${TP[tp]}.nii.gz" >> ${NameScript}
        echo "${PathSeg}/seg_maths ${PathScratch}/${ID}_${TP[tp]}_DGM* -add ${PathScratch}/Infra*${ID}_${TP[tp]}* -add ${PathScratch}/Insula*${ID}_${TP[tp]}* -bin -dil 1 -sub ${PathScratch}/WMminIn_${ID}_${TP[tp]}.nii.gz -thr 0 ${PathScratch}/Correction_${ID}_${TP[tp]}.nii.gz" >> ${NameScript}

        echo "${PathSeg}/seg_maths ${PathScratch}/LesionMahal*${ID}_${TP[tp]}* -tp 2 -uthr 4 -bin -mul ${PathScratch}/Corr*Mer*${ID}_${TP[tp]}* -mul ${PathScratch}/Correction_${ID}_${TP[tp]}.nii.gz -mul -1 -add ${PathScratch}/Corr*Mer*${ID}_${TP[tp]}* -thr 0 ${PathScratch}/Lesion_${ID}_${TP[tp]}_corr.nii.gz" >> ${NameScript}

        echo "${PathSeg}/Seg_Analysis -LesWMI ${OptWMI}  -inLesCorr ${PathScratch}/Lesion_${ID}_${TP[tp]}_corr.nii.gz -inTxt2 ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${TP[tp]}_${Opt}.txt ${PathScratch}/${ModalitiesTot}_BiASM_${ID}_${TP[tp]}_${Opt}.nii.gz -mask ${PathScratch}/GIF_${ID}_${TP[tp]}_B1.nii.gz -Package 1 -SegType 1 -WeightedSeg 3 3 1 -connect -correct -inModa 2 1 3 -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${PathScratch}/ICBM_ICSF_${ID}_${TP[tp]}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${ID}_${TP[tp]}.nii.gz -inPriorsCGM ${PathScratch}/ICBM_CGM_${ID}_${TP[tp]}.nii.gz -inPriorsECSF ${PathScratch}/ICBM_ECSF_${ID}_${TP[tp]}.nii.gz -TO 1 -juxtaCorr 1 -SP ${OptSP} -LevelCorrection ${OptCL} -inArtefact ${PathScratch}/Artefacts_${ID}_${TP[tp]}.nii.gz ${ChangePathString} -ParcellationIn ${PathScratch}/GIF_Parcellation_${ID}_${TP[tp]}.nii.gz -typeVentrSeg 1 -outWM 1 -outConnect 1 -Neigh 6" >> ${NameScript}

    fi

    echo "rm ${PathScratch}/LesionWeigh* ${PathScratch}/Binary* ${PathScratch}/WMDil* ${PathScratch}/WMCard*  ${PathScratch}/LesionInit* ${PathScratch}/DataR* ${PathScratch}/DataT* ${PathScratch}/Summ*${ID}_${TP[tp]}* ${PathScratch}/LesSegHard* ${PathScratch}/Check* ${PathScratch}/BinaryNIV*" >> ${NameScript}
    echo "cp ${PathScratch}/*Co* ${PathResults}/." >> ${NameScript}
    echo "cp ${PathScratch}/*Co*.txt ${PathResults}/." >> ${NameScript}
    echo "cp ${PathScratch}/LesionMahal* ${PathResults}/.">>${NameScript}
    echo "cp ${PathScratch}/Txt* ${PathResults}/." >> ${NameScript}
    echo "cp ${PathScratch}/Out* ${PathResults}/." >> ${NameScript}
    echo "cp ${PathScratch}/Autho* ${PathResults}/." >> ${NameScript}
    echo "cp ${PathScratch}/*Infar* ${PathResults}/." >> ${NameScript}





           
      #echo "cp ${PathScratch}/*EMfO*GC* ${PathID}/${PN}_GW${OptTP}/Results/." >> ${NameScript}
      #echo "cp ${PathScratch}/*HM_GWto* ${PathID}/${PN}_GW${OptTP}/Results/." >> ${NameScript}      
done
    



    # if [ ! -f ${PathID}/${PN}_GW${OptTP}/Results/SummarisedCorr*${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}* ] & [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/${NameTogether}_BiASM_${PN}_${TP[tp]}_${SignModa}*HM${OptTP}${TAOpt}.nii.gz ]
    # then
    #     echo "${PathSeg}/Seg_Analysis -inTxt2 ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.txt ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz -mask ${PathScratch}/${PN}_GW${OptTP}_B1.nii.gz -Package 1 -SegType 1 -WeightedSeg 1 3 1 -connect -correct -inModa 2 ${ArrayModaCode[0]} ${ArrayModaCode[1]}  -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${PathScratch}/ICBM_ICSF_${PN}_GW${OptTP}.nii.gz -inPriorsECSF ${PathScratch}/ICBM_ECSF_${PN}_GW${OptTP}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${PN}_GW${OptTP}.nii.gz -inPriorsCGM ${PathScratch}/ICBM_CGM_${PN}_GW${OptTP}.nii.gz -TO 1 -juxtaCorr 1 ${ChangePathString}" >> ${NameScript}
    # fi

    # if [ ! -f ${PathID}/${PN}_GW${OptTP}/Results/IO_${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.nii.gz ]
    # then
    #     echo "${PathSeg}/Seg_Analysis -inTxt2 ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.txt ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.nii.gz -IO 1 -mask ${NameBrainmaskGW} ${ChangePathString}" >> ${NameScript}
    # fi


    NameDataCorrectedGW=${PathScratch}/DataCorrected_${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz


    echo "rm ${PathScratch}/I*HM*" >> ${NameScript}
    echo "rm ${PathScratch}/MRF*" >> ${NameScript}
    echo "rm ${PathScratch}/I*fin*" >> ${NameScript}
    echo "rm ${PathScratch}/BinaryNIV*" >> ${NameScript}
    echo "rm ${PathScratch}/PriorsAdapted*" >>${NameScript}
    echo "rm ${PathScratch}/LesSegHard*" >> ${NameScript}
    echo "rm ${PathScratch}/LeafCaract*" >> ${NameScript}
    echo "rm ${PathScratch}/WMCard*" >> ${NameScript}
    echo "rm ${PathScratch}/ST*" >> ${NameScript}
    echo "rm ${PathScratch}/BG0*GM*" >> ${NameScript}
    echo "rm ${PathScratch}/BG*" >> ${NameScript}
    echo "rm ${PathID}/*/Results/BG0*fin*" >> ${NameScript}
    echo "rm ${PathScratch}/Summarised*HM*" >> ${NameScript}
    echo "rm ${PathScratch}/T*EMfO*GC*" >> ${NameScript}
    echo "rm ${PathID}/*/Priors/GC*" >> ${NameScript}
    echo "rm ${PathScratch}/ICBM*" >> ${NameScript}
    echo "rm ${PathScratch}/*ArtConstr*" >> ${NameScript}
    echo "rm ${PathScratch}/BF*" >> ${NameScript}
    echo "rm ${PathScratch}/ConnectO*" >> ${NameScript}
    echo "rm ${PathScratch}/*Card*" >> ${NameScript}
    echo "cp ${PathScratch}/Mean* ${PathID}/${PN}_GW${OptTP}/." >> ${NameScript}
fi

echo "rm ${PathScratch}/*Card*" >> ${NameScript}
echo "rm ${PathScratch}/ICBM*" >> ${NameScript}
echo "rm ${PathScratch}/*ArtConstr*" >> ${NameScript}
echo "cp ${PathScratch}/*EMf* ${PathResults}/." >>${NameScript}
echo "cp ${PathScratch}/*BiASM*GW${OptTP}* ${PathResults}/." >>${NameScript}
echo "cp ${PathScratch}/*BiASM*HM${OptTP}* ${PathResults}/." >>${NameScript}
echo "cp ${PathScratch}/*GW${OptTP}* ${PathResults}/." >>${NameScript}
echo "cp ${PathScratch}/Mean* ${PathResults}/." >> ${NameScript}
echo "cp -r ${PathScratch}/nrr_2 ${PathID}/${PN}_GW${OptTP}/Reg/." >> ${NameScript}
echo "rm -r ${PathScratch}" >> ${NameScript}
echo "rm ${PathResults}/*ArtConstruction*" >> ${NameScript}
echo "rm ${PathResults}/BF*" >> ${NameScript}
echo "rm ${PathResults}/I*" >> ${NameScript}

function finish	{
     rm	-rf ${PathScratch}
}
trap finish EXIT ERR

if [[ $2 == "sh" ]]
then
    sh ${NameScript}
else
    qsub -S /bin/bash -l hostname="!cheech*&!chong*&!lum*&!allen*" -l h_rt=30:0:0 -l h_vmem=7.9G -l tmem=7.9G -l tscratch=20G -N LongWMH${PN} ${NameScript}
fi
