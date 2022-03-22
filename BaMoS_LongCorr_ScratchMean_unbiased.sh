#!/bin/sh
# set -x
 


if [ $# -lt 1 ]
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
OptArt=1
#############################################################################
# read the input parameters


TAString=""
OptTP=""
TAOpt=""


. $1
ScriptGW=/home/csudre/Scripts/GroupwiseReg_BaMoSLong_unbiased.sh
if ((OptTA==1))
then
TAString="-TypicalityAtlas 1"
TAOpt="TA1"
fi

echo ${OptTP}
NameScript=${PathID}/${PN}_GW${OptTP}/Scripts/GlobalISBI_${PN}${Node}.sh
echo \#\!/bin/sh > ${NameScript}
NameScriptGroupwise=${PathID}/${PN}_GW${OptTP}/Scripts/GroupwiseScript_${PN}${Node}.sh

# Be careful here with Simu : not possible to run for different Tests
PathScratch=/scratch0/${PN}_Long${OptTP}
ChangePathString="-inChangePath ${PathScratch}/"
echo "mkdir ${PathScratch}" >> ${NameScript}
echo "cp ${PathID}/${PN}_*${OptTP}/* ${PathScratch}/." >> ${NameScript}
echo "cp ${PathID}/${PN}_*${OptTP}/Results/* ${PathScratch}/." >> ${NameScript}
echo "cp ${PathID}/${PN}_*${OptTP}/Priors/* ${PathScratch}/." >> ${NameScript}
echo "cp -r ${PathID}/${PN}_*${OptTP}/Reg/nrr_2 ${PathScratch}/." >> ${NameScript}
echo "cp ${PathID}/${PN}_*0*/* ${PathScratch}/." >> ${NameScript}
echo "cp ${PathID}/${PN}_*0*/Results/* ${PathScratch}/." >> ${NameScript}
echo "cp ${PathID}/${PN}_*/Priors/* ${PathScratch}/." >> ${NameScript}
echo "cp ${PathID}/${PN}_*0*/Brainmask/* ${PathScratch}/." >> ${NameScript}
echo "cp ${PathID}/${PN}_*0*/GIF/* ${PathScratch}/." >> ${NameScript}
echo "cp ${PathID}/${PN}_*${OptTP}/GIF/* ${PathScratch}/." >> ${NameScript}

#Create folders needed for performing task, copy the appropriate data, do the registration and stacking
#mkdir ${PathID}/Study
mkdir ${PathID}/${PN}_GW${OptTP} ${PathID}/${PN}_GW${OptTP}/Priors ${PathID}/${PN}_GW${OptTP}/Brainmask ${PathID}/${PN}_GW${OptTP}/Scripts ${PathID}/${PN}_GW${OptTP}/Results
for ((tp=0;tp<${#TP[@]};tp++))
do
#Create the needed folders
    mkdir ${PathID}/${PN}_${TP[tp]} ${PathID}/${PN}_${TP[tp]}/Priors ${PathID}/${PN}_${TP[tp]}/Brainmask ${PathID}/${PN}_${TP[tp]}/Scripts ${PathID}/${PN}_${TP[tp]}/Results


    if ((DoWD==1))
    then

#echo "${PathWD}/mri_watershed -atlas ${PathID}/${PN}_${TP[tp]}/T1* ${PathID}/${PN}_${TP[tp]}/Brainmask/${PN}_${TP[tp]}_wd.nii.gz" >> ${NameScript}

        if [ ! -f ${PathID}/${PN}_${TP[tp]}/Brainmask/${PN}_${TP[tp]}_B1.nii.gz ]
        then
            echo "${PathSeg}/seg_maths ${PathID}/${PN}_01/Brainmask/${PN}_01_wd.nii.gz -lconcomp -fill -bin ${PathID}/${PN}_${TP[tp]}/Brainmask/${PN}_${TP[tp]}_lcc.nii.gz " >> ${NameScript}
            echo "${PathSeg}/seg_maths ${PathID}/${PN}_01/Brainmask/${PN}_01_lcc.nii.gz -odt char ${PathID}/${PN}_01/Brainmask/${PN}_01_B1.nii.gz " >> ${NameScript}
            echo "${PathReg}/reg_aladin -ref ${PathID}/${PN}_${TP[tp]}/T1_${PN}_${TP[tp]} -flo ${PathID}/${PN}_01/T1* -aff ${PathID}/${PN}_${TP[tp]}/AffTransT1toT1_${PN}_${TP[tp]}.txt ">> ${NameScript}
            echo "${PathReg}/reg_f3d -ref ${PathID}/${PN}_${TP[tp]}/T1* -flo ${PathID}/${PN}_01/T1* -aff ${PathID}/${PN}_${TP[tp]}/AffTransT1toT1_${PN}_${TP[tp]}.txt -cpp ${PathID}/${PN}_${TP[tp]}/CppT1toT1_${PN}_${TP[tp]}.nii.gz -maxit 100  ">> ${NameScript}
            echo "reg_resample -ref ${PathID}/${PN}_${TP[tp]}/T1* -flo ${PathID}/*01/Brainmask/*wd.nii.gz -cpp ${PathID}/${PN}_${TP[tp]}/CppT1toT1_${PN}_${TP[tp]}.nii.gz -res ${PathID}/${PN}_${TP[tp]}/Brainmask/${PN}_${TP[tp]}_trans.nii.gz " >> ${NameScript}
#echo "${PathWD}/mri_watershed -atlas ${PathID}/${PN}_${TP[tp]}/T1* ${PathID}/${PN}_${TP[tp]}/Brainmask/${PN}_${TP[tp]}_wd.nii.gz" >> ${NameScript}
            echo "${PathSeg}/seg_maths ${PathID}/${PN}_${TP[tp]}/Brainmask/${PN}_${TP[tp]}_trans.nii.gz -lconcomp -fill -bin ${PathID}/${PN}_${TP[tp]}/Brainmask/${PN}_${TP[tp]}_lcc.nii.gz " >> ${NameScript}
            echo "${PathSeg}/seg_maths ${PathID}/${PN}_${TP[tp]}/Brainmask/${PN}_${TP[tp]}_lcc.nii.gz -odt char ${PathID}/${PN}_${TP[tp]}/Brainmask/${PN}_${TP[tp]}_B1.nii.gz " >> ${NameScript}
        fi
    fi


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
        if [ ! -f ${PathID}/${PN}_${TP[tp]}/ST_${PN}_${TP[tp]}.nii.gz ]
        then
            if [ ! -f ${PathID}/${PN}_${TP[tp]}/FLAIR*${TP[tp]}.nii.gz ]
            then
                echo "${PathReg}/reg_aladin -ref ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -flo ${PathScratch}/${FinalMod[mod]}_${PN}_${TP[tp]}_init.nii.gz -res ${PathScratch}/${FinalMod[mod]}_${PN}_${TP[tp]}.nii.gz" >> ${NameScript}
            fi
        fi
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
    if [ ! -f ${PathID}/${PN}_${TP[tp]}/Priors/${PN}_${TP[tp]}_AffTransf.txt ]
    then
		echo "${PathReg}/reg_aladin -ref ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -flo ${PathICBM}/ICBM_Template.nii.gz -aff ${PathScratch}/${PN}_${TP[tp]}_AffTransf.txt" >> ${NameScript}
        echo "cp ${PathScratch}/${PN}_${TP[tp]}_AffTransf.txt ${PathID}/${PN}_${TP[tp]}/Priors/${PN}_${TP[tp]}_AffTransf.txt " >> ${NameScript}
    fi
    if [ ! -f ${PathID}/${PN}_${TP[tp]}/Priors/${PN}_${TP[tp]}_cpp.nii.gz ]
    then
		echo "${PathReg}/reg_f3d -ref ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -flo ${PathICBM}/ICBM_Template.nii.gz -aff ${PathScratch}/${PN}_${TP[tp]}_AffTransf.txt -cpp ${PathScratch}/${PN}_${TP[tp]}_cpp.nii.gz" >> ${NameScript}
        echo "cp ${PathScratch}/${PN}_${TP[tp]}_cpp.nii.gz ${PathID}/${PN}_${TP[tp]}/Priors/${PN}_${TP[tp]}_cpp.nii.gz " >> ${NameScript}
    fi
    for p in CGM DGM ECSF ICSF Out WM
    do
        if [ ! -f ${PathID}/${PN}_${TP[tp]}/Priors/ICBM_${p}_${PN}_${TP[tp]}.nii.gz ]
        then
            echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -flo ${PathICBM}/ICBM_${p}.nii.gz -cpp ${PathScratch}/${PN}_${TP[tp]}_cpp.nii.gz -res ${PathScratch}/ICBM_${p}_${PN}_${TP[tp]}.nii.gz" >> ${NameScript}
        fi
    done

# summing to obtain CSFs and GM
    if [ ! -f ${PathID}/${PN}_${TP[tp]}/Priors/*_CSFs_${PN}_${TP[tp]}.nii.gz ]
    then
		echo " ${PathSeg}/seg_maths ${PathScratch}/ICBM_ICSF_${PN}_${TP[tp]}.nii.gz -add ${PathScratch}/ICBM_ECSF_${PN}_${TP[tp]}.nii.gz ${PathScratch}/ICBM_CSFs_${PN}_${TP[tp]}.nii.gz" >> ${NameScript}
    fi
    if [ ! -f ${PathID}/${PN}_${TP[tp]}/Priors/*_GM_${PN}_${TP[tp]}.nii.gz ]
    then
        echo "${PathSeg}/seg_maths ${PathScratch}/ICBM_CGM_${PN}_${TP[tp]}.nii.gz -add ${PathScratch}/ICBM_DGM_${PN}_${TP[tp]}.nii.gz ${PathScratch}/ICBM_GM_${PN}_${TP[tp]}.nii.gz" >> ${NameScript}
    fi

    array_PriorsICBM=("${PathScratch}/ICBM_GM_${PN}_${TP[tp]}.nii.gz" "${PathScratch}/ICBM_WM_${PN}_${TP[tp]}.nii.gz" "${PathScratch}/ICBM_CSFs_${PN}_${TP[tp]}.nii.gz" "${PathScratch}/ICBM_Out_${PN}_${TP[tp]}.nii.gz")

    array_PriorsGIF=("${PathScratch}/GIF_GM_${PN}_${TP[tp]}.nii.gz" "${PathScratch}/GIF_WM_${PN}_${TP[tp]}.nii.gz" "${PathScratch}/GIF_CSF_${PN}_${TP[tp]}.nii.gz" "${PathScratch}/GIF_Out_${PN}_${TP[tp]}.nii.gz")

    if ((UseGIF==1))
    then
        PriorsArray=("Out" "CSF" "CGM" "WMI" "DGM" "Brainstem")
        for ((p=0;p<6;p++))
        do
            echo "${PathSeg}/seg_maths ${PathScratch}/*${PN}_${TP[tp]}*prior* -tp ${p} ${PathScratch}/GIF_${PriorsArray[p]}_${PN}_${TP[tp]}.nii.gz" >> ${NameScript}
        done
        echo "${PathSeg}/seg_maths ${PathScratch}/GIF_CGM_${PN}_${TP[tp]}.nii.gz -add ${PathScratch}/GIF_DGM_${PN}_${TP[tp]}.nii.gz ${PathScratch}/GIF_GM_${PN}_${TP[tp]}.nii.gz" >> ${NameScript}
        echo "${PathSeg}/seg_maths ${PathScratch}/GIF_WMI_${PN}_${TP[tp]}.nii.gz -add ${PathScratch}/GIF_Brainstem_${PN}_${TP[tp]}.nii.gz ${PathScratch}/GIF_WM_${PN}_${TP[tp]}.nii.gz" >> ${NameScript}
        if [ -f ${PathID}/${PN}_${TP[tp]}/Priors/GIF_GM_${PN}_${TP[tp]}.nii.gz ]
        then
            array_PriorsToUse=(${array_PriorsGIF[*]})
        else
            array_PriorsToUse=(${array_PriorsICBM[*]})
        fi
    else
        array_PriorsToUse=(${array_PriorsICBM[*]})
    fi

    ProgMod=""
    if ((${#array_ToStack[@]}==2))
    then
        ProgMod="-progMod 1 2 "
    else
        ProgMod="-progMod 2 2 1"
    fi
    # Perform Seg_BiASM EMfO on each
    if [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/${NameTogether}_EMfO_${PN}_${TP[tp]}*.nii.gz ]
    then
        echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_${TP[tp]}*TIV.nii.gz -bin -odt char ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz" >> ${NameScript}
        echo "cp ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz ${PathScratch}/${PN}_${TP[tp]}_GIF_B1.nii.gz" >> ${NameScript}
		echo "${PathSeg}/Seg_BiASM -in ${#array_ToStack[@]} ${array_ToStack[*]} -priors 4 ${array_PriorsToUse[*]} -mask ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz -out 2 ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}.nii.gz ${PathScratch}/${NameTogether}_EMfOf_${PN}_${TP[tp]}.nii.gz -txt_out ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}.txt -bc_order 3 -CovPriors 8 -BFP 1 -maxRunEM 1 ${ProgMod} -AtlasSmoothing 1 1 -AtlasWeight 1 ${AW}  -SMOrder 0 -KernelSize 3 -PriorsKept 5 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -deleteUW 1 -outliersM 3 -outliersW 0.01 -init_splitUnif 0 -splitAccept 0 -unifTot 1 -MRF 1 -GMRF ${NameGMatrix} -juxtaCorr 1 -meanPriors 0 -priorDGM ${PathScratch}/ICBM_DGM_${PN}_${TP[tp]}.nii.gz -BiASM 0 ${TAString}" >> ${NameScript}
        echo "rm ${PathScratch}/BG* " >> ${NameScript}
        echo "rm ${PathScratch}/MRF* " >> ${NameScript}

	fi
# Seg_Analysis for IO

    if [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/I_${NameTogether}_EMfO_${PN}_${TP[tp]}${SignModa}.nii.gz ]
    then
		echo "${PathSeg}/Seg_Analysis -inTxt2 ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}${SignModa}.txt ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}${SignModa}.nii.gz -IO 1 -mask ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz ${ChangePathString}" >> ${NameScript}
    fi
	echo "cp ${PathScratch}/*${NameTogether}*EMfO*${PN}_${TP[tp]}* ${PathID}/${PN}_${TP[tp]}/Results/." >> ${NameScript}
done

if ((DoGW==1))
then

    if  [ ! -f ${PathID}/*GW${OptTP}/T1*GW* ] && [ ! -f ${PathID}/*GW${OptTP}/Results/T1_*GW* ]
    then
# Going to the GW part
#Create the appropriate folders 
        mkdir ${PathID}/${PN}_GW${OptTP} ${PathID}/${PN}_GW${OptTP}/Results/ ${PathID}/${PN}_GW${OptTP}/Brainmask ${PathID}/${PN}_GW${OptTP}/Scripts ${PathID}/${PN}_GW${OptTP}/Priors ${PathID}/${PN}_GW${OptTP}/Reg

# Create the file for the group wise
        NameFileGW=${PathID}/${PN}_GW${OptTP}/Reg/GroupwiseParam_${PN}.sh
        echo "" > ${NameFileGW}
        echo "#!/bin/sh" >> ${PathID}/${PN}_GW${OptTP}/Reg/GroupwiseParam_${PN}.sh
        echo "export PathSeg=$PathSeg" >> ${PathID}/${PN}_GW${OptTP}/Reg/GroupwiseParam_${PN}.sh
        echo "export PathReg=$PathReg" >> ${PathID}/${PN}_GW${OptTP}/Reg/GroupwiseParam_${PN}.sh
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

#Once done, get the priors on the new mean model 

#Destack the result
        array_DCGW=()
        NameFinalAverage=${PathScratch}/nrr_2/average_nonrigid_it_2.nii.gz
        for ((i=0;i<${#array_ToStack[@]};i++))
        do
            echo "${PathSeg}/seg_maths ${NameFinalAverage} -tp ${i} ${PathScratch}/${FinalMod[i]}_${PN}_GW${OptTP}.nii.gz" >> ${NameScript}
            array_DCGW=(${array_DCGW[*]} "${PathScratch}/${FinalMod[i]}_${PN}_GW${OptTP}.nii.gz" )
        done

        echo "cp -r ${PathScratch}/nrr_2 ${PathID}/${PN}_GW${Opt}/Reg/." >> ${NameScript}

# First do a reg_resample on all time points and all atlases using the final cpp transformation
        echo "${#TP[@]} time points "

        for p in CGM DGM GM WM Out ICSF ECSF CSFs
        do
            if [ ! -f ${PathScratch}/ICBM_${p}_${PN}_GW${OptTP}.nii.gz ]
            then
            array_priorspre=()
                for ((tp=0;tp<${#TP[@]};tp++))
                do
                    if [ ! -f ${PathScratch}/ICBM_${p}_${PN}_${TP[tp]}_preGW.nii.gz ]
                    then
                        echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_GW${OptTP}.nii.gz -flo  ${PathScratch}/ICBM*_${p}_${PN}_${TP[tp]}* -cpp  ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}*it2.nii.gz -res ${PathScratch}/ICBM_${p}_${PN}_${TP[tp]}_preGW.nii.gz" >> ${NameScript}
                    fi
                    array_priorspre=(${array_priorspre[*]} "${PathScratch}/ICBM_${p}_${PN}_${TP[tp]}_preGW.nii.gz")
                done

        # Mean over the transformed time points

                echo "${PathSeg}/Seg_Analysis -meanImages ${#TP[@]} ${array_priorspre[*]} -outMean ${PathScratch}/ICBM_${p}_${PN}_GW${OptTP}.nii.gz" >> ${NameScript}
            fi
        done
    fi


    for p in CGM DGM GM WM Out ICSF ECSF CSFs
    do
        if [ ! -f ${PathID}/${PN}_GW${OptTP}/Priors/ICBM_${p}_${PN}_GW${OptTP}.nii.gz ]
        then
            array_priorspre=()
            for ((tp=0;tp<${#TP[@]};tp++))
            do
                if [ ! -f {PathID}/${PN}_GW${OptTP}/Priors/ICBM_${p}_${PN}_${TP[tp]}_preGW.nii.gz ]
                then
                    echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_GW${OptTP}.nii.gz -flo  ${PathScratch}/ICBM*_${p}_${PN}_${TP[tp]}* -cpp  ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}*it2.nii.gz -res ${PathScratch}/ICBM_${p}_${PN}_${TP[tp]}_preGW.nii.gz" >> ${NameScript}
                fi
                array_priorspre=(${array_priorspre[*]} "${PathScratch}/ICBM_${p}_${PN}_${TP[tp]}_preGW.nii.gz")
            done

    # Mean over the transformed time points
            if [ ! -f ${PathID}/${PN}_GW${OptTP}/Priors/ICBM_${p}_${PN}_GW${OptTP}.nii.gz ]
            then
                echo "${PathSeg}/Seg_Analysis -meanImages ${#TP[@]} ${array_priorspre[*]} -outMean ${PathScratch}/ICBM_${p}_${PN}_GW${OptTP}.nii.gz" >> ${NameScript}
                echo "rm ${PathID}/${PN}_GW${OptTP}/Priors/*preGW*" >> ${NameScript}
            fi
        fi
    done


    if [ ! -f ${PathID}/${PN}_GW${OptTP}/MeanT1_${PN}.nii.gz ]
    then
        array_meanpre=();
        for ((tp=0;tp<${#TP[@]};tp++))
        do
            if [ ! -f ${PathID}/${PN}_GW${OptTP}/T1_${PN}_${TP[tp]}_preGW.nii.gz ]
            then
                echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_GW${OptTP}.nii.gz -flo  ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -cpp  ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}*it2.nii.gz -res ${PathScratch}/T1_${PN}_${TP[tp]}_preGW.nii.gz" >> ${NameScript}
            fi
            array_meanpre=(${array_meanpre[*]} "${PathScratch}/T1_${PN}_${TP[tp]}_preGW.nii.gz")
        done

        if [ -f ${PathID}/${PN}_GW${OptTP}/MeanT1_${PN}.nii.gz ]
        then
            echo "rm ${PathScratch}/*/*log*" >> ${NameScript}
            echo  "rm ${PathScratch}/*/*inliers* " >> ${NameScript}
            echo "rm ${PathScratch}/*/*res* " >> ${NameScript}
            echo "rm ${PathScratch}/*/mask* " >> ${NameScript}
        fi

# Mean over the transformed time points

        echo "${PathSeg}/Seg_Analysis -meanImages ${#TP[@]} ${array_meanpre[*]} -outMean ${PathScratch}/MeanT1_${PN}.nii.gz" >> ${NameScript}
        echo "cp ${PathScratch}/Mean* ${PathID}/${PN}_GW${OptTP}/." >> ${NameScript}

        echo "rm ${PathID}/${PN}_GW${OptTP}/Reg/*/*log*" >> ${NameScript}
        echo  "rm ${PathID}/${PN}_GW${OptTP}/Reg/*/*inliers* " >> ${NameScript}
        echo "rm ${PathID}/${PN}_GW${OptTP}/Reg/*/*res* " >> ${NameScript}
        echo "rm ${PathID}/${PN}_GW${OptTP}/Reg/*/mask* " >> ${NameScript}
        echo "rm ${PathID}/${PN}_GW${OptTP}/Reg/*/*matched* " >> ${NameScript}
        echo "rm -r ${PathID}/${PN}_GW${OptTP}/Reg/aff* " >> ${NameScript}
        echo "rm -r ${PathID}/${PN}_GW${OptTP}/Reg/nrr_1 " >> ${NameScript}
        echo "rm -r ${PathID}/${PN}_GW${OptTP}/*preGW* " >> ${NameScript}
        echo "rm -r ${PathScratch}/*preGW* " >> ${NameScript}
    fi



# if ((UseGIF==1))
#then
# for p in GM WM CSF Out
# do
# if [ ! -f ${PathScratch}/GIF_${p}_${PN}_GW${OptTP}.nii.gz ]
#then
# array_priorspre=()
# for ((tp=0;tp<${#TP[@]};tp++))
#  do
# if [ ! -f ${PathID}/${PN}_GW${OptTP}/Priors/GIF_${p}_${PN}_${TP[tp]}_preGW.nii.gz ]
# then
#  echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_GW${OptTP}.nii.gz -flo  ${PathScratch}/GIF*_${p}_*${PN}_${TP[tp]}* -cpp  ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}*it2.nii.gz -res ${PathScratch}/GIF_${p}_${PN}_${TP[tp]}_preGW.nii.gz" >> ${NameScript}
# fi
#array_priorspre=(${array_priorspre[*]} "${PathScratch}/GIF_${p}_${PN}_${TP[tp]}_preGW.nii.gz")
# done

# Mean over the transformed time points

# echo "${PathSeg}/Seg_Analysis -meanImages ${#TP[@]} ${array_priorspre[*]} -outMean ${PathScratch}/GIF_${p}_${PN}_GW${OptTP}_bt.nii.gz" >> ${NameScript}
#  echo "${PathSeg}/seg_maths ${PathScratch}/GIF_${p}_${PN}_GW${OptTP}_bt.nii.gz -thr 0.05  ${PathScratch}/GIF_${p}_${PN}_GW${OptTP}.nii.gz" >> ${NameScript}
#fi
#done
#fi

# Similar but on brain mask

#array_BMpre=()
#if [ ! -f ${PathID}/${PN}_GW${OptTP}/Brainmask/${PN}_GW${OptTP}_B1.nii.gz ]
#then
#  for ((tp=0;tp<${#TP[@]};tp++))
# do
#   if [ ! -f ${PathID}/${PN}_GW${OptTP}/Brainmask/${PN}_${TP[tp]}_preB1.nii.gz ]
#   then
#    echo "${PathReg}/reg_resample -ref ${PathID}/${PN}_GW${OptTP}/T1_${PN}_GW${OptTP}.nii.gz -flo  ${PathID}/${PN}_${TP[tp]}/Brainmask/*B1* -cpp  ${PathID}/${PN}_GW${OptTP}/Reg/nrr_2/*cpp*${PN}_${TP[tp]}*it2.nii.gz -res ${PathID}/${PN}_GW${OptTP}/Brainmask/${PN}_${TP[tp]}_preB1.nii.gz" >> ${NameScript}
#   echo "${PathSeg}/seg_maths ${PathID}/${PN}_GW${OptTP}/Brainmask/${PN}_${TP[tp]}_preB1.nii.gz -odt float ${PathID}/${PN}_GW${OptTP}/Brainmask/${PN}_${TP[tp]}_preB1.nii.gz " >> ${NameScript}
# fi
# array_BMpre=(${array_BMpre[*]} "${PathID}/${PN}_GW${OptTP}/Brainmask/${PN}_${TP[tp]}_preB1.nii.gz")
# done

# Mean over the transformed time points

#echo "${PathSeg}/Seg_Analysis -meanImages ${#TP[@]} ${array_BMpre[*]} -outMean ${PathID}/${PN}_GW${OptTP}/Brainmask/${PN}_GW${OptTP}_B1.nii.gz" >> ${NameScript}
#echo "${PathSeg}/seg_maths ${PathID}/${PN}_GW${OptTP}/Brainmask/${PN}_GW${OptTP}_B1.nii.gz -thr 0.5 -bin -dil 1 -fill -ero 1 -odt char ${PathID}/${PN}_GW${OptTP}/Brainmask/${PN}_GW${OptTP}_B1.nii.gz" >> ${NameScript}
#fi

# In case there is a false positive correction using the parcellation for the cingulate gyrus and the subcallosal area




fi

UseGIF=1
if ((UseGIF==1))
then
#if [ ! -f ${PathID}/PN_GW/Priors/${PN}_GW${OptTP}_Parcellation.nii.gz ]
#  then
    echo "cp ${PathScratch}/Mean*Parcellation* ${PathScratch}/${PN}_GW${OptTP}_Parcellation.nii.gz ">>${NameScript}
# fi
#101 102 139 140 187 188 169 170
    if [ ! -f ${PathID}/${PN}_GW${OptTP}/${PN}_GW${OptTP}_Artefacts.nii.gz ]
    then
        if ((${#ArtefactArray[@]}>0))
        then
            stringAddition="${PathScratch}/${PN}_GW${OptTP}_ArtConstruction_${i}.nii.gz "
            for ((i=0;i<${#ArtefactArray[@]};i++))
            do
                Value=${ArtefactArray[i]}
                ValueMin=`echo "$Value - 0.5" | bc -l`
                ValueMax=`echo "$Value + 0.5" | bc -l`
                echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_GW${OptTP}_Parcellation.nii.gz -thr $ValueMin -uthr $ValueMax -bin ${PathScratch}/${PN}_GW${OptTP}_ArtConstruction_${i}.nii.gz " >> ${NameScript}
                stringAddition="${stringAddition} -add ${PathScratch}/${PN}_GW${OptTP}_ArtConstruction_${i}.nii.gz "
            done
            echo "${PathSeg}/seg_maths ${stringAddition} -bin ${PathScratch}/${PN}_GW${OptTP}_Artefacts.nii.gz " >> ${NameScript}
            echo "cp ${PathScratch}/${PN}_GW${OptTP}_Artefacts.nii.gz ${PathID}/${PN}_GW${OptTP}/Priors/." >> ${NameScript}
        fi
    fi
    for ((tp=0;tp<${#TP[@]};tp++))
    do
        if [ ! -f ${PathID}/${PN}_GW${OptTP}/${PN}_${TP[tp]}_Artefacts.nii.gz ]
        then
            echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -flo  ${PathScratch}/${PN}_GW${OptTP}_Artefacts.nii.gz -res  ${PathScratch}/${PN}_${TP[tp]}_Artefacts.nii.gz -cpp ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}*it2*backward*.nii.gz " >> ${NameScript}
        fi
    done


    array_artefacts=()
    for ((tp=0;tp<${#TP[@]};tp++))
    do
        if [ ! -f ${PathScratch}/ICBM_${p}_${PN}_${TP[tp]}_preGW.nii.gz ]
        then
            echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_GW${OptTP}.nii.gz -flo  ${PathScratch}/${PN}_${TP[tp]}_Artefacts.nii.gz -cpp  ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}*it2.nii.gz -res ${PathScratch}/${PN}_${TP[tp]}_preGW_Artefacts.nii.gz" >> ${NameScript}
        fi
        array_artefacts=(${array_artefacts[*]} "${PathScratch}/${PN}_${TP[tp]}_preGW_Artefacts.nii.gz")
    done



    if [ ! -f ${PathID}/${PN}_GW${OptTP}/Priors/${PN}_GW${OptTP}_Artefacts.nii.gz ]
    then
        echo "${PathSeg}/Seg_Analysis -meanImages ${#TP[@]} ${array_artefacts[*]} -outMean ${PathScratch}/${PN}_GW${OptTP}_Artefacts.nii.gz" >> ${NameScript}
        echo "rm ${PathID}/${PN}_GW${OptTP}/Priors/*preGW*" >> ${NameScript}
    fi


    for p in CGM DGM GM WM Out CSF
    do
        if [ ! -f ${PathScratch}/GIF_${p}_${PN}_GW${OptTP}.nii.gz ]
        then
            array_priorspre=()
            for ((tp=0;tp<${#TP[@]};tp++))
            do
                if [ ! -f ${PathScratch}/ICBM_${p}_${PN}_${TP[tp]}_preGW.nii.gz ]
                then
                    echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_GW${OptTP}.nii.gz -flo  ${PathScratch}/GIF*_${p}_${PN}_${TP[tp]}.nii.gz -cpp  ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}*it2.nii.gz -res ${PathScratch}/GIF_${p}_${PN}_${TP[tp]}_preGW.nii.gz" >> ${NameScript}
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
        if [ ! -f ${PathScratch}/GIF_B1_${PN}_${TP[tp]}_preGW.nii.gz ]
        then
            echo "cp ${PathScratch}/${PN}_${TP[tp]}_GIF_TIV.nii.gz ${PathScratch}/${PN}_${TP[tp]}_GIF_B1.nii.gz" >> ${NameScript}
            echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_GW${OptTP}.nii.gz -flo  ${PathScratch}/${PN}_${TP[tp]}_GIF_B1.nii.gz -cpp  ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}*it2.nii.gz -res ${PathScratch}/GIF_B1_${PN}_${TP[tp]}_preGW.nii.gz" >> ${NameScript}
        fi
        array_brainmask=(${array_brainmask[*]} "${PathScratch}/GIF_B1_${PN}_${TP[tp]}_preGW.nii.gz")
    done

# Mean over the transformed time points

    echo "${PathSeg}/Seg_Analysis -meanImages ${#TP[@]} ${array_brainmask[*]} -outMean ${PathScratch}/GIF_B1_${PN}_GW${OptTP}.nii.gz" >> ${NameScript}
    echo "${PathSeg}/seg_maths ${PathScratch}/GIF_B1_${PN}_GW${OptTP}.nii.gz -thr 0.5 -bin -odt char ${PathScratch}/${PN}_GW${OptTP}_B1.nii.gz" >> ${NameScript}
    echo "cp ${PathScratch}/GIF_B1_${PN}_GW${OptTP}.nii.gz ${PathID}/${PN}_GW${OptTP}/Results/${PN}_GW${OptTP}_GIF_B1.nii.gz" >> ${NameScript}
    echo "cp ${PathScratch}/GIF_B1_${PN}_GW${OptTP}.nii.gz ${PathScratch}/${PN}_GW${OptTP}_GIF_B1.nii.gz" >> ${NameScript}
    array_brainmask=()






# First do a reg_resample on all time points and all atlases using the final cpp transformation
    echo "${#TP[@]} time points "

    for p in CGM DGM GM WM Out ICSF ECSF CSFs
    do
        if [ ! -f ${PathScratch}/ICBM_${p}_${PN}_GW${OptTP}.nii.gz ]
        then
            array_priorspre=()
            for ((tp=0;tp<${#TP[@]};tp++))
            do
                if [ ! -f ${PathScratch}/ICBM_${p}_${PN}_${TP[tp]}_preGW.nii.gz ]
                then
                    echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_GW${OptTP}.nii.gz -flo  ${PathScratch}/ICBM*_${p}_${PN}_${TP[tp]}.nii.gz -cpp  ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}*it2.nii.gz -res ${PathScratch}/ICBM_${p}_${PN}_${TP[tp]}_preGW.nii.gz" >> ${NameScript}
                fi
                array_priorspre=(${array_priorspre[*]} "${PathScratch}/ICBM_${p}_${PN}_${TP[tp]}_preGW.nii.gz")
            done

        # Mean over the transformed time points

            echo "${PathSeg}/Seg_Analysis -meanImages ${#TP[@]} ${array_priorspre[*]} -outMean ${PathScratch}/ICBM_${p}_${PN}_GW${OptTP}.nii.gz" >> ${NameScript}
        fi
    done


    if [ ! -f ${PathID}/${PN}_GW${OptTP}/MeanT1_${PN}.nii.gz ]
    then
        array_meanpre=();
        for ((tp=0;tp<${#TP[@]};tp++))
        do
            if [ ! -f ${PathID}/${PN}_GW${OptTP}/T1_${PN}_${TP[tp]}_preGW.nii.gz ]
            then
                echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_GW${OptTP}.nii.gz -flo  ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -cpp  ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}*it2.nii.gz -res ${PathScratch}/T1_${PN}_${TP[tp]}_preGW.nii.gz" >> ${NameScript}
            fi
            array_meanpre=(${array_meanpre[*]} "${PathScratch}/T1_${PN}_${TP[tp]}_preGW.nii.gz")
        done

        if [ -f ${PathID}/${PN}_GW${OptTP}/MeanT1_${PN}.nii.gz ]
        then
            echo "rm ${PathScratch}/*/*log*" >> ${NameScript}
            echo  "rm ${PathScratch}/*/*inliers* " >> ${NameScript}
            echo "rm ${PathScratch}/*/*res* " >> ${NameScript}
            echo "rm ${PathScratch}/*/mask* " >> ${NameScript}

        fi
    fi
fi

if ((DoMS==1))
then
# Model selection
    echo "We have to do MS seg"

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
    if ((UseGIF==1))
    then
        if [ -f ${PathID}/${PN}_GW${OptTP}/GIF/Mean*prior* ]
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
        elif [ -f ${PathID}/${PN}*GW${OptTP}/Priors/GIF_GM_${PN}_${TP[0]}.nii.gz ]
        then
            array_PriorsToUse=(${array_PriorsGIF[*]})
            NameBrainmaskGW=${PathScratch}/${PN}_GW${OptTP}_GIF_B1.nii.gz
        elif [ -f ${PathID}/${PN}*GW${OptTP}/Priors/GIF_GM_${PN}_GW${OptTP}.nii.gz ] || [ -f ${PathID}/${PN}*${TP[0]}/GIF/*prior* ]
        then
            array_PriorsToUse=(${array_PriorsGIF[*]})
            NameBrainmaskGW=${PathScratch}/${PN}_GW${OptTP}_GIF_B1.nii.gz
        else
            array_PriorsToUse=(${array_PriorsICBM[*]})
            NameBrainmaskGW=${PathScratch}/${PN}_GW${OptTP}_B1.nii.gz
        fi
    else
        array_PriorsToUse=(${array_PriorsICBM[*]})
    fi



    if [ ! -f ${PathID}/${PN}_GW${OptTP}/Results/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.txt ]
    then
        if ((Do2==0))
        then
            echo "Do2=0"
            echo "${PathSeg}/Seg_BiASM -in ${#array_DCGW[@]} ${array_DCGW[*]} -priors 4 ${array_PriorsToUse[*]} -mask ${NameBrainmaskGW} -out 2 ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz ${PathScratch}/${NameTogether}_BiASMG_${PN}_GW${OptTP}.nii.gz -txt_out ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.txt -bc_order 0 -inDC_flag 1 -CovPriors 8 -BFP 1 -maxRunEM ${MRE} ${ProgMod} -AtlasSmoothing 1 1 -AtlasWeight 1 ${AW}  -SMOrder 0 -KernelSize 3 -PriorsKept 5 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -deleteUW 1 -outliersM 3 -outliersW ${OW} -init_splitUnif 0 -splitAccept 0 -unifTot 1 -MRF 1 -GMRF ${NameGMatrix} -juxtaCorr 1 -priorDGM ${PathScratch}/ICBM_DGM_${PN}_GW${OptTP}.nii.gz ${TAString}" >> ${NameScript}

        else
            if [ ! -f ${PathID}/${PN}_GW${OptTP}/Results/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.txt ]
            then
                echo "Segmentation 2 to do"
                echo "${PathSeg}/Seg_BiASM -in 2 ${array_DCGW[0]} ${array_DCGW[1]} -priors 4 ${array_PriorsToUse[*]} -mask ${NameBrainmaskGW} -out 2 ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz ${PathScratch}/${NameTogether}_BiASMG_${PN}_GW${OptTP}.nii.gz -txt_out ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.txt -bc_order 0 -inDC_flag 1 -CovPriors 8 -BFP 1 -maxRunEM ${MRE} -AtlasSmoothing 1 1 -AtlasWeight 1 ${AW}  -SMOrder 0 -KernelSize 3 -PriorsKept 5 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -deleteUW 1 -outliersM 3 -outliersW ${OW} -init_splitUnif 0 -splitAccept 0 -unifTot 1 -MRF 1 -GMRF ${NameGMatrix} -juxtaCorr 1 -progMod 0 -priorDGM ${PathScratch}/ICBM_DGM_${PN}_GW${OptTP}.nii.gz ${TAString} " >> ${NameScript}
            fi

        fi

        if ((Do2==1))
        then
            echo "2 moda to consider first"
            if [ ! -f ${PathID}/${PN}_GW${OptTP}/Results/SummarisedCorr*${NameTogether}_BiASM_${PN}_GW${OptTP}* ]
            then
                echo "${PathSeg}/Seg_Analysis -inTxt2 ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.txt ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz -mask ${NameBrainmaskGW} -Package 1 -SegType 1 -WeightedSeg 1 3 1 -connect -correct -inModa 2 ${ArrayModaCode[0]} ${ArrayModaCode[1]}  -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${PathScratch}/ICBM_ICSF_${PN}_GW${OptTP}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${PN}_GW${OptTP}.nii.gz -TO 1 -juxtaCorr 1 ${ChangePathString}" >> ${NameScript}
            fi
# EMf using summarised Corrected for atlases at each time point
#read -a NameSummarisedCorr <<< $(ls ${PathID}/${PN}_GW${OptTP}/Results/SummarisedCorr*${NameTogether}_BiASM_${PN}_GW${OptTP}.nii.gz)
NameSummarisedCorr=${PathScratch}/SummarisedCorrected_WS3WT1WC1JC1ST1CL2CIV1_TO22${NameTogether}_BiASM_${PN}_GW${OptTP}.nii.gz
        fi
    fi
#fi

# Obtention of final segmentation

    if ((Do3==1))
    then
        if [ ! -f ${PathID}/${PN}_GW${OptTP}/Results/LesionCorr*TO22${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}*${SignModa}.nii.gz ]
        then

            if [ ! -f ${PathID}/${PN}_GW${OptTP}/Results/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}*${SignModa}.nii.gz ]
            then
                if((Do3==1))
                then
                    echo "${PathSeg}/Seg_Analysis -inTxt2 ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.txt ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz -mask ${NameBrainmaskGW} -Package 1 -SegType 1 -WeightedSeg 1 3 1 -connect -correct -inModa 2  ${ArrayModaCode[0]} ${ArrayModaCode[1]}  -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${PathScratch}/ICBM_ICSF_${PN}_GW${OptTP}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${PN}_GW${OptTP}.nii.gz -juxtaCorr 1 -TO 1 ${ChangePathString}" >> ${NameScript}


                    echo "${PathSeg}/Seg_BiASM -SMOrder 0 -averageGC 0.5 ${PathScratch}/SummarisedCorrected_WS3WT1WC1JC1ST1CL2CIV1_TO22${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz  ${PathScratch}/SummarisedCorrected_WS3WT1WC1JC1ST1CL2CIV1_TO22${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz -in 3 ${PathScratch}/T1_${PN}_GW${OptTP}.nii.gz  ${PathScratch}/FLAIR_${PN}_GW${OptTP}.nii.gz  ${PathScratch}/T2_${PN}_GW${OptTP}.nii.gz  -averageIO 0.5 ${PathScratch}/IO_${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz ${PathScratch}/IO_${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz -inDC ${PathScratch}/DataCorrected_${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz -txt_in ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.txt -meanPriors 0 -txt_out ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.txt -out 2 ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.nii.gz ${PathScratch}/${NameTogether}_BiASM_${PN}${TAOpt}${SignModa}.nii.gz -meanPriors 0 -juxtaCorr 0 -mask ${PathScratch}/${PN}_GW${OptTP}_B1.nii.gz -MRF 1 -GMRF ${NameGMatrix} -bc_order 0 -CovPriors 8 -priorDGM ${PathScratch}/ICBM_DGM_${PN}_GW${OptTP}.nii.gz -BiASM 0 -maxRunEM 50 -outliersM 3 -BiASM 1 -PriorsKept 5 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4" >> ${NameScript}
                fi
            fi

            echo "${PathSeg}/Seg_Analysis -inTxt2 ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.txt ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.nii.gz -mask ${PathScratch}/${PN}_GW${OptTP}_B1.nii.gz -Package 1 -SegType 1 -WeightedSeg 1 3 1 -connect -correct -inModa ${#ArrayModaCode[@]} ${ArrayModaCode[*]} -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${PathScratch}/ICBM_ICSF_${PN}_GW${OptTP}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${PN}_GW${OptTP}.nii.gz -TO 1 -juxtaCorr 1 ${ChangePathString}" >> ${NameScript}

# EMf using summarised Corrected for atlases at each time point
#read -a NameSummarisedCorr <<< $(ls ${PathID}/${PN}_GW${OptTP}/Results/SummarisedCorr*${NameTogether}_BiASM_${PN}_GW${OptTP}${SignModa}.nii.gz)
            NameSummarisedCorr=${PathScratch}/SummarisedCorrected_WS3WT1WC1JC1ST1CL2CIV1_TO22${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.nii.gz
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
                echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -flo ${NameSummarisedCorr} -cpp ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}*it2*backward*.nii.gz -res  ${PathScratch}/GC_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz " >> ${NameScript}
                if ((UseGIF==1))
                then
                    echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -flo ${PathScratch}/${PN}_GW${OptTP}*B1.nii.gz -cpp ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}*it2*backward*.nii.gz -res  ${PathScratch}/${PN}_${TP[tp]}_GIF.nii.gz " >> ${NameScript}
                    echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_${TP[tp]}_GIF.nii.gz -bin -fill -odt char ${PathScratch}/${PN}_${TP[tp]}_GIF_B1.nii.gz " >> ${NameScript}
                    array_Brainmask=("${array_Brainmask[*]}" ${PathScratch}/${PN}_${TP[tp]}_GIF_B1.nii.gz )
                else
                    array_Brainmask=("${array_Brainmask[*]}" ${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz )
                fi


# Destack it into priors
                array_Priors=()
                for i in 0 1 2 3
                do
                    echo "${PathSeg}/seg_maths ${PathScratch}/GC_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz -tp ${i} ${PathScratch}/GC_${FinalPriors[i]}_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz " >> ${NameScript}
                    array_Priors=(${array_Priors[*]} "${PathScratch}/GC_${FinalPriors[i]}_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz ")
                    echo ${PathScratch}
                done

# EMf using these atlases with progressive and no BF
                if [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.nii.gz ] & [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/${NameTogether}_BiASM_${PN}_${TP[tp]}_${TAOpt}${SignModa}*HM${OptTP}${TAOpt}.nii.gz ]
				then
                    echo "${PathSeg}/Seg_BiASM -in ${#array_DC[@]} ${array_DC[*]} -priors 4 ${array_Priors[*]} -mask ${array_Brainmask[TP[tp]]} -out 2 ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}.nii.gz ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}.nii.gz -txt_out ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}.txt -bc_order 0 -CovPriors 8 -BFP 1 -maxRunEM 1 ${ProgMod} -AtlasSmoothing 1 1 -AtlasWeight 1 ${AW}  -SMOrder 0 -KernelSize 3 -PriorsKept 5 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -deleteUW 1 -outliersM 3 -outliersW 0.01 -init_splitUnif 0 -splitAccept 0 -unifTot 1 -MRF 1 -GMRF ${NameGMatrix} -juxtaCorr 1 -BiASM 0 -inDC_flag 1 -priorDGM ${PathScratch}/ICBM_DGM_${PN}_${TP[tp]}.nii.gz -BiASM 0 -maxRunEM 1 -juxtaCorr 1 ${TAString}" >> ${NameScript}
                    echo "rm ${PathScratch}/BG*" >> ${NameScript}
                    echo "rm ${PathScratch}/GC_*" >> ${NameScript}
                fi

# Get IO atlases from that one
                if [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/IO_${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${SignModa}.nii.gz ] & [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/${NameTogether}_BiASM_${PN}_${TP[tp]}_${TAOpt}${SignModa}*HM${OptTP}${TAOpt}.nii.gz ]
                then
                    echo "${PathSeg}/Seg_Analysis -inTxt2 ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.txt ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.nii.gz -IO 1 -mask ${array_Brainmask[TP[tp]]} ${ChangePathString}" >> ${NameScript}
                fi

            done
        fi
    fi



    if [ ! -f ${PathID}/${PN}_GW${OptTP}/Results/SummarisedCorr*${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}* ] & [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/${NameTogether}_BiASM_${PN}_${TP[tp]}_${SignModa}*HM${OptTP}${TAOpt}.nii.gz ]
    then
        echo "${PathSeg}/Seg_Analysis -inTxt2 ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.txt ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz -mask ${PathScratch}/${PN}_GW${OptTP}_B1.nii.gz -Package 1 -SegType 1 -WeightedSeg 1 3 1 -connect -correct -inModa 2 ${ArrayModaCode[0]} ${ArrayModaCode[1]}  -inRuleTxt ${RuleFileName} -WMCard 1 -inPriorsICSF ${PathScratch}/ICBM_ICSF_${PN}_GW${OptTP}.nii.gz -inPriorsECSF ${PathScratch}/ICBM_ECSF_${PN}_GW${OptTP}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${PN}_GW${OptTP}.nii.gz -inPriorsCGM ${PathScratch}/ICBM_CGM_${PN}_GW${OptTP}.nii.gz -TO 1 -juxtaCorr 1 ${ChangePathString}" >> ${NameScript}
    fi


    NameDataCorrectedGW=${PathScratch}/DataCorrected_${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}.nii.gz
    for ((tp=0;tp<${#TP[@]};tp++))
    do
        NameBrainmask=""
        if ((UseGIF==1))
        then
            echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -flo ${PathScratch}/${PN}_GW${OptTP}_GIF_B1.nii.gz -cpp ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}*it2*backward*.nii.gz -res  ${PathScratch}/${PN}_${TP[tp]}_GIF.nii.gz " >> ${NameScript}
            echo "${PathSeg}/seg_maths ${PathScratch}/${PN}_${TP[tp]}_GIF.nii.gz -thr 0.5 -bin -fill -odt char ${PathScratch}/${PN}_${TP[tp]}_GIF_B1.nii.gz " >> ${NameScript}
            NameBrainmask=${PathScratch}/${PN}_${TP[tp]}_GIF_B1.nii.gz
        else
            NameBrainmask=${PathScratch}/${PN}_${TP[tp]}_B1.nii.gz
        fi


        NameDC=${PathScratch}/DataCorrected_${NameTogether}_EMfO_${PN}_${TP[tp]}.nii.gz
        array_DC=()
        for ((i=0;i<${#array_ToStack[@]};i++))
        do
            echo "${PathSeg}/seg_maths ${NameDC} -tp ${i} ${PathScratch}/${FinalMod[i]}_${PN}_${TP[tp]}_DC.nii.gz" >> ${NameScript}
            array_DC=(${array_DC[*]} "${PathScratch}/${FinalMod[i]}_${PN}_${TP[tp]}_DC.nii.gz" )
        done

        NameSummarisedCorr=${PathScratch}/SummarisedCorrected_WS3WT1WC1JC1SP1ST1CL2CIV1_TO22${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.nii.gz

# resample SummarisedCorr into initial space
        echo "${PathReg}/reg_resample -ref ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -flo ${NameSummarisedCorr} -cpp ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}*it2*backward*.nii.gz -res  ${PathScratch}/GC_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz " >> ${NameScript}


# Destack it into priors
        array_Priors=()
        for i in 0 1 2 3
        do
            echo "${PathSeg}/seg_maths ${PathScratch}/GC_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz -tp ${i} ${PathScratch}/GC_${FinalPriors[i]}_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz " >> ${NameScript}
            array_Priors=(${array_Priors[*]} "${PathScratch}/GC_${FinalPriors[i]}_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz ")
        done

# EMf using these atlases with progressive and no BF
        if [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.nii.gz ] || [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.txt ] || [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/DataCorrected_HM_GWto${TP[tp]}_${NameTogether}_${PN}_0.nii.gz ]
        then
            if [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.nii.gz ]
            then
                echo "Have to redo EMfO absolutely !!!"
                echo "${PathSeg}/Seg_BiASM -in ${#array_DC[@]} ${array_DC[*]} -priors 4 ${array_Priors[*]} -mask ${NameBrainmask} -out 2 ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}.nii.gz ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}.nii.gz -txt_out ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}.txt -bc_order 0 -CovPriors 8 -BFP 1 -maxRunEM 1 ${ProgMod} -AtlasSmoothing 1 1 -AtlasWeight 1 ${AW}  -SMOrder 0 -KernelSize 3 -PriorsKept 5 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4 -deleteUW 1 -outliersM 3 -outliersW 0.01 -init_splitUnif 0 -splitAccept 0 -unifTot 1 -MRF 1 -GMRF ${NameGMatrix} -juxtaCorr 1 -BiASM 0 -inDC_flag 1 -priorDGM ${PathScratch}/ICBM_DGM_${PN}_${TP[tp]}.nii.gz -BiASM 0 -maxRunEM 1 -juxtaCorr 1 ${TAString} " >> ${NameScript}
                echo "rm ${PathID}/${PN}_${TP[tp]}/Results/BG*" >> ${NameScript}
                echo "rm ${PathID}/${PN}_${TP[tp]}/Priors/GC*" >> ${NameScript}
                echo "Have to redo EMfO GC"
            fi

# Get IO atlases from that one
            if [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/IO_${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.nii.gz ]
            then
                echo "${PathSeg}/Seg_Analysis -inTxt2 ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.txt ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.nii.gz -IO 1 -mask ${NameBrainmask} ${ChangePathString}" >> ${NameScript}
            fi

            if [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/DataCorrected_GWto${TP[tp]}_${NameTogether}_${PN}.nii.gz ]
            then
                echo "${PathReg}/reg_resample -ref ${PathScratch}/DataCorrected_${NameTogether}_EMfO_${PN}_${TP[tp]}.nii.gz -flo ${NameDataCorrectedGW} -cpp ${PathScratch}/nrr_2/*cpp*${PN}_${TP[tp]}*it2*backward*.nii.gz -res  ${PathScratch}/DataCorrected_GWto${TP[tp]}_${NameTogether}_${PN}.nii.gz " >> ${NameScript}
            fi
            if [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/DataCorrected_HM_GWto${TP[tp]}_${NameTogether}_${PN}_0.nii.gz ]
            then
                if [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/I_${NameTogether}_EMfO_${PN}_${TP[tp]}_GC.nii.gz ]
                then
                    echo "${PathSeg}/Seg_Analysis -inTxt2 ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.txt ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.nii.gz -IO 1 -mask ${NameBrainmask} ${ChangePathString}" >> ${NameScript}
                fi
                echo "${PathSeg}/Seg_Analysis -maskMatch ${PathScratch}/I_${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}.nii.gz -orderFit 2 -matchRef 1 ${PathScratch}/DataCorrected_GWto${TP[tp]}_${NameTogether}_${PN}.nii.gz -matchFloat 1 ${PathScratch}/DataCorrected_${NameTogether}_EMfO_${PN}_${TP[tp]}.nii.gz -outFit ${PathScratch}/DataCorrected_HM_GWto${TP[tp]}_${NameTogether}_${PN}.nii.gz" >> ${NameScript}
            fi
        fi
        if [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.nii.gz ] || [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.txt ] [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/Data*${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.nii.gz ]
        then
            if [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/IO_${NameTogether}_EMfO_${PN}_${TP[tp]}_GC.nii.gz ]
            then
                echo "${PathSeg}/Seg_Analysis -inTxt2 ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.txt ${PathScratch}/${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.nii.gz -IO 1 -mask ${NameBrainmask} ${ChangePathString}" >> ${NameScript}
            fi

            echo "${PathSeg}/Seg_BiASM -averageGC 0.5 ${PathScratch}/GC_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz  ${PathScratch}/GC_${PN}_${TP[tp]}${OptTP}${TAOpt}.nii.gz -averageIO 0.5 ${PathScratch}/IO_${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.nii.gz ${PathScratch}/IO_${NameTogether}_EMfO_${PN}_${TP[tp]}_GC${OptTP}${TAOpt}${SignModa}.nii.gz -inDC ${PathScratch}/DataCorrected_HM_GWto${TP[tp]}_${NameTogether}_${PN}_0.nii.gz -txt_in ${PathScratch}/${NameTogether}_BiASM_${PN}_GW${OptTP}${TAOpt}${SignModa}.txt -meanPriors 1 -txt_out ${PathScratch}/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.txt -out 2 ${PathScratch}/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.nii.gz ${PathScratch}/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.nii.gz -meanPriors 1 -juxtaCorr 1 -mask ${NameBrainmask} -MRF 1 -GMRF ${NameGMatrix} -bc_order 0 -CovPriors 8 -priorDGM ${PathScratch}/ICBM_DGM_${PN}_${TP[tp]}.nii.gz -BiASM 0 -maxRunEM 1 -outliersM 3 -PriorsKept 5 -unifSplitW 0.5 -varInitUnif 1 -uniformTC 4" >> ${NameScript}

            echo "rm ${PathScratch}/MRF*" >> ${NameScript}
            echo "rm ${PathScratch}/PriorsAdapted*" >> ${NameScript}
            echo "rm ${PathScratch}/BG*" >> ${NameScript}

            if [ ! -f ${PathID}/${PN}_${TP[tp]}/${PN}_${TP[tp]}_to_ICBM.txt ]
            then
                echo "${PathReg}/reg_aladin -ref ${PathICBM}/ICBM.nii.gz -flo ${PathScratch}/T1_${PN}_${TP[tp]}.nii.gz -rigOnly -aff ${PathScratch}/${PN}_${TP[tp]}_to_ICBM.txt ">> ${NameScript}
                echo "cp ${PathScratch}/${PN}_${TP[tp]}_to_ICBM.txt ${PathID}/${PN}_${TP[tp]}/${PN}_${TP[tp]}_to_ICBM.txt" >> ${NameScript}
                echo "" >> ${NameScript}
            fi
            if [ ! -f ${PathID}/${PN}_${TP[tp]}/ICBM_${PN}_${TP[tp]}_RigTrans.txt ]
            then
                echo "${PathReg}/reg_transform -invAff ${PathScratch}/${PN}_${TP[tp]}_to_ICBM.txt ${PathScratch}/ICBM_${PN}_${TP[tp]}_RigTrans.txt" >> ${NameScript}
                echo "cp ${PathScratch}/ICBM_${PN}_${TP[tp]}_RigTrans.txt ${PathID}/${PN}_${TP[tp]}/ICBM_${PN}_${TP[tp]}_RigTrans.txt" >> ${NameScript}
                echo "" >> ${NameScript}
            fi


            if [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/LesionCorr*SP1ST1*_TO22${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.nii.gz ]
            then
                echo "${PathSeg}/Seg_Analysis -inMNI2 ${PathScratch}/ICBM_${PN}_${TP[tp]}_RigTrans.txt ${PathICBM}/ICBM.nii.gz -inTxt2 ${PathScratch}/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.txt ${PathScratch}/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.nii.gz -mask ${NameBrainmask} -Package 1 -SegType 1 -WeightedSeg 1 3 1 -connect -correct -inModa ${#ArrayModaCode[@]} ${ArrayModaCode[*]} -inRuleTxt ${RuleFileName} -WMCard 1  -inPriorsICSF ${PathScratch}/ICBM_ICSF_${PN}_${TP[tp]}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${PN}_${TP[tp]}.nii.gz -juxtaCorr 1 -TO 1 -SP 1 -oldLes 0 -inPriorsCGM ${PathScratch}/ICBM_CGM_${PN}_${TP[tp]}.nii.gz -inPriorsECSF ${PathScratch}/ICBM_ECSF_${PN}_${TP[tp]}.nii.gz ${ChangePathString}" >> ${NameScript}
            fi
        elif [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/LesionCorr*SP1ST1*_TO22${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.nii.gz ]
        then
            echo "${PathSeg}/Seg_Analysis -inMNI2 ${PathScratch}/ICBM_${PN}_${TP[tp]}_RigTrans.txt ${PathICBM}/ICBM.nii.gz -inTxt2 ${PathScratch}/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.txt ${PathScratch}/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.nii.gz -mask ${NameBrainmask} -Package 1 -SegType 1 -WeightedSeg 1 3 1 -connect -correct -inModa ${#ArrayModaCode[@]} ${ArrayModaCode[*]} -inRuleTxt ${RuleFileName} -WMCard 1  -inPriorsICSF ${PathScratch}/ICBM_ICSF_${PN}_${TP[tp]}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${PN}_${TP[tp]}.nii.gz -inPriorsCGM ${PathScratch}/ICBM_CGM_${PN}_${TP[tp]}.nii.gz -inPriorsECSF ${PathScratch}/ICBM_ECSF_${PN}_${TP[tp]}.nii.gz -juxtaCorr 1 -TO 1 -SP 1 -oldLes 0 ${ChangePathString}" >> ${NameScript}

        fi
        if((OptArt==1))
        then
            if [ ! -f ${PathID}/${PN}_${TP[tp]}/Results/LesionCorr*AC1SP1ST1*_TO22${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.nii.gz ]
            then
                echo "${PathSeg}/Seg_Analysis -inMNI2 ${PathScratch}/ICBM_${PN}_${TP[tp]}_RigTrans.txt ${PathICBM}/ICBM.nii.gz -inTxt2 ${PathScratch}/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.txt ${PathScratch}/${NameTogether}_BiASM_${PN}_${TP[tp]}_HM${OptTP}${TAOpt}.nii.gz -mask ${NameBrainmask} -Package 1 -SegType 1 -WeightedSeg 1 3 1 -connect -correct -inModa ${#ArrayModaCode[@]} ${ArrayModaCode[*]} -inRuleTxt ${RuleFileName} -WMCard 1  -inPriorsICSF ${PathScratch}/ICBM_ICSF_${PN}_${TP[tp]}.nii.gz -inPriorsDGM ${PathScratch}/ICBM_DGM_${PN}_${TP[tp]}.nii.gz -inPriorsCGM ${PathScratch}/ICBM_CGM_${PN}_${TP[tp]}.nii.gz -inPriorsECSF ${PathScratch}/ICBM_ECSF_${PN}_${TP[tp]}.nii.gz -juxtaCorr 1 -TO 1 -SP 1 -oldLes 0 -inArtefact ${PathScratch}/${PN}_${TP[tp]}_Artefacts.nii.gz " >> ${NameScript}
            fi
    
            echo "cp ${PathScratch}/LesionCorr*${PN}_${TP[tp]}* ${PathID}/${PN}_${TP[tp]}/Results/." >> ${NameScript}
            echo "cp ${PathScratch}/${NameTogether}*BiASM*${PN}_${TP[tp]}* ${PathID}/${PN}_${TP[tp]}/Results/." >> ${NameScript}
            echo "cp ${PathScratch}/${NameTogether}*EMfO*${PN}_${TP[tp]}* ${PathID}/${PN}_${TP[tp]}/Results/." >> ${NameScript}
            echo "cp ${PathScratch}/DataCorr*${NameTogether}*BiASM*${PN}_${TP[tp]}* ${PathID}/${PN}_${TP[tp]}/Results/." >> ${NameScript}
            echo "cp ${PathScratch}/DataCorr*${NameTogether}*GW*${PN}_${TP[tp]}* ${PathID}/${PN}_${TP[tp]}/Results/." >> ${NameScript}
        fi

    done

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
echo "cp ${PathScratch}/*BiASM*GW${OptTP}* ${PathID}/${PN}_GW${OptTP}/Results/." >>${NameScript}
echo "cp ${PathScratch}/*BiASM*HM${OptTP}* ${PathID}/${PN}_GW${OptTP}/Results/." >>${NameScript}
echo "cp ${PathScratch}/*GW${OptTP}* ${PathID}/${PN}_GW${OptTP}/Results/." >>${NameScript}
echo "cp ${PathScratch}/Mean* ${PathID}/${PN}_GW${OptTP}/." >> ${NameScript}
echo "cp -r ${PathScratch}/nrr_2 ${PathID}/${PN}_GW${OptTP}/Reg/." >> ${NameScript}
echo "rm -r ${PathScratch}" >> ${NameScript}

if [[ $2 == "sh" ]]
then
    sh ${NameScript}
else
    qsub -S /bin/bash -l h_rt=20:0:0 -l h_vmem=7.9G -l tmem=7.9G -l tscratch=20G -N GlobalISBI${PN} ${NameScript}
fi
