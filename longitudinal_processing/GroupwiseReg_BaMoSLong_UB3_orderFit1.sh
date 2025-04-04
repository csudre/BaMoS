#!/bin/sh

set -x
#### What could be done ##################################################################
# - add a preprocessing step in order to intensity normalise all the input images ???
# - Any other ?
##########################################################################################

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


#############################################################################
# read the input parameters
. $1
NameScriptGroupwise="${2}"
echo ${NameScriptGroupwise}
#############################################################################
## the argument value are checked
if [ ${#IMG_INPUT[@]} -lt 2 ]
then
	echo "Less than 2 images have been specified"
	echo "Exit ..."
	exit
fi

if [ ! -e ${TEMPLATE} ]
then
	echo "The template image (${TEMPLATE}) does not exist"
	echo "Exit ..."
	exit
fi

if [ "${TEMPLATE_MASK}" != "" ] && [ ! -f ${TEMPLATE_MASK} ]
then
	echo "The template image mask (${TEMPLATE_MASK}) does not exist"
	echo "Exit ..."
fi

IMG_NUMBER=${#IMG_INPUT[@]}
MASK_NUMBER=${#IMG_INPUT_MASK[@]}
INLIER_NUMBER=${#IMG_INPUT_INLIERS[@]}
if [ ${MASK_NUMBER} -gt 0 ] && [ ! -f ${IMG_INPUT_MASK[0]} ] \
	&& [ ${MASK_NUMBER} != ${IMG_NUMBER} ]
then
	echo "The number of images is different from the number of floating masks"
	echo "Exit ..."
	exit
fi

#############################################################################
## SET UP THE NIFTYREG EXECUTABLES
AFFINE=${PathReg}/reg_aladin
NRR=${PathReg}/reg_f3d
RES=${PathReg}/reg_resample
AVERAGE=${PathReg}/reg_average
TRANS=${PathReg}/reg_transform
TOOLS=${PathReg}/reg_tools
SEGMATHS=${PathSeg}/seg_maths
SEGANALYSIS=${PathBaMoS}/Seg_Analysis
SEGBIASM=${PathBaMoS}/Seg_BiASM

#############################################################################
echo ""
echo "************************************************************"
echo ">>> There are ${IMG_NUMBER} input images to groupwise register <<<"
echo ">>> The template image to initialise the registration is ${TEMPLATE} <<<"
echo "************************************************************"
echo ""
#############################################################################
# CREATE THE RESULT FOLDER
if [ ! -d ${RES_FOLDER} ]
then
	echo "The output image folder (${RES_FOLDER}) does not exist"
	mkdir ${RES_FOLDER}
	if [ ! -d ${RES_FOLDER} ]
	then
		echo "Unable to create the ${RES_FOLDER} folder"
		echo "Exit ..."
#exit
	else
		echo "The output image folder (${RES_FOLDER}) has been created"
	fi
fi

#############################################################################
#############################################################################
echo \#\!/bin/sh > ${NameScriptGroupwise}

echo "echo "Beginning Groupwise from ${NameScriptGroupwise}" " >> ${NameScriptGroupwise}
# PERFORM THE RIGID/AFFINE REGISTRATION
echo "Performing the rigid/affine registration"
# The initial average image is as specified by the user
averageImage=${TEMPLATE}
echo "${PathSeg}/seg_maths ${averageImage} -tp 0 ${RES_FOLDER}/average1moda.nii.gz" >> ${NameScriptGroupwise}
echo ${averageImage} $AFF_IT_NUM

# Loop over all iterations
for (( CUR_IT=1; CUR_IT<=${AFF_IT_NUM}; CUR_IT++ ))
do
echo "$CUR_IT" "current iteration"
	# Check if the iteration has already been performed
	if [ ! -f ${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}.nii.gz ]
	then
		#############################
		# Create a folder to store the result
		if [ ! -d ${RES_FOLDER}/aff_${CUR_IT} ]
		then
			mkdir ${RES_FOLDER}/aff_${CUR_IT}
		fi

		#############################
		# Run the rigid or affine registration
        echo "Run the rigid"
			# All registration are performed serially
        echo "Performing registration serially "
        for (( i=0 ; i<${IMG_NUMBER}; i++ ))
        do
            name=`basename ${IMG_INPUT[${i}]} .gz`
            name=`basename ${name} .nii`
            name=`basename ${name} .hdr`
            name=`basename ${name} .img`
				# Check if the registration has already been performed
            if [ ! -f ${RES_FOLDER}/aff_${CUR_IT}/aff_mat_${name}_it${CUR_IT}.txt ]
            then
                aladin_args=""
					# Registration is forced to be rigid for the first step
                if [ ${CUR_IT} == 1 ]
                then
                    aladin_args="-rigOnly"
                else
                    # Check if a previous affine can be use for initialisation
                    if [ -f ${RES_FOLDER}/aff_`expr ${CUR_IT} - 1`/aff_mat_${name}_it`expr ${CUR_IT} - 1`.txt ]
                    then
                        aladin_args="-inaff 	${RES_FOLDER}/aff_`expr ${CUR_IT} - 1`/aff_mat_${name}_it`expr ${CUR_IT} - 1`.txt"
                    fi
                fi
					# Check if a mask has been specified for the reference image
                if [ "${TEMPLATE_MASK}" != "" ]
                then
                    aladin_args="${aladin_args} -rmask ${TEMPLATE_MASK}"
                fi
                # Check if a mask has been specified for the floating image
                if [ ${MASK_NUMBER} == ${IMG_NUMBER} ]
                then
                    aladin_args="${aladin_args} -fmask ${IMG_INPUT_MASK[${i}]}"
                fi
                result="/dev/null"
                if [ "${CUR_IT}" == "${AFF_IT_NUM}" ]
                then
                    result="${RES_FOLDER}/aff_${CUR_IT}/aff_res_${name}_it${CUR_IT}.nii.gz"
                fi
                # Perform the registration
                echo "${PathReg}/reg_aladin ${AFFINE_args} ${aladin_args} \
                    -ref ${averageImage} \
                    -flo ${IMG_INPUT[${i}]} \
                    -aff ${RES_FOLDER}/aff_${CUR_IT}/aff_mat_${name}_it${CUR_IT}.txt \
                    -res ${result} > ${RES_FOLDER}/aff_${CUR_IT}/aff_log_${name}_it${CUR_IT}.txt" >> ${NameScriptGroupwise}
                echo "${PathReg}/reg_transform -ref ${averageImage} -def ${RES_FOLDER}/aff_${CUR_IT}/aff_mat_${name}_it${CUR_IT}.txt ${RES_FOLDER}/aff_${CUR_IT}/affdef_mat_${name}_it${CUR_IT}.nii.gz" >> ${NameScriptGroupwise}
		if [ ! -f ${RES_FOLDER}/aff_${CUR_IT}/aff_mat_${name}_it${CUR_IT}.txt ]
                then
                    echo "Error when creating \
                        ${RES_FOLDER}/aff_${CUR_IT}/aff_mat_${name}_it${CUR_IT}.txt"
#exit
                fi
            fi
        done


		#############################
        echo ${IMG_INPUT[0]}
		if [ "${CUR_IT}" != "${AFF_IT_NUM}" ]
		then
		# The transformation are demean'ed to create the average image
		# Note that this is not done for the last iteration step
			list_average=""
			list_names=""
			for img in ${IMG_INPUT[@]}
			do
				name=`basename ${img} .gz`
				name=`basename ${name} .nii`
				name=`basename ${name} .hdr`
				name=`basename ${name} .img`
				list_names="${list_names} ${name}"
				list_average="${list_average} \
					${RES_FOLDER}/aff_${CUR_IT}/affdef_mat_${name}_it${CUR_IT}.nii.gz ${img}"
			done

# Creation of the list of transformation to apply on the Inlier masks as well
			array_names=()
			#echo ${list_names}
			read -a array_names <<< ${list_names}
			#echo "array_names is " ${array_names[*]}
			list_inlier=""
           		list_average_inlier=""
			for ((i=0;i<${#IMG_INPUT_INLIERS[@]};i++))
			do
			list_inlier="${list_inlier} ${RES_FOLDER}/aff_${CUR_IT}/aff_mat_${array_names[i]}_it${CUR_IT}.txt ${IMG_INPUT_INLIERS[i]}"
list_average_inlier="${list_average_inlier} ${RES_FOLDER}/aff_${CUR_IT}/affdef_mat_${array_names[i]}_it${CUR_IT}.nii.gz ${IMG_INPUT_INLIERS[i]}"
			done
            echo "list of inliers is " ${list_average_inlier}


                # The average is created on the host
echo "${PathReg}/reg_average ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}.nii.gz  -demean1 ${RES_FOLDER}/average1moda.nii.gz ${list_average_inlier}" >> ${NameScriptGroupwise}
#echo "${PathReg}/reg_average -demean ${RES_FOLDER}/average_inliers_affine_it_${CUR_IT}.nii.gz ${list_average_inlier}" >> ${NameScriptGroupwise}
# The average is created on the host
echo "${PathReg}/reg_average  ${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}_temp.nii.gz  -demean1 ${averageImage}  ${list_average}" >> ${NameScriptGroupwise}



########### IF MEAN AND HISTOGRAM MATCHING, WOULD HAVE TO BE INSERTED THERE BUT THEN NOT POSSIBLE TO USE DEMEAN1: NEED TO SEPARATE BETWEEN REG_RESAMPLE AND REG_AVERAGE
#echo ${IMG_INPUT[0]}
            echo "Going to perform the histogram matching "
# First do the resampling and create the lists for possible averaging of masks and resampled images
            array_resampled=()
            list_resampled=""
            array_average=()
            array_inlier=()


            list_resampledMaskInliers=""
            read -a array_average <<< ${list_average}
            echo "length of list_average is ${#array_average[@]}"
            read -a array_inlier <<< ${list_inlier}
            echo "list of inliers is " ${array_inlier[*]}
            echo "${array_average[*]} is array average"
            for ((i=0;i<${#array_average[@]};i=i+2))
            do

                echo "${PathReg}/reg_resample -ref ${averageImage} -flo ${array_average[i+1]} -aff ${array_average[i]} -res ${RES_FOLDER}/aff_${CUR_IT}/resampled_affine_it_${CUR_IT}_${i}.nii.gz" >> ${NameScriptGroupwise}
                list_resampled="${list_resampled} \ ${RES_FOLDER}/aff_${CUR_IT}/resampled_affine_it_${CUR_IT}_${i}.nii.gz"

                echo "${array_inlier[i]}"


                echo "${PathReg}/reg_resample -ref ${averageImage} -flo ${array_inlier[i+1]} -aff ${array_inlier[i]} -res ${RES_FOLDER}/aff_${CUR_IT}/inliers_affine_it_${CUR_IT}_${i}.nii.gz" >> ${NameScriptGroupwise}
                list_resampledMaskInliers="${list_resampledMaskInliers} \ ${RES_FOLDER}/aff_${CUR_IT}/inliers_affine_it_${CUR_IT}_${i}.nii.gz "
            done

            echo ${list_resampledMaskInliers}
            read -a array_resampled <<< ${list_resampled}
            array_Maskresampled=()
            read -a array_Maskresampled <<< ${list_resampledMaskInliers}
            echo ${array_Maskresampled[*]}

#echo "${PathReg}/reg_average ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}.nii.gz -avg ${array_Maskresampled[*]}"

#echo "Doing seg_maths averaging"
echo "${PathBaMoS}/Seg_Analysis -meanImages ${#array_Maskresampled[@]} ${array_Maskresampled[*]} -outMean ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}.nii.gz" >> ${NameScriptGroupwise}


            let "NM=${#array_Maskresampled[@]}-1"

            echo "${PathSeg}/seg_maths ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}.nii.gz -thr 0.9 -bin -odt -char ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}_bin.nii.gz" >> ${NameScriptGroupwise}

# Second, perform the intensity matching based on inlier masks and get the images out. Take the first image in the list as reference point.
            list_matched=""
            echo "Doing the matching... for"

#echo "${PathBaMoS}/Seg_Analysis -meanImages ${#array_resampled[@]} ${array_resampled[*]} -outMean ${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}_temp.nii.gz" >> ${NameScriptGroupwise}
             averageImage=${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}_temp.nii.gz
            echo "${PathSeg}/seg_maths ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}_bin.nii.gz -dil 10 -fill ${RES_FOLDER}/aff_${CUR_IT}/average_mask.nii.gz " >>${NameScriptGroupwise}

            for ((i=0;i<${#array_resampled[@]};i++))
            do
#echo "Seg_Analysis -maskMatch ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}_bin.nii.gz -orderFit 1 -matchRef 1${array_resampled[0]} -matchFloat 1 ${array_resampled[i]} -outFit ${RES_FOLDER}/aff_${CUR_IT}/average_matched_affine_it_${CUR_IT}_${i}.nii.gz"

                echo "${PathBaMoS}/Seg_Analysis -mask ${RES_FOLDER}/aff_${CUR_IT}/average_mask.nii.gz -maskMatch ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}_bin.nii.gz -orderFit 1 -matchRef 1 ${averageImage} -matchFloat 1 ${array_resampled[i]} -outFit ${RES_FOLDER}/aff_${CUR_IT}/average_matched_affine_it_${CUR_IT}_${i}.nii.gz" >> ${NameScriptGroupwise}

#Seg_Analysis -maskMatch ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}_bin.nii.gz -orderFit 2 -matchRef ${array_resampled[0]} -matchFloat ${array_resampled[i]} -outFit ${RES_FOLDER}/aff_${CUR_IT}/average_matched_affine_it_${CUR_IT}_${i}.nii.gz

                list_matched="${list_matched} \ ${RES_FOLDER}/aff_${CUR_IT}/average_matched_affine_it_${CUR_IT}_${i}_0.nii.gz "
            done

# Realise average on intensity matched images
            array_matched=()
            read -a array_matched <<< ${list_matched}
            echo "${PathBaMoS}/Seg_Analysis -meanImages ${#array_matched[@]} ${array_matched[*]} -outMean ${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}.nii.gz" >> ${NameScriptGroupwise}


#################################################################

# The average is created on the host reg_average \ ${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}.nii.gz \ -demean1 ${averageImage} \ ${list_average}


            if [ ! -f ${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}.nii.gz ]
            then
                echo "Error when creating \
                    ${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}.nii.gz"
#exit
            fi
	
		else # if [ "${CUR_IT}" == "${AFF_IT_NUM}" ]
			# All the result images are directly averaged during the last step

		# We still need list average and everything for the transformation of the masks...
			list_average=""
			list_names=""
			for img in ${IMG_INPUT[@]}
			do
				name=`basename ${img} .gz`
				name=`basename ${name} .nii`
				name=`basename ${name} .hdr`
				name=`basename ${name} .img`
				list_names="${list_names} ${name}"
				list_average="${list_average} \
					${RES_FOLDER}/aff_${CUR_IT}/aff_mat_${name}_it${CUR_IT}.txt ${img}"
			done

# Creation of the list of transformation to apply on the Inlier masks as well
			array_names=()
			#echo ${list_names}
			read -a array_names <<< ${list_names}
			#echo "array_names is " ${array_names[*]}
			list_inlier=""
            list_resampled=""
			for ((i=0;i<${#IMG_INPUT_INLIERS[@]};i++))
			do
				list_inlier="${list_inlier} \
					${RES_FOLDER}/aff_${CUR_IT}/aff_mat_${array_names[i]}_it${CUR_IT}.txt ${IMG_INPUT_INLIERS[i]}"
                list_resampled="${list_resampled} ${RES_FOLDER}/aff_${CUR_IT}/aff_res_${array_names[i]}_it${CUR_IT}.nii.gz"
			done
            echo "list of inliers is " ${list_inlier}


########### IF MEAN AND HISTOGRAM MATCHING, WOULD HAVE TO BE INSERTED THERE BUT THEN NOT POSSIBLE TO USE DEMEAN1: NEED TO SEPARATE BETWEEN REG_RESAMPLE AND REG_AVERAGE
#echo ${IMG_INPUT[0]}
            echo "Going to perform the histogram matching "
# First do the resampling and create the lists for possible averaging of masks and resampled images
            array_resampled=()
# list_resampled=`ls ${RES_FOLDER}/aff_${CUR_IT}/aff_res_*_it${CUR_IT}.nii*`
            array_average=()
            array_inlier=()


            list_resampledMaskInliers=""
            read -a array_average <<< ${list_average}
#echo ${array_average[1]}
            echo "length of list_average is ${#array_average[@]}"
            read -a array_inlier <<< ${list_inlier}
            echo "list of inliers is " ${array_inlier[*]}

            for ((i=0;i<${#array_average[@]};i=i+2))
            do
#echo "reg_resample -ref ${averageImage} -flo ${array_average[i+1]} -aff ${array_average[i]} -res ${RES_FOLDER}/aff_${CUR_IT}/resampled_affine_it_${CUR_IT}_${i}.nii.gz"
                echo "${PathReg}/reg_resample -ref ${averageImage} -flo ${array_average[i+1]} -aff ${array_average[i]} -res ${RES_FOLDER}/aff_${CUR_IT}/resampled_affine_it_${CUR_IT}_${i}.nii.gz" >> ${NameScriptGroupwise}
#list_resampled="${list_resampled} \ ${RES_FOLDER}/aff_${CUR_IT}/resampled_affine_it_${CUR_IT}_${i}.nii.gz"

#echo "list of inliers is in loop " ${array_inlier[*]}
                echo "${array_inlier[i]}"



                echo "${PathReg}/reg_resample -ref ${averageImage} -flo ${array_inlier[i+1]} -aff ${array_inlier[i]} -res ${RES_FOLDER}/aff_${CUR_IT}/inliers_affine_it_${CUR_IT}_${i}.nii.gz" >> ${NameScriptGroupwise}
                list_resampledMaskInliers="${list_resampledMaskInliers} \ ${RES_FOLDER}/aff_${CUR_IT}/inliers_affine_it_${CUR_IT}_${i}.nii.gz "
            done

            echo ${list_resampledMaskInliers}
            read -a array_resampled <<< ${list_resampled}
            array_Maskresampled=()
            read -a array_Maskresampled <<< ${list_resampledMaskInliers}
            echo ${array_Maskresampled[*]}

            echo "reg_average ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}.nii.gz -avg ${array_Maskresampled[*]}"

            echo "Doing seg_maths averaging"
            echo "${PathBaMoS}/Seg_Analysis -meanImages ${#array_Maskresampled[@]} ${array_Maskresampled[*]} -outMean ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}.nii.gz" >> ${NameScriptGroupwise}

            let "NM=${#array_Maskresampled[@]}-1"

            echo "${PathSeg}/seg_maths ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}.nii.gz -thr 0.9 -bin -odt -char ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}_bin.nii.gz" >> ${NameScriptGroupwise}

# Second, perform the intensity matching based on inlier masks and get the images out. Take the first image in the list as reference point.
            list_matched=" "
            echo "Doing the matching... for"

            echo "${PathBaMoS}/Seg_Analysis -meanImages ${#array_resampled[@]} ${array_resampled[*]} -outMean ${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}_temp.nii.gz" >> ${NameScriptGroupwise}
             averageImage=${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}_temp.nii.gz
             echo "${PathSeg}/seg_maths ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}_bin.nii.gz -dil 10 -fill ${RES_FOLDER}/aff_${CUR_IT}/average_mask.nii.gz " >>${NameScriptGroupwise}

            for ((i=0;i<${#array_resampled[@]};i++))
            do
                echo "Seg_Analysis -mask ${RES_FOLDER}/aff_${CUR_IT}/average_mask.nii.gz -maskMatch ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}_bin.nii.gz -orderFit 1 -matchRef 1 ${array_resampled[0]} -matchFloat 1 ${array_resampled[i]} -outFit ${RES_FOLDER}/aff_${CUR_IT}/average_matched_affine_it_${CUR_IT}_${i}.nii.gz"

                echo "${PathBaMoS}/Seg_Analysis -mask ${RES_FOLDER}/aff_${CUR_IT}/average_mask.nii.gz -maskMatch ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}_bin.nii.gz -orderFit 1 -matchRef 1 ${averageImage} -matchFloat 1 ${array_resampled[i]} -outFit ${RES_FOLDER}/aff_${CUR_IT}/average_matched_affine_it_${CUR_IT}_${i}.nii.gz" >> ${NameScriptGroupwise}

                list_matched="${list_matched} \ ${RES_FOLDER}/aff_${CUR_IT}/average_matched_affine_it_${CUR_IT}_${i}_0.nii.gz "
            done

# Realise average on intensity matched images
            array_matched=()
            read -a array_matched <<< ${list_matched}
            echo "${PathBaMoS}/Seg_Analysis -meanImages ${#array_matched[@]} ${array_matched[*]} -outMean ${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}.nii.gz" >> ${NameScriptGroupwise}

            echo "Final affine averaging done"
#################################################################


				#reg_average \
					${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}.nii.gz \
					-avg \
					`ls ${RES_FOLDER}/aff_${CUR_IT}/aff_res_*_it${CUR_IT}.nii*`
				if [ ! -f ${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}.nii.gz ]
				then
					echo "Error when creating \
						${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}.nii.gz"
#exit
				fi
        fi
#fi # if [ "${CUR_IT}" != "${AFF_IT_NUM}" ]
	else # if [ ! -f ${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}.nii.gz ]
		echo "${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}.nii.gz already exists"
	fi # if [ ! -f ${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}.nii.gz ]
	# Update the average image used as a reference
	averageImage=${RES_FOLDER}/aff_${CUR_IT}/average_affine_it_${CUR_IT}.nii.gz
	echo "${PathSeg}/seg_maths ${averageImage} -tp 0 ${RES_FOLDER}/average1moda.nii.gz" >> ${NameScriptGroupwise}
done # Loop over affine iteration


#############################################################################
#############################################################################
### Non rigid registration loop

echo "Beginning non rigid registration loop"

for (( CUR_IT=1; CUR_IT<=${NRR_IT_NUM}; CUR_IT++ ))
do
    echo "Treating $CUR_IT"
	#############################
	# Check if the current average image has already been created
	if [ ! -f ${RES_FOLDER}/nrr_${CUR_IT}/average_nonrigid_it_${CUR_IT}.nii.gz ]
	then

		#############################
		# Create a folder to store the current results
		if [ ! -d ${RES_FOLDER}/nrr_${CUR_IT} ]
		then
			mkdir ${RES_FOLDER}/nrr_${CUR_IT}
		fi

		#############################
		# Run the nonrigid registrations
        echo "Run non rigid"
        echo "${IMG_NUMBER}"
			for (( i=0 ; i<${IMG_NUMBER}; i++ ))
			do
				name=`basename ${IMG_INPUT[i]} .gz`
				name=`basename ${name} .nii`
				name=`basename ${name} .hdr`
				name=`basename ${name} .img`
				# Check if the registration has already been performed
				if [ ! -f ${RES_FOLDER}/nrr_${CUR_IT}/nrr_cpp_${name}_it${CUR_IT}.nii* ]
                    echo "need to do it"
				then
					f3d_args=""
					# Check if a mask has been specified for the reference image
					if [ "${TEMPLATE_MASK}" != "" ]
					then
						f3d_args="${f3d_args} -rmask ${TEMPLATE_MASK}"
					fi
					# Check if a mask has been specified for the floating image
					if [ ${MASK_NUMBER} == ${IMG_NUMBER} ]
					then
						f3d_args="${f3d_args} -fmask ${IMG_INPUT_MASK[i]}"
					fi
					if [ ${AFF_IT_NUM} -gt 0 ]
					then
						f3d_args="${f3d_args} -aff \
							${RES_FOLDER}/aff_${AFF_IT_NUM}/aff_mat_${name}_it${AFF_IT_NUM}.txt"
					fi
					result="/dev/null"
					if [ "${CUR_IT}" == "${NRR_IT_NUM}" ]
					then
						result="${RES_FOLDER}/nrr_${CUR_IT}/nrr_res_${name}_it${CUR_IT}.nii.gz"
					fi
					# Perform the registration
                        echo "${PathReg}/reg_f3d ${NRR_args} ${f3d_args} -ref ${averageImage} \ -flo ${IMG_INPUT[${i}]} \ -cpp ${RES_FOLDER}/nrr_${CUR_IT}/nrr_cpp_${name}_it${CUR_IT}.nii.gz \ -res ${result} > ${RES_FOLDER}/nrr_${CUR_IT}/nrr_log_${name}_it${CUR_IT}.txt"
                        echo " ${PathReg}/reg_f3d ${NRR_args} ${f3d_args} \
						-ref ${averageImage} \
						-flo ${IMG_INPUT[${i}]} \
						-cpp ${RES_FOLDER}/nrr_${CUR_IT}/nrr_cpp_${name}_it${CUR_IT}.nii.gz \
						-res ${result} > ${RES_FOLDER}/nrr_${CUR_IT}/nrr_log_${name}_it${CUR_IT}.txt " >> ${NameScriptGroupwise}
				fi
			done
				#############################


# Get the inlier list for the mask. We need it for all of the transformation even the last one in the loop
			list_inlier=""
			if [ "${INLIER_NUMBER}" == "${IMG_NUMBER}" ]
			then
				for ((i=0;i<${#IMG_INPUT[@]};i++))
				do
					name=`basename ${IMG_INPUT[i]} .gz`
					name=`basename ${name} .nii`
					name=`basename ${name} .hdr`
					name=`basename ${name} .img`
					list_inlier="${list_inlier} \
					${RES_FOLDER}/aff_${AFF_IT_NUM}/aff_mat_${name}_it${AFF_IT_NUM}.txt \
					${RES_FOLDER}/nrr_${CUR_IT}/nrr_cpp_${name}_it${CUR_IT}.nii.gz ${IMG_INPUT_INLIERS[i]}"
				done
            else
			echo "not same number inlier ing"
            fi
			list_average=""
			for img in ${IMG_INPUT[@]}
			do
				name=`basename ${img} .gz`
				name=`basename ${name} .nii`
				name=`basename ${name} .hdr`
				name=`basename ${name} .img`
				list_average="${list_average} \
					${RES_FOLDER}/aff_${AFF_IT_NUM}/aff_mat_${name}_it${AFF_IT_NUM}.txt \
					${RES_FOLDER}/nrr_${CUR_IT}/nrr_cpp_${name}_it${CUR_IT}.nii.gz ${img}"
			done



		# The transformation are demean'ed to create the average image
		# Note that this is not done for the last iteration step
		if [ "${CUR_IT}" != "${NRR_IT_NUM}" ]
		then
			list_average=""
			for img in ${IMG_INPUT[@]}
			do
				name=`basename ${img} .gz`
				name=`basename ${name} .nii`
				name=`basename ${name} .hdr`
				name=`basename ${name} .img`
				list_average="${list_average} \
					${RES_FOLDER}/aff_${AFF_IT_NUM}/aff_mat_${name}_it${AFF_IT_NUM}.txt \
					${RES_FOLDER}/nrr_${CUR_IT}/nrr_cpp_${name}_it${CUR_IT}.nii.gz ${img}"
			done
			
				# The average is created on the host
				echo " ${PathReg}/reg_average \
					${RES_FOLDER}/nrr_${CUR_IT}/average_nonrigid_it_${CUR_IT}.nii.gz \
					-demean3 ${averageImage} \
					${list_average} " >> ${NameScriptGroupwise}
				if [ ! -f ${RES_FOLDER}/nrr_${CUR_IT}/average_nonrigid_it_${CUR_IT}.nii.gz ]
				then
					echo "Error when creating \
						${RES_FOLDER}/nrr_${CUR_IT}/average_nonrigid_it_${CUR_IT}.nii.gz"
#exit
				fi
			
        # if [ "${CUR_IT}" == "${NRR_IT_NUM}" ]
			# All the result images are directly averaged during the last step

        else
###################### INCORPORATION OF THE HISTOGRAM MATCHING PROCESS AGAIN AT THIS STAGE BASED ON THE INLIER SEGMENTATION ############################

########### IF MEAN AND HISTOGRAM MATCHING, WOULD HAVE TO BE INSERTED THERE BUT THEN NOT POSSIBLE TO USE DEMEAN1: NEED TO SEPARATE BETWEEN REG_RESAMPLE AND REG_AVERAGE

# First do the resampling and create the lists for possible averaging of masks and resampled images

            list_resampled=""
            for ((i=0;i<${#IMG_INPUT[@]};i++))
            do
                name=`basename ${IMG_INPUT[i]} .gz`
                name=`basename ${name} .nii`
                name=`basename ${name} .hdr`
                name=`basename ${name} .img`
                list_resampled="${list_resampled} \ ${RES_FOLDER}/nrr_${CUR_IT}/nrr_res_${name}_it${CUR_IT}.nii.gz "
            done
#list_resampled=`ls ${RES_FOLDER}/nrr_${CUR_IT}/nrr_res_*_it${CUR_IT}.nii*`
            list_mask=""
            array_resampled=()
            read -a array_resampled <<< ${list_resampled}
            array_average=()
            array_inlier=()
            list_resampledMaskInliers=""
            read -a array_average <<< ${list_average}
#echo ${array_average[1]}
            echo "length of list_average is ${#array_average[@]}"
            read -a array_inlier <<< ${list_inlier}
            echo "list of inliers is " ${array_inlier[*]}

#Create a mask for all of the resampled image 
            echo "doing mask creation"
            for ((i=0;i<${IMG_NUMBER};i++))
            do
                echo "${PathSeg}/seg_maths ${array_resampled[i]} -tp 0 -thr 0 -bin -dil 1 -fill -ero 1 -odt char ${RES_FOLDER}/nrr_${CUR_IT}/mask_f3d_it_${CUR_IT}_${i}.nii.gz" >> ${NameScriptGroupwise}
                list_mask="${list_mask} \ ${RES_FOLDER}/nrr_${CUR_IT}/mask_f3d_it_${CUR_IT}_${i}.nii.gz"
            done

            array_mask=()
            read -a array_mask <<< ${list_mask}

#Create the resampled inlier masks for the further histogram matching
            list_resampledMaskInliers=""

            for ((i=0;i<${#array_inlier[@]};i=i+3))
            do
                echo "${PathReg}/reg_resample -ref ${averageImage} -flo ${array_inlier[i+2]} -cpp ${array_inlier[i+1]} -res ${RES_FOLDER}/nrr_${CUR_IT}/inliers_f3d_it_${CUR_IT}_${i}.nii.gz" >> ${NameScriptGroupwise}
                list_resampledMaskInliers="${list_resampledMaskInliers} \ ${RES_FOLDER}/nrr_${CUR_IT}/inliers_f3d_it_${CUR_IT}_${i}.nii.gz "
            done

            array_Maskresampled=()
            read -a array_Maskresampled <<< ${list_resampledMaskInliers}

            echo "Doing seg_maths averaging"
            echo "${PathBaMoS}/Seg_Analysis -meanImages ${#array_Maskresampled[@]} ${array_Maskresampled[*]} -outMean ${RES_FOLDER}/nrr_${CUR_IT}/average_inliers_f3d_it_${CUR_IT}.nii.gz" >> ${NameScriptGroupwise}

            let "NM=${#array_Maskresampled[@]}-1"
#echo "seg_maths ${array_Maskresampled[0]} -merge ${NM} 4 ${array_Maskresampled[@]:1} ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}_merged.nii.gz"
#seg_maths ${array_Maskresampled[0]} -merge ${NM} 4 ${array_Maskresampled[@]:1} ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}_merged.nii.gz
#seg_maths ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}_merged.nii.gz -tmean ${RES_FOLDER}/aff_${CUR_IT}/average_inliers_affine_it_${CUR_IT}.nii.gz





#reg_average ${RES_FOLDER}/nrr_${CUR_IT}/average_inliers_f3d_it_${CUR_IT}.nii.gz -avg ${list_resampledMaskInliers}
            echo "${PathSeg}/seg_maths ${RES_FOLDER}/nrr_${CUR_IT}/average_inliers_f3d_it_${CUR_IT}.nii.gz -thr 0.9 -bin -odt -char ${RES_FOLDER}/nrr_${CUR_IT}/average_inliers_f3d_it_${CUR_IT}_bin.nii.gz" >> ${NameScriptGroupwise}
echo "${PathSeg}/seg_maths ${RES_FOLDER}/nrr_${CUR_IT}/average_inliers_f3d_it_${CUR_IT}_bin.nii.gz -dil 10 -fill ${RES_FOLDER}/nrr_${CUR_IT}/average_mask.nii.gz " >>${NameScriptGroupwise}

            echo "${PathBaMoS}/Seg_Analysis -meanImages ${#array_resampled[@]} ${array_resampled[*]} -outMean ${RES_FOLDER}/nrr_${CUR_IT}/average_nonrigid_it_${CUR_IT}_temp.nii.gz  " >> ${NameScriptGroupwise}
            averageImage=${RES_FOLDER}/nrr_${CUR_IT}/average_nonrigid_it_${CUR_IT}_temp.nii.gz

# Second, perform the intensity matching based on inlier masks and get the images out. Take the first image in the list as reference point.
            list_matched="${array_resampled[0]}"
            list_matched=""

            for ((i=0;i<${#array_resampled[@]};i++))
            do
                echo "${PathBaMoS}/Seg_Analysis -mask ${RES_FOLDER}/nrr_${CUR_IT}/average_mask.nii.gz  -maskMatch ${RES_FOLDER}/nrr_${CUR_IT}/average_inliers_f3d_it_${CUR_IT}_bin.nii.gz -orderFit 1 -matchRef 1 ${averageImage} -matchFloat 1 ${array_resampled[i]} -outFit ${RES_FOLDER}/nrr_${CUR_IT}/average_matched_f3d_it_${CUR_IT}_${i}.nii.gz" >> ${NameScriptGroupwise}
                list_matched="${list_matched} \ ${RES_FOLDER}/nrr_${CUR_IT}/average_matched_f3d_it_${CUR_IT}_${i}_0.nii.gz  "
            done

# Realise average on intensity matched images
#reg_average \	${RES_FOLDER}/nrr_${CUR_IT}/average_nonrigid_it_${CUR_IT}.nii.gz \ -avg \ ${list_matched}

# Realise average on intensity matched images
            array_matched=()
            echo "${list_matched} is list matched"
            read -a array_matched <<< ${list_matched}
            echo "${PathBaMoS}/Seg_Analysis -meanImages ${#array_matched[@]} ${array_matched[*]} -outMean ${RES_FOLDER}/nrr_${CUR_IT}/average_nonrigid_it_${CUR_IT}.nii.gz " >> ${NameScriptGroupwise}
            echo "echo "Normally average done for this iteration non rigid ${CUR_IT}" ">> ${NameScriptGroupwise}


#################################################################



				#reg_average \
					${RES_FOLDER}/nrr_${CUR_IT}/average_nonrigid_it_${CUR_IT}.nii.gz \
					-avg \
					`ls ${RES_FOLDER}/nrr_${CUR_IT}/nrr_res_*_it${CUR_IT}.nii*`
				if [ ! -f ${RES_FOLDER}/nrr_${CUR_IT}/average_nonrigid_it_${CUR_IT}.nii.gz ]
				then
					echo "Error when creating \
						${RES_FOLDER}/nrr_${CUR_IT}/average_nonrigid_it_${CUR_IT}.nii.gz"
#exit
				fi
			
		fi # if [ "${CUR_IT}" != "${NRR_IT_NUM}" ]
	else # if [ ! -f ${RES_FOLDER}/nrr_${CUR_IT}/average_nonrigid_it_${CUR_IT}.nii.gz ]
		echo "${RES_FOLDER}/nrr_${CUR_IT}/average_nonrigid_it_${CUR_IT}.nii.gz already exists"
	fi # if [ ! -f ${RES_FOLDER}/nrr_${CUR_IT}/average_nonrigid_it_${CUR_IT}.nii.gz ]
	# Update the average image
	averageImage=${RES_FOLDER}/nrr_${CUR_IT}/average_nonrigid_it_${CUR_IT}.nii.gz
    echo "Trying to go to next iteration"
done
#############################################################################
