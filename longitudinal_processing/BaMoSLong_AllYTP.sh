Name=${1}
NameR=${2}
Name2=${3}
j=0
echo $Name $NameR $Name2
ArrayT1=($(ls ${Name}/*[0-9]/T1*${Name2}.nii.gz))
ArrayFLAIR=($(ls ${Name}/*[0-9]/FLAIR*${Name2}.nii.gz))
ArrayGIF=()
ArrayTP=()
for k in ${ArrayT1[@]}
do
ArrayGIF=(${ArrayGIF[@]} ${k%.nii.gz})
TP=${k%${Name2}.nii.gz}
TP=${TP%_}
TP=${TP##*${Name}_}
ArrayTP=(${ArrayTP[@]} $TP)
done

echo ${#ArrayFLAIR[@]} numb of FLAIR
Les=($(ls ${Name}/GWResults/Correct*1Les*))
if [[ ${#ArrayFLAIR[@]} -gt ${#Les[@]} ]] && [ -f ${ArrayGIF[0]}/*Parc* ] && [[ ${#ArrayFLAIR[@]} == ${#ArrayT1[@]} ]] && [[ ${#ArrayFLAIR[@]} -gt 1 ]]
then


StringFLAIR=""
StringT1=""
StringGIF=""
StringTP=""
for ((p=0;p<${#ArrayFLAIR[@]};p++))
do
StringFLAIR="${StringFLAIR} $(pwd)/${ArrayFLAIR[p]}"
StringT1="$StringT1 $(pwd)/${ArrayT1[p]}"
StringGIF="$StringGIF $(pwd)/${ArrayGIF[p]}"
StringTP="$StringTP ${ArrayTP[p]}"
done

sh ~/Scripts/BaMoSLong_CompleteCL.sh ${Name}${NameR} "\($StringFLAIR\)" "\($StringT1\)" "\($StringGIF\)" "\($StringTP\)" $(pwd)/${Name}/GWResults${NameR} 1 0 GW GW 1

#echo sh ~/Scripts/BaMoSLong_CompleteCL.sh ${Name}_${NameR} "\($(pwd)/${ArrayFLAIR[0]} $(pwd)/${ArrayFLAIR[1]} $(pwd)/${ArrayFLAIR[2]} $(pwd)/${ArrayFLAIR[3]}\)" "\($(pwd)/${ArrayT1[0]} $(pwd)/${ArrayT1[1]} $(pwd)/${ArrayT1[2]} $(pwd)/${ArrayT1[3]}\)" "\($(pwd)/${ArrayGIF[0]} $(pwd)/${ArrayGIF[1]} $(pwd)/${ArrayGIF[2]} ${ArrayGIF[3]}\)" "\(${ArrayTP[0]} ${ArrayTP[1]} ${ArrayTP[2]} ${ArrayTP[3]}\)" $(pwd)/${Name}/GWResults${NameR} 1 0 GW GW 1 

else

echo "Impossible processing"

if [[ ${#ArrayFLAIR[@]} -lt ${#ArrayT1[@]} ]]

then
ArrayGIF=()
ArrayT1=()
ArrayTP=()
for k in ${ArrayFLAIR[@]}
do
Array=(${k//// })
ArrayGIF=(${ArrayGIF[@]} ${k%/*}/T1_${Array[1]})
ArrayT1=(${ArrayT1[@]} ${k%/*}/T1_${Array[1]}.nii.gz)
TP=${Array[1]##*${Name}_}
ArrayTP=(${ArrayTP[@]} $TP)


done


for ((p=0;p<${#ArrayFLAIR[@]};p++))
do
StringFLAIR="${StringFLAIR} $(pwd)/${ArrayFLAIR[p]}"
StringT1="$StringT1 $(pwd)/${ArrayT1[p]}"
StringGIF="$StringGIF $(pwd)/${ArrayGIF[p]}"
StringTP="$StringTP ${ArrayTP[p]}"
done

sh ~/Scripts/BaMoSLong_CompleteCL.sh ${Name}${NameR} "\($StringFLAIR\)" "\($StringT1\)" "\($StringGIF\)" "\($StringTP\)" $(pwd)/${Name}/GWResults${NameR} 1 0 GW GW 1
fi
#echo $(ls ${ArrayGIF[0]})
#echo $(ls ${ArrayGIF[1]})
 
fi;


#echo $IDmin > LastTreatedBaMoS2.txt; j=0;
