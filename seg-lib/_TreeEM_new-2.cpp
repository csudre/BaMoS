#include "_TreeEM_new.h"

//
//  TreeEM.cpp
//  TreeSeg
//
//  Created by Carole Sudre on 08/04/2013.
//  Copyright (c) 2013 Carole Sudre. All rights reserved.
//


TreeEM::TreeEM(){

    this->Parent = NULL;
    this->DataImage=NULL;
    this->Mask = NULL;
    this->Priors = NULL;
    this->PriorsAdapted=NULL;
    //this->Short_2_Long_Indices = NULL;
    //this->Long_2_Short_Indices;
    this->L2S=NULL;
    this->S2L=NULL;
    this->NormResp=NULL;
    this->NormWeight=1;
    this->IndFactor=0;
    this->NumberMaskedElements=0;
    this->DPChildren=NULL;
    this->FlagDistClassInd=0;
    this->FlagOutliers=0;
    this->ParametersDistribution=NULL;
    this->SplitCheck=NULL;
    this->MergeCheck=NULL;
    // this->BFBasisFunctions=NULL;
    //    this->BFCorrection=NULL;
    this->DataBFCorrected=NULL;
    this->BFCoeffs=NULL;
    this->HardSeg=NULL;
    this->MRF=NULL;
    this->GMatrix=NULL;

}

static int numbbins=100;
static int BForder=4;
static bool BFFlag=1;
static int KernelSize=3;
static int NumbMaxLeaves=15;
static int MaxNumbLeavesperClass=15;
static PrecisionTYPE Threshold=1E-5;
static int MaxIteration=50;
static int MinIteration=6;


TreeEM::TreeEM(nifti_image * DataInput){

    this->Parent=NULL;
    this->DataImage=DataInput;
    this->Mask=NULL;
    this->L2S=NULL;
    this->S2L=NULL;
    this->Priors=NULL;
    this->HardSeg=NULL;

}

float * TreeEM::LogGaussBlur(float * GaussianBlur,int TotalSize){
    float * LogGaussBlur=new float[TotalSize];
    for(int i=0;i<TotalSize;i++){
        if(GaussianBlur[i]<=0){
            LogGaussBlur[i]=-10E6;
        }
        else{
            LogGaussBlur[i]=logf(GaussianBlur[i]);
        }
    }
    return LogGaussBlur;
}

inline float TreeEM::pow_int(const float base,
                             int exp){
    if(exp==0){return 1;}
    float result = base;
    while (--exp){result *= base;}
    return result;
}

float * TreeEM::TransposeMatrix(float * MatrixToTranspose, int SizeM,int SizeN){
    float * TransposedMatrix =new float [SizeN*SizeM];
    for(int n=0;n<SizeN;n++){
        for(int m=0;m<SizeM;m++){
            TransposedMatrix[m*SizeN+n]=MatrixToTranspose[n*SizeM+m];
        }
    }
    return TransposedMatrix;
}

float * TreeEM::ProductMatrix(float * MatrixL, float * MatrixR,int * Size){
    if(Size[1]!=Size[2]){
        cout<<"Size incompatibilities in the matrices"<<endl;
        return NULL;
    }
    float * MatrixProduct=new float[Size[0]*Size[3]];
    for(int i=0;i<Size[0]*Size[3];i++){
        MatrixProduct[i]=0;
    }
    for(int i=0;i<Size[0];i++){
        for(int k=0;k<Size[3];k++){
            for(int j=0;j<Size[1];j++){
                MatrixProduct[i+k*Size[0]]+=MatrixL[i+j*Size[0]]*MatrixR[j+k*Size[1]];
            }
        }
    }
    return MatrixProduct;
}

void TreeEM::invertMatrix(float * MatrixToInvert, int SizeMatrix)  {
//  if (size!=size2){
//      cout << "Matrix in not square" << endl;
//      return;
//    }
  if (SizeMatrix <= 1) {
      if(MatrixToInvert[0]!=0){
      MatrixToInvert[0]=1.0/MatrixToInvert[0];
      }
      return;
  }
  for (int i=1; i < SizeMatrix; i++) {
      MatrixToInvert[i] /= MatrixToInvert[0]; // normalize row 0
    }
  for (int i=1; i < SizeMatrix; i++)  {
      for (int j=i; j < SizeMatrix; j++)  { // do a column of L
          PrecisionTYPE sum = 0.0;
          for (int k = 0; k < i; k++){
              sum += MatrixToInvert[j*SizeMatrix+k] * MatrixToInvert[k*SizeMatrix+i];
            }
          MatrixToInvert[j*SizeMatrix+i] -= sum;
        }
      if (i == SizeMatrix-1) continue;
      for (int j=i+1; j < SizeMatrix; j++)  {  // do a row of U
          PrecisionTYPE sum = 0.0;
          for (int k = 0; k < i; k++){
              sum +=(PrecisionTYPE) MatrixToInvert[i*SizeMatrix+k]*MatrixToInvert[k*SizeMatrix+j];
            }
          MatrixToInvert[i*SizeMatrix+j] =(MatrixToInvert[i*SizeMatrix+j]-sum) / MatrixToInvert[i*SizeMatrix+i];
        }
    }
  for ( int i = 0; i < SizeMatrix; i++ )  {// invert L
    for ( int j = i; j < SizeMatrix; j++ )  {
        float x = 1.0;
        if ( i != j ) {
            x = 0.0;
            for ( int k = i; k < j; k++ ){
                x -= MatrixToInvert[j*SizeMatrix+k]*MatrixToInvert[k*SizeMatrix+i];
              }
          }
        MatrixToInvert[j*SizeMatrix+i] = x / MatrixToInvert[j*SizeMatrix+j];
      }
  }
  for ( int i = 0; i < SizeMatrix; i++ ) {  // invert U
    for ( int j = i; j < SizeMatrix; j++ )  {
        if ( i == j ){continue;}
        PrecisionTYPE sum = 0.0;
        for ( int k = i; k < j; k++ ){
            sum += (PrecisionTYPE)MatrixToInvert[k*SizeMatrix+j]*( (i==k) ? 1.0 : MatrixToInvert[i*SizeMatrix+k] );
          }
        MatrixToInvert[i*SizeMatrix+j] = -sum;
      }
  }
  for ( int i = 0; i < SizeMatrix; i++ ){   // final inversion
    for ( int j = 0; j < SizeMatrix; j++ )  {
        PrecisionTYPE sum = 0.0;
        for ( int k = ((i>j)?i:j); k < SizeMatrix; k++ ){
            sum +=(PrecisionTYPE) ((j==k)?1.0:MatrixToInvert[j*SizeMatrix+k])*MatrixToInvert[k*SizeMatrix+i];
          }
        MatrixToInvert[j*SizeMatrix+i] = (float)sum;
      }
  }
  return;
}

float TreeEM::determinant(float * data, int size){

  int i,j,j1,j2;
  float det = 0;
  float **m = NULL;
//  cout<<"Size is "<<size<<endl;
  if (size < 1) {
return 0;
    } else if (size == 1) { /* Shouldn't get used */
      det = data[0];
    } else if (size == 2) {
      det = data[0+size*0] * data[1+size*1] - data[1+size*0] * data[0+size*1];
    } else {
      det = 0;
      for (j1=0;j1<size;j1++) {
          m = (float **) malloc((size-1)*sizeof(float *));
          for (i=0;i<size-1;i++)
            m[i] = (float *) malloc((size-1)*sizeof(float));
          for (i=1;i<size;i++) {
              j2 = 0;
              for (j=0;j<size;j++) {
                  if (j == j1)
                    continue;
                  m[i-1][j2] = data[i+size*j];
                  j2++;
                }
            }
          det += (PrecisionTYPE)pow(-1.0,1.0+j1+1.0) * data[0+size*j1] * Determinant_lib(m,size-1);
          for (i=0;i<size-1;i++)
            free(m[i]);
          free(m);
        }
    }
  return (float)det;
}

float TreeEM::Determinant_lib(float **a,int n)
{
  int i,j,j1,j2;
  float det = 0;
  float **m = NULL;

  if (n < 1) { /* Error */

    } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
    } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
          m = (float **) malloc((n-1)*sizeof(float *));
          for (i=0;i<n-1;i++)
            m[i] = (float *) malloc((n-1)*sizeof(float));
          for (i=1;i<n;i++) {
              j2 = 0;
              for (j=0;j<n;j++) {
                  if (j == j1)
                    continue;
                  m[i-1][j2] = a[i][j];
                  j2++;
                }
            }
          det += (PrecisionTYPE)pow(-1.0,1.0+j1+1.0) * a[0][j1] * Determinant_lib(m,n-1);
          for (i=0;i<n-1;i++)
            free(m[i]);
          free(m);
        }
    }
  return((float)det);
}



/************* COPY METHODS ***********************/

TreeEM * TreeEM::CopyTree(TreeEM * ParentNew){
    TreeEM * CopiedTree=new TreeEM();

    // Copy all elements in the node according to preestablished copying methods
    CopiedTree->Parent=ParentNew;
    if (this->GetDataDirect()==NULL && ParentNew==NULL) {
        cout<<"Improper beginning of copy"<<endl;
        CopiedTree->DataImage=this->FindRoot()->CopyImage();
        CopiedTree->Mask=this->FindRoot()->CopyMask();
    }
    else{
        //cout<<"Proper beginning of copy"<<endl;
        CopiedTree->DataImage=this->CopyImage();
        CopiedTree->Mask=this->CopyMask();
    }
    CopiedTree->Priors=this->CopyPriors();
    CopiedTree->PriorsAdapted=this->CopyPriorsAdapted();
    CopiedTree->L2S=this->CopyL2S();
    CopiedTree->S2L=this->CopyS2L();
    //    CopiedTree->Distribution=this->CopyDistribution();
    //    CopiedTree->NonNormResp=this->CopyNonNormResp();
    CopiedTree->NormResp=this->CopyNormResp();
    CopiedTree->NonNormWeight=this->GetNonNormWeight();
    CopiedTree->NormWeight=this->GetNormWeight();
    CopiedTree->ParametersDistribution=this->CopyParameters();
    CopiedTree->SplitCheck=this->CopySplitCheck();
    CopiedTree->MergeCheck=this->CopyMergeCheck();
    //    CopiedTree->BFCorrection=this->CopyBFCorrection();
    CopiedTree->DataBFCorrected=this->CopyDataBFCorrected();
    CopiedTree->BFCoeffs=this->CopyBFCoeffs();
    CopiedTree->IndFactor=this->CopyIndFactor();
    CopiedTree->DPChildren=this->CopyDPChildren();
    CopiedTree->FlagDistClassInd=this->FlagDistClassInd;
    CopiedTree->FlagOutliers=this->FlagOutliers;
    CopiedTree->NumberMaskedElements=this->NumberMaskedElements;
    CopiedTree->HardSeg=this->CopyHardSeg();
    CopiedTree->MRF=this->CopyMRF();

    // Copying for all children of this but with the proper parent pointer
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        CopiedTree->Children.push_back((this->GetChild(c))->CopyTree(CopiedTree));
    }
    //cout<<"Node copied"<<endl;
    return CopiedTree;
}

// Returns a pointer to a nifti_image that is a strict copy from the data image directly obtained in this
nifti_image * TreeEM::CopyImage(){
    if (this->GetDataDirect()==NULL) {// Meaning that it is not the root
        //cout<<"Direct Data is NULL"<<endl;
        return NULL;
    }

    else{ // DataImage is directly available so we have to copy it;
        if(this->GetDataImage()->datatype!=DT_FLOAT){
            this->ConvertDataImageToFloat();
        }
        //cout<<this->GetDataImage()<<" and "<<this->GetDataDirect()<<endl;
        nifti_image * CopiedImage=nifti_copy_nim_info(this->GetDataImage());
        // CopiedImage is at the moment with a NULL pointer to the data. According memory must then be allocated and filled with the values
        int numbvox=this->GetDataImage()->nvox;
        CopiedImage->data=(void *)calloc(this->GetDataImage()->nvox,sizeof(float));
        float * CopiedImageData=static_cast<float *>(CopiedImage->data);
        float * CopiedImageData_PTR=CopiedImageData;
        float * Data_PTR=static_cast<float *>(this->GetDataImage()->data);
        for (int i=0; i<numbvox; i++,Data_PTR++,CopiedImageData_PTR++) {
            *CopiedImageData_PTR=*Data_PTR;
        }
        //cout<<"Image copied"<<endl;
        return CopiedImage;
    }
}

// Returns as a nifti_image binary converted the mask used if a mask is used. Similar to CopyImage since the Mask can only be at the root also
nifti_image * TreeEM::CopyMask(){
    if (this->GetMaskDirect()==NULL) {
        return NULL;
    }
    else{
        if(this->GetMask()->datatype!=DT_BINARY){
            this->MakeMaskBinary();
        }
        nifti_image * CopiedMask=nifti_copy_nim_info(this->GetMask());
        int numbvox=this->GetMask()->nvox;
        CopiedMask->data=(void *)calloc(this->GetMask()->nvox,sizeof(bool));
        bool * CopiedMaskData=static_cast<bool *>(CopiedMask->data);
        bool * CopiedMaskData_PTR=CopiedMaskData;
        bool * Mask_PTR=static_cast<bool *>(this->GetMask()->data);
        for (int i=0; i<numbvox; i++,Mask_PTR++,CopiedMaskData_PTR++) {
            *CopiedMaskData_PTR=*Mask_PTR;
        }
        //cout<<" Mask copied "<<endl;
        return CopiedMask;
    }
}

// Returns a pointer to a nifti_image float converted that is the copy of the priors used.
nifti_image * TreeEM::CopyPriors(){
    if (this->GetPriorsDirect()==NULL) {
        return NULL;
    }
    else{
        if (this->GetPriorsDirect()->datatype!=DT_FLOAT) {
            this->MakePriorsFloat();
        }
        nifti_image * CopiedPriors=nifti_copy_nim_info(this->GetPriorsDirect());
        int numbvox=this->GetPriors()->nvox;
        CopiedPriors->data=(void *)calloc(this->GetPriors()->nvox,sizeof(float));
        float * CopiedPriorsData=static_cast<float *>(CopiedPriors->data);
        float * CopiedPriorsData_PTR=CopiedPriorsData;
        float * Priors_PTR=static_cast<float *>(this->GetPriors()->data);
        for (int i=0; i<numbvox; i++,Priors_PTR++,CopiedPriorsData_PTR++) {
            *CopiedPriorsData_PTR=*Priors_PTR;
        }
        //cout<<" Priors copied"<< endl;
        return CopiedPriors;
    }
}

// Only stored at general classes. If exist makes copy of PriorsAdapted and returns corresponding float pointer
float * TreeEM::CopyPriorsAdapted(){
    if (this->GetPriorsAdaptedDirect()==NULL) {
        return NULL;
    }
    else {
        int numel=this->GetNumberElements();
        float * CopiedPriorsAdapted=new float[numel];
        float * CopiedPriorsAdapted_PTR=CopiedPriorsAdapted;
        float * PriorsAdapted_PTR=this->GetPriorsAdapted();
        for (int i=0; i<numel; i++,CopiedPriorsAdapted_PTR ++,PriorsAdapted_PTR++) {
            *CopiedPriorsAdapted_PTR=*PriorsAdapted_PTR;
        }
        return CopiedPriorsAdapted;
    }
}

// Similar to previous ones since it is only stored at the Root. If exists makes the copy of L2S and returns a int pointer
int * TreeEM::CopyL2S(){
    if (this->GetL2SDirect()==NULL) {
        return NULL;
    }
    else {
        int numel=this->GetNumberElements();
        int * CopiedL2S=new int[numel];
        int * CopiedL2S_PTR=CopiedL2S;
        int * L2S_PTR=this->GetL2S();
        for (int i=0; i<numel; i++,CopiedL2S_PTR ++,L2S_PTR++) {
            *CopiedL2S_PTR=*L2S_PTR;
        }
        return CopiedL2S;
    }
}

int * TreeEM::CopyS2L(){
    if (this->GetS2LDirect()==NULL) {
        return NULL;
    }
    else {
        int numelmasked=this->GetNumberMaskedElements();
        int * CopiedS2L=new int[numelmasked];
        int * CopiedS2L_PTR=CopiedS2L;
        int * S2L_PTR=this->GetS2L();
        for (int i=0; i<numelmasked; i++,CopiedS2L_PTR++,S2L_PTR++) {
            *CopiedS2L_PTR=*S2L_PTR;
        }
        return CopiedS2L;
    }
}

int * TreeEM::CopyHardSeg(){
    if(this->GetHardSegDirect()==NULL){
        return NULL;
    }
    else{
        int numelmasked=this->GetNumberMaskedElements();
        int * CopiedHardSeg=new int[numelmasked];
        int * CopiedHardSeg_PTR=CopiedHardSeg;
        int * HardSeg_PTR=this->GetHardSeg();
        for (int i=0; i<numelmasked; i++,CopiedHardSeg_PTR++,HardSeg_PTR++) {
            *CopiedHardSeg_PTR=*HardSeg_PTR;
        }
        return CopiedHardSeg;
    }
}

//float * TreeEM::CopyDistribution(){
//    if (this->GetS2L()==NULL) {
//        return NULL;
//    }
//    else {
//        int numelmasked=this->GetNumberMaskedElements();
//        float * CopiedDistribution=new float[numelmasked];
//        float * CopiedDistribution_PTR=CopiedDistribution;
//        float * Distribution_PTR=this->GetDistribution();
//        for (int i=0; i<numelmasked; i++,CopiedDistribution_PTR++,Distribution_PTR++) {
//            *CopiedDistribution_PTR=*Distribution_PTR;
//        }
//        return CopiedDistribution;
//    }
//}

//float * TreeEM::CopyNonNormResp(){
//    if (this->GetS2L()==NULL) {
//        return NULL;
//    }
//    else {
//        int numelmasked=this->GetNumberMaskedElements();
//        float * CopiedNonNormResp=new float[numelmasked];
//        float * CopiedNonNormResp_PTR=CopiedNonNormResp;
//        float * NonNormResp_PTR=this->GetNonNormResp();
//        for (int i=0; i<numelmasked; i++,CopiedNonNormResp_PTR++,NonNormResp_PTR++) {
//            *CopiedNonNormResp_PTR=*NonNormResp_PTR;
//        }
//        return CopiedNonNormResp;
//    }
//}

float * TreeEM::CopyNormResp(){
    if (this->GetS2L()==NULL) {
        return NULL;
    }
    else {
        int numelmasked=this->GetNumberMaskedElements();
        float * CopiedNormResp=new float[numelmasked];
        float * CopiedNormResp_PTR=CopiedNormResp;
        float * NormResp_PTR=this->GetNormResp();
        for (int i=0; i<numelmasked; i++,CopiedNormResp_PTR++,NormResp_PTR++) {
            *CopiedNormResp_PTR=*NormResp_PTR;
        }
        return CopiedNormResp;
    }
}

float TreeEM::CopyIndFactor(){
    return this->IndFactor;
}

Parameters * TreeEM::CopyParameters(){
    if (this->GetParameters()==NULL) {
        return NULL;
    }
    else{
        Parameters * CopiedParameters=new Parameters;
        CopiedParameters->DistributionType=this->GetDistributionType();
        CopiedParameters->SizeParameters=this->GetSizeParameters();
        int sizeParams=CopiedParameters->SizeParameters;
        if (sizeParams==0) { // Distinction done if the size of the Parameters is more than 0 or not to allocate proper amount of memory
            CopiedParameters->ValueParameters=NULL;
        }
        else{
            CopiedParameters->ValueParameters=new float[CopiedParameters->SizeParameters];
            float * CopiedValueParameters_PTR=CopiedParameters->ValueParameters;

            float * ValueParameters_PTR=this->GetParametersValue();
            for (int i=0; i<sizeParams; i++,CopiedValueParameters_PTR++,ValueParameters_PTR++) {
                *CopiedValueParameters_PTR=*ValueParameters_PTR;
            }
        }
        return CopiedParameters;
    }
}

bool * TreeEM::CopyMergeCheck(){
    if (this->GetMergeCheck()==NULL) {
        return NULL;
    }
    else{
        int numbDirectLeaves=this->GetNumberDirectLeaves();
        int numbchild=this->GetNumberChildren();
        if (numbDirectLeaves==0) {
            return NULL;
        }
        bool * MergeCheck_PTR=this->GetMergeCheck();
//        bool * CopiedMergeCheck=new bool[numbDirectLeaves*numbDirectLeaves];
        bool * CopiedMergeCheck=new bool[numbchild * numbchild];
//        for(int i=0;i<numbDirectLeaves*numbDirectLeaves;i++){
//            CopiedMergeCheck[i]=false;
//        }
        for(int i=0;i<numbchild*numbchild;i++){
            CopiedMergeCheck[i]=false;
        }
        bool * CopiedMergeCheck_PTR=CopiedMergeCheck;
//        for (int i=0; i<numbDirectLeaves*numbDirectLeaves; i++,CopiedMergeCheck_PTR++,MergeCheck_PTR++) {
//            *CopiedMergeCheck_PTR=*MergeCheck_PTR;
//        }
        for (int i=0; i<numbchild*numbchild; i++,CopiedMergeCheck_PTR++,MergeCheck_PTR++) {
            *CopiedMergeCheck_PTR=*MergeCheck_PTR;
        }
        return CopiedMergeCheck;
    }
}

bool * TreeEM::CopySplitCheck(){
    if (this->GetSplitCheck()==NULL) {
        return NULL;
    }
    else{
        int numbDirectLeaves=this->GetNumberDirectLeaves();
        int numbchild=this->GetNumberChildren();
        if (numbDirectLeaves==0) {
            return NULL;
        }
        bool * SplitCheck_PTR=this->GetSplitCheck();
        bool * CopiedSplitCheck=new bool[numbchild];
        for(int i=0;i<numbchild;i++){
            CopiedSplitCheck[i]=false;
        }
        bool * CopiedSplitCheck_PTR=CopiedSplitCheck;
        for (int i=0; i<numbchild; i++,CopiedSplitCheck_PTR++,SplitCheck_PTR++) {
            *CopiedSplitCheck_PTR=*SplitCheck_PTR;
        }
        return CopiedSplitCheck;
    }
}

//// Returns a float pointer to an array where the BFBasisFunctions have been copied
//float * TreeEM::CopyBFCorrection(){
//    // Normally the BF are at the root and only stored there, elsewhere it is a NULL pointer
//    if (this->GetBFCorrectionDirect()==NULL) {
//        return NULL;
//    }
//    else{
//        int numelmasked=this->GetNumberMaskedElements();
//        int numbmodal=this->GetNumberModalities();
//        //int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
//        float * CopiedBFCorrection=new float[numbmodal*numelmasked];//{0};
//        for (int i=0; i<numbmodal*numelmasked; i++) {
//            CopiedBFCorrection[i]=0;
//        }
//        float * CopiedBFCorrection_PTR=CopiedBFCorrection;
//        float * BFCorrection_PTR=this->GetBFCorrectionDirect();
//        for (int i=0; i<numbmodal*numelmasked; i++,CopiedBFCorrection_PTR++,BFCorrection_PTR++) {
//            *CopiedBFCorrection_PTR=*BFCorrection_PTR;
//        }
//        return CopiedBFCorrection;
//    }
//}

float * TreeEM::CopyDataBFCorrected(){
    // Normally the BF are at the root and only stored there, elsewhere it is a NULL pointer
    if (this->GetDataBFCorrectedDirect()==NULL) {
        return NULL;
    }
    else{
        int numelmasked=this->GetNumberMaskedElements();
        int numbmodal=this->GetNumberModalities();
        //int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
        float * CopiedDataBFCorrected=new float[numbmodal*numelmasked];//{0};
        for (int i=0; i<numbmodal*numelmasked; i++) {
            CopiedDataBFCorrected[i]=0;
        }
        float * CopiedDataBFCorrected_PTR=CopiedDataBFCorrected;
        float * DataBFCorrected_PTR=this->GetDataBFCorrected();
        for (int i=0; i<numbmodal*numelmasked; i++,CopiedDataBFCorrected_PTR++,DataBFCorrected_PTR++) {
            *CopiedDataBFCorrected_PTR=*DataBFCorrected_PTR;
        }
        return CopiedDataBFCorrected;
    }
}

// Returns a vector of pointers to the arrays containing the coefficients corresponding to the different modalities
float * TreeEM::CopyBFCoeffs(){
    //Declaration of the element to return;
    int numbmodal=this->GetNumberModalities();
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    float * BFCoeffsToCopy=this->GetBFCoeffsDirect();


    // 1st case where there is no available coeffs, the size of the vector to copy is 0;
    if (this->GetBFCoeffsDirect()==NULL) {
        return NULL;
    }

    else{
        float * CopiedBFCoeffs=new float[numbmodal*numbBF];
        for(int i=0;i<numbmodal*numbBF;i++){
            CopiedBFCoeffs[i]=0;
        }
        for(int l=0;l<numbBF*numbmodal;l++){
            CopiedBFCoeffs[l]=BFCoeffsToCopy[l];
        }
        return CopiedBFCoeffs;
    }
}

float* TreeEM::CopyDPChildren(){
    float* DPChildrenToCopy=this->GetDPChildrenDirect();
    float* CopiedDPChildren=NULL;
//    int DPsize=DPChildrenToCopy.size();
    if(DPChildrenToCopy==NULL){
        return CopiedDPChildren;
    }
    else{
        int numbchild=this->GetNumberChildren();
        if(this->GetFlagDistClassInd()){

            CopiedDPChildren=new float[MaxSupport*numbchild];
            for(int i=0;i<MaxSupport*numbchild;i++){
                CopiedDPChildren[i]=DPChildrenToCopy[i];
            }
        }
        else{
            int SizeHistogram=(int)pow_int(MaxSupport,numbchild);
            CopiedDPChildren=new float[SizeHistogram];
            for(int i=0;i<SizeHistogram;i++){
                CopiedDPChildren[i]=DPChildrenToCopy[i];
            }
        }
        return CopiedDPChildren;
    }
}

TreeEM::~TreeEM (){
    MakeEmpty ();
}

void TreeEM::MakeEmpty (){
    int numbchild=this->GetNumberChildren();
    for (int c = 0; c < numbchild; c++)
    {
        if(this->GetChild(c)!=NULL){
            delete this->Children [c];
            this->Children[c] = NULL;
        }
    }
    this->Children.clear();
    if (this->GetBFCoeffsDirect()!=NULL) {
        delete[] this->BFCoeffs;
        this->BFCoeffs=NULL;
    }
    //    if (this->GetBFCorrectionDirect()!=NULL) {
    //        delete[] this->BFCorrection;
    //        this->BFCorrection=NULL;
    //    }
    if(this->GetDataBFCorrectedDirect()!=NULL){
        delete [] this->DataBFCorrected;
        this->DataBFCorrected=NULL;
    }
    if (this->GetDataDirect()!=NULL) {
        nifti_image_free(this->GetDataDirect());
        this->DataImage=NULL;
    }
    if (this->GetMaskDirect()!=NULL) {
        nifti_image_free(this->GetMaskDirect());
        this->Mask=NULL;
    }
    if (this->GetL2SDirect()!=NULL) {
        delete [] this->GetL2SDirect();
        this->L2S=NULL;
    }
    if (this->GetS2LDirect()!=NULL) {
        delete [] this->GetS2LDirect();
        this->S2L=NULL;
    }
    if(this->GetHardSegDirect()!=NULL){
        delete [] this->GetHardSegDirect();
        this->HardSeg=NULL;
    }
    if(this->GetMRF()!=NULL){
        delete [] this->GetMRF();
        this->MRF=NULL;
    }
    if(this->GetDPChildrenDirect()!=NULL){
        delete[] this->GetDPChildrenDirect();
        this->DPChildren=NULL;
    }
    if (this->GetPriorsDirect()!=NULL) {
        nifti_image_free(this->GetPriorsDirect());
        this->Priors=NULL;
    }
    if(this->GetPriorsAdaptedDirect()!=NULL){
        delete [] PriorsAdapted;
        this->PriorsAdapted=NULL;
    }
    //    if (this->NonNormResp!=NULL){
    //        delete this->NonNormResp;
    //        this->NonNormResp=NULL;
    //    }
    if (this->NormResp!=NULL){
        delete [] this->NormResp;
        this->NormResp=NULL;
    }
    //    if (this->Distribution!=NULL){
    //        delete []this->Distribution;
    //        this->Distribution=NULL;
    //    }
    if (this->ParametersDistribution!=NULL){
        delete this->ParametersDistribution;
        this->ParametersDistribution=NULL;
    }
    if (this->GetSplitCheck()!=NULL) {
        delete [] this->GetSplitCheck();
        this->SplitCheck=NULL;
    }
    if (this->GetMergeCheck()!=NULL) {
        delete [] this->GetMergeCheck();
        this->MergeCheck=NULL;
    }

    this->Parent=NULL;
    return;
}

// Get Functions
TreeEM * TreeEM::GetParent(){
    return this->Parent;
}

vector<TreeEM *> TreeEM::GetChildren(){
    return this->Children;
}

// Returns the number of children
int TreeEM::GetNumberChildren(){
    return this->GetChildren().size();
}

// Returns the number of anatomical classes : rules for anatomical classes, if without outlier class, number of children of the root if with outliers, supposed to be at second level under first child
int TreeEM::GetNumberGeneralClasses(){
    int OutlierType = this->GetFlagOutliers();
    switch(OutlierType){
    case 0:{ // case no outlier model
      return this->GetNumberChildren();
      break;
    }
    case 1:{ // case outlier model in three level
      return this->GetChild(0)->GetNumberChildren();
      break;
    }
    case 2:{ // case outlier model in horizontal model
      return this->GetNumberChildren()-1;
      break;
    }
    default :{ // no outlier model;
        return this->GetNumberChildren();
        break;
    }
    }
}

// Gets the child of index c after checking that we have at least c+1 children
TreeEM * TreeEM::GetChild(int c){
    if (c>=this->GetNumberChildren()) {
        cout<<"Getting over bounds of number of children"<<endl;
        cout<<"c and the address of the parent are "<<c<<" "<<this<<endl;
        return NULL;
    }
    else {
        return this->GetChildren()[c];
    }
}

nifti_image * TreeEM::GetDataDirect(){
    return this->DataImage;
}

nifti_image * TreeEM::GetMaskDirect(){
    return this->Mask;
}

nifti_image * TreeEM::GetMask(){
    if (this->GetMaskDirect()!=NULL) {
        //cout<<"Has direct Mask reference"<<endl;
        return this->GetMaskDirect();
    }
    else if (this->GetParent()==NULL){ // Means it is a root or an initial node
        //cout<<"Root or IN"<<endl;
        return NULL;
    }
    else {
        return this->GetParent()->GetMask();
    }
}

nifti_image * TreeEM::GetPriors(){
    if (this->Priors !=NULL){
        return this->Priors;
    }
    else if(this->GetParent()==NULL){// Means it is a root or an initial node
        return NULL;
    }
    else {
        return this->GetParent()->GetPriors();
    }
}

float * TreeEM::GetPriorsAdapted(){
    if(this->PriorsAdapted!=NULL){
        return this->PriorsAdapted;
    }
    else if(this->GetParent()==NULL){// Means it is a root or an initial node
        return NULL;
    }
    else {
        return this->GetParent()->GetPriorsAdapted();
    }
}

nifti_image * TreeEM::GetPriorsDirect(){
    return this->Priors;
}

float * TreeEM::GetPriorsAdaptedDirect(){
    return this->PriorsAdapted;
}

//float * TreeEM::GetNonNormResp(){
//    return this->NonNormResp;
//}

float * TreeEM::GetNormResp(){
    return this->NormResp;
}

float TreeEM::GetNonNormWeight(){
    return this->NonNormWeight;
}

float TreeEM::GetNormWeight(){
    return this->NormWeight;
}

float TreeEM::GetPartNormWeight(){
    if (this->Parent==NULL || this->GetParent()->GetPriorsDirect()!=NULL) {
        return this->GetNormWeight();
    }
    else return this->GetNormWeight()*(this->GetParent())->GetPartNormWeight();
}

float TreeEM::GetPartNormWeightAbovePriors(){
    TreeEM* NodePrior=this->FindGeneralClassPriors();
    if(NodePrior==NULL){
        return 1;
    }
    else return NodePrior->GetParent()->GetPartNormWeight();
}

Parameters * TreeEM::GetParameters(){
    return this->ParametersDistribution;
}

int TreeEM::GetDistributionType(){
    return this->GetParameters()->DistributionType;
}

//float * TreeEM::GetDistribution(){
//    return this->Distribution;
//}

int * TreeEM::GetL2SDirect(){
    return this->L2S;
}

int * TreeEM::GetS2LDirect(){
    return this->S2L;
}

int * TreeEM::GetHardSegDirect(){
    return this->HardSeg;
}

bool * TreeEM::GetSplitCheck(){
    return this->SplitCheck;
}

bool * TreeEM::GetMergeCheck(){
    return this->MergeCheck;
}

/*Returns the number of bins used for the histogram. Has to be the same for each member. Therefore only stored at the root and put to 0 everywhere else*/
int TreeEM::GetNumbbins(){
    return numbbins;
}

// Returns the pointer to the beginning of the values for the means in the parameters values
/* WARNING : do not create a copy of the mean values in another array*/
float * TreeEM::GetMean(){
    if (this->GetNumberChildren()!=0) {
        return NULL;
    }
    else {
        return this->GetParameters()->ValueParameters;
    }
}

// Creates the PrecisionTYPE array and returns the pointer to it which will hold the mean for the considered NormResp given the values (no need for the considered tree to be a leaf)
float * TreeEM::GetMeanDirect(){
    int numbmodal=this->GetNumberModalities();
    int numel=this->GetNumberElements();
    float * MeanResult=new float[numbmodal];//{0};
    for (int i=0; i<numbmodal; i++) {
        MeanResult[i]=0;
    }
    float * Data=static_cast<float*>(this->GetDataImage()->data);
    float * Data_PTR=Data;
    int * L2S_PTR=this->GetL2S();

    for (int m=0; m<numbmodal; m++) {
        Data_PTR=&Data[m*numel];
        PrecisionTYPE MeanResult_tmp=0;
        L2S_PTR=this->GetL2S();
        float * NormResp_PTR=this->GetNormResp();
        PrecisionTYPE SumNormResp=0;
        for (int i=0; i<numel; i++,Data_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
                MeanResult_tmp+=(PrecisionTYPE)*Data_PTR*(*NormResp_PTR);
                SumNormResp+=(PrecisionTYPE)*NormResp_PTR;
                NormResp_PTR++;
            }
        }
        MeanResult[m]=(float)MeanResult_tmp/SumNormResp;
    }
    return MeanResult;
}

float * TreeEM::GetVarianceDirect(){
    //    int numel=this->GetNumberElements();
        int numbmodal=this->GetNumberModalities();
        int numelmasked=this->GetNumberMaskedElements();
        PrecisionTYPE sumResp=0; // Initialisation of the denominator
        //float * PointerToDataBegin = static_cast<float *>(this->GetDataImage()->data);
        //    float * PointerToDataBegin=this->MakeDataBFCorrected();
        float * PointerToDataBegin=this->GetDataBFCorrected();
        float * PointerToDataBegin_PTR1=PointerToDataBegin;
        float * PointerToDataBegin_PTR2=PointerToDataBegin;
        float * NormalisedResponsabilities_PTR=this->GetNormResp();
        //    int * L2S_PTR=this->GetL2S();

        // Calculation of the denominator (sum over the active voxels of the normalised responsabilities)
        for (int i=0; i<numelmasked; i++,NormalisedResponsabilities_PTR++) {
            sumResp+=(PrecisionTYPE)(*NormalisedResponsabilities_PTR);
        }

        float * VarianceDirect=new float[numbmodal*numbmodal];
        float * MeanDirect=this->GetMeanDirect();
        for(int m1=0;m1<numbmodal;m1++){
            for(int m2=0;m2<numbmodal;m2++){
                    VarianceDirect[m1+m2*numbmodal]=0;
            }
        }

        for(int m1=0;m1<numbmodal;m1++){
            // First data pointer to the beginning of the modality m1 considered
            for(int m2=0;m2<numbmodal;m2++){
                PointerToDataBegin_PTR1=&PointerToDataBegin[m1*numelmasked];
                PointerToDataBegin_PTR2=&PointerToDataBegin[m2*numelmasked]; // Second Data pointer to the beginning of the modality m2 considered
                NormalisedResponsabilities_PTR=this->GetNormResp(); // Reinitialisation of the responsabilities pointer to the beginning
                PrecisionTYPE VarianceToUpdate_tmp=0;
                for(int i=0;i<numelmasked;i++,PointerToDataBegin_PTR1++,PointerToDataBegin_PTR2++,NormalisedResponsabilities_PTR++){
                    // Update of the numerator of the Variance Calculation only if in the case of an active voxel
                    VarianceToUpdate_tmp+=(PrecisionTYPE)(*NormalisedResponsabilities_PTR)*((*PointerToDataBegin_PTR1)-MeanDirect[m1])*((*PointerToDataBegin_PTR2)-MeanDirect[m2]);
                }
                if (sumResp !=0) {
                    VarianceDirect[m1+m2*numbmodal]=(float)VarianceToUpdate_tmp/sumResp;
                    if (m1==m2) {
                        VarianceDirect[m1+m2*numbmodal]=VarianceDirect[m1+m2*numbmodal]<=1E-6?1E-6:VarianceDirect[m1+m2*numbmodal]; // in order to avoid going to 0 if too sharp distribution but not changing non diagonal of variance
                    }


                }
                else{
                    VarianceDirect[m1+m2*numbmodal]=VarianceToUpdate_tmp/numelmasked;// this->GetNumberMaskedElements();
                }
                // Use of the symmetry property of the Variance matrix
                VarianceDirect[m2+m1*numbmodal]=VarianceDirect[m1+m2*numbmodal];
            }
        }
        // Clearing memory
        delete [] MeanDirect;
        MeanDirect=NULL;
        return VarianceDirect;
}

float * TreeEM::GetDiagVarianceDirect(){
    int numbmodal=this->GetNumberModalities();
    int numel=this->GetNumberElements();
    float * MeanUsedB=this->GetMeanDirect();
    float MeanUsed[MaxNumbModal];
    for(int m=0;m<MaxNumbModal;m++){
        if(m<numbmodal){
            MeanUsed[m]=MeanUsedB[m];
        }
        else{
            MeanUsed[m]=0;
        }
    }
    delete [] MeanUsedB;
    float * Data=static_cast<float *>(this->GetDataImage()->data);
    float * Data_PTR=Data;
    int * L2S_PTR=this->GetL2S();
    float * NormResp_PTR=this->GetNormResp();
    float * DiagVarianceResult=new float[numbmodal];//{0};
    for (int i=0; i<numbmodal; i++) {
        DiagVarianceResult[i]=0;
    }
    PrecisionTYPE SumNormResp=0;
    PrecisionTYPE DiagVarianceResult_tmp=0;
    for (int m=0; m<numbmodal; m++) {
        Data_PTR=&Data[m*numel];
        L2S_PTR=this->GetL2S();
        NormResp_PTR=this->GetNormResp();
        SumNormResp=0;
        DiagVarianceResult_tmp=0;
        for (int i=0; i<numel; i++,Data_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
                DiagVarianceResult_tmp+=(PrecisionTYPE)*NormResp_PTR*pow_int((*Data_PTR-MeanUsed[m]), 2);
                SumNormResp+=(PrecisionTYPE)*NormResp_PTR;
                NormResp_PTR++;
            }
        }
        DiagVarianceResult[m]=(float)DiagVarianceResult_tmp/SumNormResp;
    }

    return DiagVarianceResult;
}

// Returns the pointer to the beginning of the parameters values concerned by the variance
/*WARNING : do not create a copy of the values pointed as belonging to the variance matrix*/
float * TreeEM::GetVariance(){
    if (this->GetNumberChildren()!=0) { /* If the considered element is not a leaf, there is no parameters since the distribution is of mixture type*/
        return NULL;
    }
    else {
        switch (this->GetDistributionType()) { /* For the moment only Gaussian distribution = default case might change in the future*/

        default:
            return this->GetParameters()->ValueParameters+this->GetNumberModalities();
            break;
        }
    }
}


int TreeEM::GetSizeParameters(){
    if (this->GetDistributionType()==1 && this->GetParameters()->SizeParameters==(int)((this->GetNumberModalities()+1)*this->GetNumberModalities())) {
        return this->GetParameters()->SizeParameters;
    }
    //cout<<"Not Gaussian or not appropriate size of parameters"<< this->GetDistributionType()<< endl;

    return this->GetParameters()->SizeParameters;
}

float * TreeEM::GetParametersValue(){
    return this->GetParameters()->ValueParameters;
}

//// Returns in a vector of pointers to nifti_images the Priors of the children
//vector<nifti_image *> TreeEM::GetPriorsVector(){
//    vector<nifti_image *> PriorsVector;
//    int numbchild=this->GetNumberChildren();
//    for(int c=0;c<numbchild;c++){
//        PriorsVector.push_back(this->GetChild(c)->GetPriors());
//        //cout<<PriorsVector[c];
//    }
//    return PriorsVector;
//}

vector<nifti_image*> TreeEM::GetPriorsVector(){
    vector<nifti_image*> PriorsVector;
    if(this->GetPriorsDirect()!=NULL){
        PriorsVector.push_back(this->GetPriorsDirect());
    }
    else{
        int numbchild=this->GetNumberChildren();
        for(int c=0;c<numbchild;c++){
            vector<nifti_image*> PriorsVectorChild=this->GetChild(c)->GetPriorsVector();
            for(int i=0;i<PriorsVectorChild.size();i++){
                PriorsVector.push_back(PriorsVectorChild[i]);
            }
        }
    }
    return PriorsVector;
}

// Return the vector of nodes in which the statistical atlases are stored (may not be first level if outlier model is chosen)
vector<TreeEM*> TreeEM::GetPriorsNodeVector(){
    vector<TreeEM*> PriorNodeVector;
    if(this->GetPriorsDirect()!=NULL){
        PriorNodeVector.push_back(this);
    }
    else{
        int numbchild=this->GetNumberChildren();
        for(int c=0;c<numbchild;c++){
            vector<TreeEM*> PriorNodeVectorChild=this->GetChild(c)->GetPriorsNodeVector();
            for(int i=0;i<PriorNodeVectorChild.size();i++){
                PriorNodeVector.push_back(PriorNodeVectorChild[i]);
            }
        }
    }
    return PriorNodeVector;
}

vector<TreeEM*> TreeEM::GetGeneralClassesVector(){
    vector<TreeEM*> GeneralClassesVector;
    TreeEM * Root=this->FindRoot();
    switch(this->GetFlagOutliers()){
    case 0:{ // no outlier model considered
        GeneralClassesVector= Root->GetChildren();
        break;
    }
    case 1:{
        GeneralClassesVector=Root->GetChild(0)->GetChildren();
        break;
    }
    case 2:{
        int numbchild=Root->GetNumberChildren();
        for(int c=0;c<numbchild-1;c++){
            GeneralClassesVector.push_back(Root->GetChild(c));
        }
        break;
    }
    default:{
        GeneralClassesVector= Root->GetChildren();
        break;
    }
    }
    return GeneralClassesVector;
}

// Returns the vector of trees with the main general classes and the outlier node
vector<TreeEM*> TreeEM::GetMainNodesVector(){
    if(this->GetFlagOutliers()!=1){ // no outlier class or outlier class at same level with evolving priors
        return this->FindRoot()->GetChildren();
    }
    else{// priors at level 2 + Node outlier
        vector<TreeEM *> GeneralClassesVector=this->FindRoot()->GetChild(0)->GetChildren();
        GeneralClassesVector.push_back(this->GetNodeOutlier());
        return GeneralClassesVector;
    }
}

// Return the pointer to the node from which the outlier classes will be derived
TreeEM* TreeEM::GetNodeOutlier(){
    if(!this->GetFlagOutliers()){
        return NULL;
    }
    else{ // the outlier node is the last of the children from the root Rule
        TreeEM* Root=this->FindRoot();
        int numbchild=Root->GetNumberChildren();
        return Root->GetChild(numbchild-1);
    }
}

//// Return in a vector of pointers to float arrays the PriorsAdapted of the children
//vector<float *> TreeEM::GetPriorsAdaptedVector(){
//    vector<float *> PriorsAdaptedVector;
//    int numbchild=this->GetNumberChildren();
//    for(int c=0;c<numbchild;c++){
//        PriorsAdaptedVector.push_back(this->GetChild(c)->GetPriorsAdapted());
//        //cout<<PriorsVector[c];
//    }
//    return PriorsAdaptedVector;
//}

vector<float*> TreeEM::GetPriorsAdaptedVector(){
    vector<float*> PriorsAdaptedVector;
    vector<TreeEM*> NodePriorVector=this->GetPriorsNodeVector();
    int numbclasses=NodePriorVector.size();
    for(int c=0;c<numbclasses;c++){
        PriorsAdaptedVector.push_back(NodePriorVector[c]->GetPriorsAdapted());
    }
    return PriorsAdaptedVector;
}

nifti_image * TreeEM::GetDataImage(){
    //if (this->WhatTypeOfTreeComponent()==ROOT || this->WhatTypeOfTreeComponent()==INITIALNODE){
    if(this->GetParent()==NULL){
        //cout<< "The considered tree is a root or an initial node"<<endl;
        if( this->DataImage==NULL) {
            cout<<"This tree is not valid : no input image is given !"<<endl;

            return NULL;
        }
        else{
            return this->DataImage;
        }
    }
    else{
        return this->GetParent()->GetDataImage();
    }
}

// Returns the max of modality m in data image after checking it can be found
float TreeEM::GetMaxDataModal(int m){
    float MaxResult=-1E32;
    if (this->GetDataImage()==NULL) { // no data image attributed in the tree
        cout<< "No data attributed to the tree, tree not valid"<< endl;
        return MaxResult;
    }
    else{
        // Check if modality wanted is compatible with Data
        if (m>=this->GetNumberModalities()){// modality wanted is out of bounds
            cout<< "Modality wanted is out of bounds"<<endl;
            return MaxResult;
        }
        // Everything has been checked so we look for the maximum
        float * Data_PTR=static_cast<float *>(this->GetDataImage()->data);
        int numel=this->GetNumberElements();
        float * Data_PTRm=&Data_PTR[m*numel];

        for (int i=0; i<numel; i++,Data_PTRm++) {
            if((*Data_PTRm)>MaxResult){ // update of the maximum
                MaxResult=(*Data_PTRm);
            }
        }
        //cout << "the max is "<<MaxResult<<endl;
        return MaxResult;
    }
}

// Returns the min of modality m in data image after checking it can be found
float TreeEM::GetMinDataModal(int m){
    float MinResult=1E32;
    if (this->GetDataImage()==NULL) { // no data image attributed in the tree
        cout<< "No data attributed to the tree, tree not valid"<< endl;
        return MinResult;
    }
    else{
        // Check if modality wanted is compatible with Data
        if (m>=this->GetNumberModalities()){// modality wanted is out of bounds
            cout<< "Modality wanted is out of bounds"<<endl;
            return MinResult;
        }
        // Everything has been checked so we look for the maximum
        float * Data_PTR=static_cast<float *>(this->GetDataImage()->data);
        int numel=this->GetNumberElements();
        float * Data_PTRm=&Data_PTR[m*numel];

        for (int i=0; i<numel; i++,Data_PTRm++) {
            if((*Data_PTRm)<MinResult){ // update of the maximum
                MinResult=(*Data_PTRm);
            }
        }
        //cout<<"the min is "<<MinResult<<endl;
        return MinResult;
    }
}

// Copy the data of modality m according to the mask in a float array and returns the pointer to the beginning
float * TreeEM::GetDataModalMasked(int m){
    if (this->GetDataImage()==NULL) { // no data image attributed in the tree
        cout<< "No data attributed to the tree, tree not valid"<< endl;
        return NULL;
    }
    else{
        // Check if modality wanted is compatible with Data
        if (m>=this->GetNumberModalities()){// modality wanted is out of bounds
            cout<< "Modality wanted is out of bounds"<<endl;
            return NULL;
        }
        int numel=this->GetNumberElements();
        float * Data_PTR=static_cast<float *>(this->GetDataImage()->data);
        float * Data_PTRm=&Data_PTR[m*numel];
        float * DataMaskedModal=new float[numel];//{0};
        for (int i=0; i<numel; i++) {
            DataMaskedModal[i]=0;
        }
        float * DataMaskedModal_PTR=DataMaskedModal;
        int * L2S_PTR=this->GetL2S();
        for (int i=0; i<numel; i++,L2S_PTR++,Data_PTRm++,DataMaskedModal_PTR++) {
            if (*L2S_PTR>=0) {
                *DataMaskedModal_PTR=*Data_PTRm;
            }
        }
        return DataMaskedModal;
    }
}

int TreeEM::GetNumberElements(){
    if (this->GetDataImage()==NULL){
        cout<<"No image to look at"<<endl;
        return 0;
    }
    //cout<< this->GetDataImage()->nx<<endl;
    return this->GetDataImage()->nx*this->GetDataImage()->ny*this->GetDataImage()->nz;
}

int TreeEM::GetNumberModalities(){
    if (this->GetDataImage()==NULL){
        cout<<"No image to look at"<<endl;
        return 0;
    }
    return this->GetDataImage()->nu*this->GetDataImage()->nt;
}

int TreeEM::GetNumberMaskedElements(){
    if(this->GetParent()==NULL){
        if(this->NumberMaskedElements==0){
            this->SetNumberMaskedElements(this->MakeNumberMaskedElements());
        }
        return this->NumberMaskedElements;
    }
    else{
        return this->GetParent()->GetNumberMaskedElements();
    }
}

int TreeEM::MakeNumberMaskedElements(){
    int numel=this->GetDataImage()->nx*this->GetDataImage()->ny*this->GetDataImage()->nz;
    if(this->GetMask()==NULL || this->GetMask()->datatype!=DT_BINARY){
        //cout<<"Mask is NULL or not Binary "<<this->GetNumberElements()<<endl;
        if (this->GetMask()!=NULL && this->GetMask()->datatype!=DT_BINARY) {
            cout<<"Mask is not binary "<< this->GetMask()<<endl;
        }
        return this->GetNumberElements();
    }
    else{

        //cout<<"Mask is correct and binary"<<endl;
        int CountMaskedVoxels=0;
        bool * MaskPointer=static_cast<bool *>(this->GetMask()->data);
        for (int i=0; i<numel; i++) {
            if (MaskPointer[i]>0) {
                CountMaskedVoxels++;
            }
        }
        return CountMaskedVoxels;
    }
}

int * TreeEM::GetL2S(){
    if(this->GetMask()==this->GetMaskDirect()){
        if (this->L2S==NULL) {
            this->MakeL2S();
        }
        return this->L2S;
    }
    else{ // we are at a place where L2S should be NULL and we get it from the parent
        //cout<<"It must be a child"<<endl;
        if (this->L2S!=NULL) {
            delete [] this->L2S;
            this->L2S=NULL;
        }
        return this->GetParent()->GetL2S();
    }
}

void TreeEM::MakeL2S(){
    if(this->GetMask()==this->GetMaskDirect()){
//        cout<<"we are at the basis level"<<endl;
        if (this->L2S!=NULL) {
            delete [] this->L2S;
        }
        int numel=this->GetNumberElements();
        this->L2S= new int [numel];
        if(this->GetMask()!=NULL && this->GetMask()->datatype==DT_BINARY){
            /* The compatibility between the dimensions of DataImage and MaskImage has to be checked when setting the MaskImage pointer beforehand*/
            bool * Maskptr = static_cast<bool *> (this->GetMask()->data);
            bool * Maskptrtmp = Maskptr;
            int * L2S_PTR = this->GetL2S();
            Maskptrtmp = Maskptr;
            int tempindex=0;


            for (int i=0; i<numel; i++,Maskptrtmp++,L2S_PTR++) {
                if ((*Maskptrtmp)>0) {
                    (*L2S_PTR)=tempindex;
                    tempindex++;

                }
                else{
                    (*L2S_PTR)=-1;
                }
            }
        }
        else {
            if(this->GetMask()==NULL){
                cout<< "No mask used : all active"<<endl;
            }
            else{
                cout<< "Mask not binary : all active" << endl;
            }
            int * L2S_PTR = this->L2S;
            int numel=this->GetNumberElements();
            for (int i=0; i<numel; i++,L2S_PTR++) {

                (*L2S_PTR)=i;
            }
        }
    }
    else{
        if(this->L2S!=NULL){
            delete [] this->L2S;
            this->L2S=NULL;
        }
    }
}

int * TreeEM::GetS2L(){
    if (this->GetMask()==this->GetMaskDirect()) { // we are at the level where S2L must be calculated
        // if the pointer to the Mask is NULL, it means that the mask is the image itself and the full image must be considered
        if (this->S2L==NULL) {
            this->MakeS2L();
        }
        return this->S2L;

    }
    else{// S2L must be NULL and we get S2L from the parent
        if (this->S2L!=NULL) {
            delete [] S2L;
            this->S2L=NULL;
        }
        return this->GetParent()->GetS2L();
    }
}

void TreeEM::MakeS2L(){
    if (this->GetMask()==this->GetMaskDirect()) { // we are at the level where S2L must be calculated
        // if the pointer to the Mask is NULL, it means that the mask is the image itself and the full image must be considered
        if(this->GetMask()!=NULL && this->GetMask()->datatype==DT_BINARY){
            bool * Maskptr = static_cast<bool *> (this->GetMask()->data);
            bool * Maskptrtmp = Maskptr;
            int numelmasked=this->GetNumberMaskedElements();

            // reinitialise S2L
            if(this->S2L!=NULL){
                delete [] this->S2L;
            }
            this->S2L= new int[numelmasked];
            int * S2L_PTR = (int *)(this->S2L);
            Maskptrtmp = Maskptr;
            int tempindex=0;
            int numel=this->GetNumberElements();
            for (int i=0; i<numel; i++) {
                if ((*Maskptrtmp)>0) {
                    S2L_PTR[tempindex]=i;
                    tempindex++;
                }
                Maskptrtmp++;
            }
        }
        else {
            printf("the datatype is not binary, all the voxels are considered as active");
            // reinitialise S2L
            if(this->S2L!=NULL){
                delete [] this->S2L;
            }
            int numelmasked=this->GetNumberMaskedElements();
            this->S2L= new int[numelmasked];
            int * S2L_PTR=this->S2L;
            int numel=this->GetNumberElements();
            for(int i=0;i<numel;i++, S2L_PTR++){
                (*S2L_PTR)=i;
            }
        }
    }
    else{
        if(this->S2L!=NULL){
            delete [] this->S2L;
            this->S2L=NULL;
        }
    }
}

int * TreeEM::GetHardSeg(){
    if (this->GetMask()==this->GetMaskDirect()) { // we are at the level where HardSeg must be calculated
        // if the pointer to the Mask is NULL, it means that the mask is the image itself and the full image must be considered
        if (this->HardSeg==NULL) {
            this->MakeHardSeg();
        }
        return this->HardSeg;

    }
    else{// HardSeg must be NULL and we get HardSeg from the parent
        if (this->HardSeg!=NULL) {
            delete [] HardSeg;
            this->HardSeg=NULL;
        }
        return this->GetParent()->GetHardSeg();
    }
}

void TreeEM::MakeHardSeg(){
    if (this->GetMask()==this->GetMaskDirect()) { // we are at the level where HardSeg must be calculated
        // if the pointer to the Mask is NULL, it means that the mask is the image itself and the full image must be considered
            int numelmasked=this->GetNumberMaskedElements();

            // reinitialise HardSeg
            if(this->HardSeg!=NULL){
                delete [] this->HardSeg;
            }
            this->HardSeg= new int[numelmasked];
            int * HardSeg_PTR = (int *)(this->HardSeg);
           vector<float*> NormRespChildrenVector;
           int numbchild=this->GetNumberChildren();
           for(int c=0;c<numbchild;c++){
               NormRespChildrenVector.push_back(this->GetChild(c)->GetNormResp());
           }
           // find for each active voxel index of max NormResp and put it in result of HardSeg;
           for(int i=0;i<numelmasked;i++,HardSeg_PTR++){
               int maxInd=0;
               float maxNormResp=0;
               for(int c=0;c<numbchild;c++){
                   if (NormRespChildrenVector[c][i]>maxNormResp){
                       maxInd=c;
                       maxNormResp=NormRespChildrenVector[c][i];
                   }
               }
               *HardSeg_PTR=maxInd;
           }
    }
    else{
        if(this->HardSeg!=NULL){
            delete [] this->HardSeg;
            this->HardSeg=NULL;
        }
    }
}




void TreeEM::GetLeaves(vector<TreeEM*> & LeavesVector){
    if (this->GetNumberChildren()==0) {
        LeavesVector.push_back(this);
    }
    else{
        int numbchild=this->GetNumberChildren();
        for (int c=0; c<numbchild; c++) {
            this->GetChild(c)->GetLeaves(LeavesVector);
        }
    }
}

vector<TreeEM *> TreeEM::GetAllLeaves(){
    vector<TreeEM *> LeavesVector;
    if(this->GetNumberChildren()==0){
        return LeavesVector;
    }
    else{
        int numbchild=this->GetNumberChildren();
        for (int c=0; c<numbchild; c++) {
            if (this->GetChild(c)->GetNumberChildren()==0) {
                LeavesVector.push_back(this->GetChild(c));
            }
            else {
                vector<TreeEM *> ChildLeaveVector=this->GetChild(c)->GetAllLeaves();
                for (int i=0; i<ChildLeaveVector.size(); i++) {
                    LeavesVector.push_back(ChildLeaveVector[i]);
                }
            }
        }
    }
    return LeavesVector;
}

int TreeEM::GetNumberAllLeaves(){
    int numbAllLeaves=0;
    if (this->GetNumberChildren()==0) {
        numbAllLeaves=0;
        return numbAllLeaves;
    }
    else{
        int numbchild=this->GetNumberChildren();
        for (int c=0; c<numbchild; c++) {
            if (this->GetChild(c)->GetNumberChildren()==0) {
                numbAllLeaves+=1;
            }
            else {
                numbAllLeaves=numbAllLeaves+this->GetChild(c)->GetNumberAllLeaves();
            }
        }
    }
    return numbAllLeaves;
}

// Return the number of free parameters in the model
int TreeEM::GetNumberFreeParameters(){
    int numbFreeParameters=0;
    int numbmodal=this->GetNumberModalities();
    if (this->GetNumberChildren()==0) {
        switch (this->GetDistributionType()) {
        case 0: // Should not be possible that a leaf is considered as a mixture but taken care of anyway
            return numbFreeParameters;
            break;

        default: // Gaussian case is default for leaves
            numbFreeParameters=(int)numbmodal+((numbmodal+1)*numbmodal)/2;
            break;
        }
    }
    else{
        int numbchild=this->GetNumberChildren();
        for (int c=0; c<numbchild; c++) {
            numbFreeParameters=numbFreeParameters+this->GetChild(c)->GetNumberFreeParameters(); // Recursive part
        }
    }
    return numbFreeParameters;
}

//Returns the value of the loglikelihood for the considered model
float TreeEM::GetLogLikelihood(){
    int numelmasked=this->GetNumberMaskedElements();
    PrecisionTYPE LL=0;
    float * Distribution=new float[numelmasked];
    float * Distribution_PTR=Distribution;
    for(int i=0;i<numelmasked;i++){
        Distribution[i]=0;
    }
    this->MakeWeightedDist(Distribution);
    int CountNeedChangeValue=0;


    for (int i=0; i<numelmasked; i++,Distribution_PTR++) {
        if (*Distribution_PTR>=1E-6) {
            LL+=(PrecisionTYPE)log(*Distribution_PTR);
        }
        else {
            CountNeedChangeValue++;
            LL+=(PrecisionTYPE)log(1E-6);
        }
    }
    if(Distribution!=NULL){
        delete [] Distribution;
        Distribution = NULL;
    }
    return (float)LL;
}

float TreeEM::GetLogLikelihoodCEM(int c){
    int numelmasked=this->GetNumberMaskedElements();
    PrecisionTYPE LL=0;
    float * Distribution=new float[numelmasked];
    float * Distribution_PTR=Distribution;
    for(int i=0;i<numelmasked;i++){
        Distribution[i]=0;
    }
    this->MakeWeightedDist(Distribution);
    int CountNeedChangeValue=0;
    int * HardSeg_PTR=this->GetHardSeg();

    for (int i=0; i<numelmasked; i++,Distribution_PTR++,HardSeg_PTR++) {
        if (*Distribution_PTR>=1E-6 && *HardSeg_PTR==c) {
            LL+=(PrecisionTYPE)log(*Distribution_PTR);
        }
        else if(*HardSeg_PTR==c) {
            CountNeedChangeValue++;
            LL+=(PrecisionTYPE)logf(1E-6);
        }
    }
    if(Distribution!=NULL){
        delete [] Distribution;
        Distribution = NULL;
    }
    return (float)LL;
}

vector<TreeEM *> TreeEM::GetAllDirectLeaves(){
    vector<TreeEM *> DirectLeavesVector;
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        if (this->GetChild(c)->GetNumberChildren()==0) {
            DirectLeavesVector.push_back(this->GetChild(c));
        }
    }
    return DirectLeavesVector;
}

int TreeEM::GetNumberDirectLeaves(){
    int numbDirectLeaves=0;
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        if (this->GetChild(c)->GetNumberChildren()==0) {
            numbDirectLeaves++;
        }
    }
    return numbDirectLeaves;
}

void TreeEM::GetNumberLeaves(int & numbleaves){

    if (this->GetNumberChildren()==0) {
        numbleaves++;
    }
    else{
        int numbchild=this->GetNumberChildren();
        for (int c=0; c<numbchild; c++) {
            this->GetChild(c)->GetNumberLeaves(numbleaves);
        }
    }
}

// Returns in a PrecisionTYPE array the shifted values according to the dimension of the shift and its direction
float * TreeEM::GetDataShifted(int dim,int dir){
    // First make sure that the input values are ok or change them accordingly
    dir=dir>0?1:-1;
    dim=dim>2?2:dim;
    dim=dim<0?0:dim;

    // Initialisation of the needed values
    int * Dimensions=new int[3];
    Dimensions[0]=this->GetDataImage()->nx;
    Dimensions[1]=this->GetDataImage()->ny;
    Dimensions[2]=this->GetDataImage()->nz;
    int * Shift=new int[3];
    Shift[0]=1;
    Shift[1]=this->GetDataImage()->nx;
    Shift[2]=this->GetDataImage()->nx*this->GetDataImage()->ny;
    int numbx=this->GetDataImage()->nx;
    int numby=this->GetDataImage()->ny;
    int numbz=this->GetDataImage()->nz;
    int indexShift=Shift[dim]*dir;
    int numbmodal=this->GetNumberModalities();
    int numel=this->GetNumberElements();
    int nvox=this->GetDataImage()->nvox;
    int newIndex=0;
    int oldIndex=0;
    // Allocation of memory for shifted matrix
    float * DataShifted= new float[nvox];//{0};
    for (int i=0; i<nvox; i++) {
        DataShifted[i]=0;
    }
    float * DataPointed= static_cast<float*>(this->GetDataImage()->data);
    for (int m=0; m<numbmodal; m++) {
        for (int i=0; i<numbx; i++) {
            for (int j=0; j<numby; j++) {
                for (int k=0; k<numbz; k++) {
                    newIndex=i+j*Shift[1]+k*Shift[2]+indexShift;
                    oldIndex=i+j*Shift[1]+k*Shift[2];
                    if (newIndex>0 && newIndex<numel) {
                        DataShifted[m*numel+newIndex]=DataPointed[m*numel+oldIndex];
                    }
                }
            }
        }
    }
    delete Dimensions;
    delete Shift;
    return DataShifted;
}

float * TreeEM::GetDataShiftedCorrected(int dim, int dir){
    // First make sure that the input values are ok or change them accordingly
    if(dir!=0){
    dir=dir>0?1:-1;
    }
    dim=dim>2?2:dim;
    dim=dim<0?0:dim;

    // Initialisation of the needed values
    int Dimensions[3];
    Dimensions[0]=this->GetDataImage()->nx;
    Dimensions[1]=this->GetDataImage()->ny;
    Dimensions[2]=this->GetDataImage()->nz;
    int Shift[3];
    Shift[0]=1;
    Shift[1]=this->GetDataImage()->nx;
    Shift[2]=this->GetDataImage()->nx*this->GetDataImage()->ny;
    int numbx=this->GetDataImage()->nx;
    int numby=this->GetDataImage()->ny;
    int numbz=this->GetDataImage()->nz;
    int indexShift=Shift[dim]*dir;
    int numbmodal=this->GetNumberModalities();
    int numel=this->GetNumberElements();
    int nvox=this->GetDataImage()->nvox;
    int numelmasked=this->GetNumberMaskedElements();
    int newIndex=0;
    int oldIndex=0;
    // Allocation of memory for shifted matrix
    float * DataShifted= new float[nvox];//{0};
    int * L2S_PTR=this->GetL2S();
    float * DataPointedNC=this->GetDataBFCorrected();
    float * DataPointed = new float[numel*numbmodal];
    for(int m=0;m<numbmodal;m++){
        L2S_PTR=this->GetL2S();
        float * DataPointedNC_PTR = &DataPointedNC[m*numelmasked];
        for (int i=0; i<numel; i++,L2S_PTR++) {

                DataPointed[i+m*numel]=0;
                if(*L2S_PTR>=0){
                    DataPointed[i+m*numel]=*DataPointedNC_PTR;
                    DataPointedNC_PTR++;
                }
                DataShifted[i+m*numel]=0;
        }
    }

//    nifti_image * Shift0=SavePartialResult(DataPointed, this->GetDataImage(), "/Users/Carole/Documents/PhD/TestShift0");
//    nifti_image_write(Shift0);
//    nifti_image_free(Shift0);
//    Shift0=NULL;


//    float * DataPointed= static_cast<float*>(this->GetDataImage()->data);
    for (int m=0; m<numbmodal; m++) {
        for (int i=0; i<numbx; i++) {
            for (int j=0; j<numby; j++) {
                for (int k=0; k<numbz; k++) {
                    newIndex=i+j*Shift[1]+k*Shift[2]+indexShift;
                    oldIndex=i+j*Shift[1]+k*Shift[2];
                    if (newIndex>0 && newIndex<numel) {
                        DataShifted[m*numel+newIndex]=DataPointed[m*numel+oldIndex];
                    }
                }
            }
        }
    }
//    delete [] Dimensions;
//    Dimensions=NULL;
//    delete [] Shift;
//    Shift = NULL;
    delete DataPointed;
    return DataShifted;
}

float * TreeEM::GetDataCorrelation(int dim){
    int numbmodal=this->GetNumberModalities();
    int numel=this->GetNumberElements();
    float * CorrelationResult=new float[numbmodal];//{0};
    for (int i=0; i<numbmodal; i++) {
        CorrelationResult[i]=0;
    }
    PrecisionTYPE Mean[MaxNumbModal];//{0};
    for (int i=0; i<MaxNumbModal; i++) {
        Mean[i]=0;
    }
    PrecisionTYPE Var[MaxNumbModal];//{0};
    for (int i=0; i<MaxNumbModal; i++) {
        Var[i]=0;
    }
    PrecisionTYPE MeanShifted[MaxNumbModal];//{0};
    for (int i=0; i<MaxNumbModal; i++) {
        MeanShifted[i]=0;
    }
    PrecisionTYPE VarShifted[MaxNumbModal];//{0};
    for (int i=0; i<MaxNumbModal; i++) {
        VarShifted[i]=0;
    }
    dim=dim<0?0:dim;
    dim=dim>2?2:dim;
    int numelmasked=this->GetNumberMaskedElements();

    // Handling the + part same afterwards for the - part
//    float * DataShifted=this->GetDataShifted(dim, 1);
    float * DataShifted=this->GetDataShiftedCorrected(dim,1);
//    nifti_image * ShiftPlus=SavePartialResult(DataShifted, this->GetDataImage(), "/Users/Carole/Documents/PhD/TestShiftPlus");
//    nifti_image_write(ShiftPlus);
//    nifti_image_free(ShiftPlus);
//    ShiftPlus=NULL;
    float * Data=this->GetDataBFCorrected();
//    float * Data=static_cast<float*>(this->GetDataImage()->data);
    float * DataShifted_PTR=DataShifted;
    float * Data_PTR=Data;
    int* L2S_PTR=this->GetL2S();


    // First calculate the mean and the standard deviation for shifted and not shifted
    for (int m=0; m<numbmodal; m++) {
        DataShifted_PTR=&DataShifted[m*numel];
        Data_PTR=&Data[m*numelmasked];
        L2S_PTR=this->GetL2S();
        float tmpMeanShifted=0;
        float tmpMean=0;
        float tmpStd=0;
        float tmpStdShifted=0;
        for (int i=0; i<numel; i++,DataShifted_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
                tmpMeanShifted+=*DataShifted_PTR;
                tmpMean+=*Data_PTR;
                tmpStdShifted+=*DataShifted_PTR*(*DataShifted_PTR);
                tmpStd+=*Data_PTR*(*Data_PTR);
                Data_PTR++;
            }
        }
        MeanShifted[m]=tmpMeanShifted/numelmasked;
        Mean[m]=tmpMean/numelmasked;
//        Var[m]=powf(tmpStd/numelmasked-Mean[m]*Mean[m],0.5);
//        Var[m]=powf(tmpStd/(numelmasked-1)-(float)(numelmasked+1)/(float)(numelmasked-1)*Mean[m]*Mean[m], 0.5);
        Var[m]=tmpStd/numelmasked-Mean[m]*Mean[m];
//        VarShifted[m]=powf(tmpStdShifted/numelmasked-MeanShifted[m]*MeanShifted[m], 0.5);
//        VarShifted[m]=powf(tmpStdShifted/(numelmasked-1)-(float)(numelmasked+1)/(float)(numelmasked-1)*MeanShifted[m]*MeanShifted[m], 0.5);
        VarShifted[m]=tmpStdShifted/numelmasked-MeanShifted[m]*MeanShifted[m];
    }

    // Then calculate the correlation

    for (int m=0; m<numbmodal; m++) {
        L2S_PTR=this->GetL2S();
        DataShifted_PTR=&DataShifted[m*numel];
        Data_PTR=&Data[m*numelmasked];
        float tmpCorrelation=0;
        for (int i=0; i<numel; i++,DataShifted_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
//                tmpCorrelation+=((*Data_PTR)-Mean[m])*((*DataShifted_PTR)-MeanShifted[m]);
                tmpCorrelation+=*Data_PTR*(*DataShifted_PTR);
                Data_PTR++;
            }
        }
//        CorrelationResult[m]=tmpCorrelation/(StdShifted[m]*Std[m]);
        tmpCorrelation=tmpCorrelation-numelmasked*Mean[m]*MeanShifted[m];
        CorrelationResult[m]=tmpCorrelation/sqrt(VarShifted[m]*Var[m]);
    }

    delete [] DataShifted;
    DataShifted=NULL;
    DataShifted_PTR=NULL;

//    float * DataShiftedMoins=this->GetDataShifted(dim, -1);
    float * DataShiftedMoins=this->GetDataShiftedCorrected(dim,-1);
    L2S_PTR=this->GetL2S();
//    nifti_image * ShiftMoins=SavePartialResult(DataShiftedMoins, this->GetDataImage(), "/Users/Carole/Documents/PhD/TestShiftMoins");
//    nifti_image_write(ShiftMoins);
//    nifti_image_free(ShiftMoins);
//    ShiftMoins=NULL;
    float * DataShiftedMoins_PTR=DataShiftedMoins;

    // Then calculate the mean and the standard deviation for shifted in other direction
    for (int m=0; m<numbmodal; m++) {
        DataShiftedMoins_PTR=&DataShiftedMoins[m*numel];
        L2S_PTR=this->GetL2S();
        float tmpMeanShiftedMoins=0;
        float tmpStdShiftedMoins=0;
        int CountActive=0;
        for (int i=0; i<numel; i++,DataShiftedMoins_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
                tmpMeanShiftedMoins+=*DataShiftedMoins_PTR;
                CountActive++;
                tmpStdShiftedMoins+=*DataShiftedMoins_PTR*(*DataShiftedMoins_PTR);
            }
        }
//        cout<< "Numelmasked is "<<CountActive<<endl;
        MeanShifted[m]=tmpMeanShiftedMoins/(float)numelmasked;
//        StdShifted[m]=powf(tmpStdShifted/numelmasked-MeanShifted[m]*MeanShifted[m], 0.5);
//        VarShifted[m]=powf(tmpStdShifted/(numelmasked-1)-(float)(numelmasked+1)/(float)(numelmasked-1)*MeanShifted[m]*MeanShifted[m], 0.5);
        VarShifted[m]=tmpStdShiftedMoins/(float)numelmasked-MeanShifted[m]*MeanShifted[m];
    }
//float test=0;
//float test2=0;
    // Then calculate the correlation
    for (int m=0; m<numbmodal; m++) {
        DataShiftedMoins_PTR=&DataShiftedMoins[m*numel];
        Data_PTR=&Data[m*numelmasked];
        L2S_PTR=this->GetL2S();
        float tmpCorrelationMoins=0;
//        float tmpCorrelationtest=0;
        for (int i=0; i<numel; i++,DataShiftedMoins_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
//                tmpCorrelationMoins+=((*Data_PTR)-Mean[m])*((*DataShiftedMoins_PTR)-MeanShifted[m]);
                tmpCorrelationMoins+=*Data_PTR*(*DataShiftedMoins_PTR);
                Data_PTR++;
            }
        }
//        CorrelationResult[m]=(float)(CorrelationResult[m]+tmpCorrelation/(VarShifted[m]*Var[m]))/(2*(numelmasked-1));
//        CorrelationResult[m]=(float)(CorrelationResult[m]+tmpCorrelation/(VarShifted[m]*Var[m]))/(2*(numelmasked));
//        test=tmpCorrelationMoins/sqrt(Var[m]*VarShifted[m]);
        tmpCorrelationMoins=tmpCorrelationMoins-numelmasked*Mean[m]*MeanShifted[m];
//        test2=test2/sqrt(VarShifted[m]*Var[m]);
        CorrelationResult[m]=(float)(CorrelationResult[m]+tmpCorrelationMoins/sqrt(VarShifted[m]*Var[m]))/(2.0*(float)numelmasked);
    }
    delete [] DataShiftedMoins;
    DataShiftedMoins=NULL;
    DataShiftedMoins_PTR=NULL;
    for(int m=0;m<numbmodal;m++){
        if(fabs(CorrelationResult[m])>1){
            cout<<"Pb with the calculation of the correlation in dim "<<dim<<endl;
        }
    }
    return CorrelationResult;
}

float TreeEM::GetDataCorrelationMod(int dim,int * modalChoice){
    int numbmodal=this->GetNumberModalities();
//    cout<< modalChoice[0]<<modalChoice[1];
    // Check if the choice of modalities are available
    for(int i =0;i<2;i++){
        if(modalChoice[i]>=numbmodal || modalChoice[i]<0){
            cout<<"Pb in the modalities to correlate with "<<modalChoice[0]<<" and "<<modalChoice[1]<<endl;
            return 0;
        }
    }

    int numel=this->GetNumberElements();
    float CorrelationResult=0;
    float Mean=0;
    float Var=0;
    float MeanShifted=0;
    float VarShifted=0;

    dim=dim<0?0:dim;
    dim=dim>2?2:dim;
    int numelmasked=this->GetNumberMaskedElements();

    // Handling the + part same afterwards for the - part
//    float * DataShifted=this->GetDataShifted(dim, 1);
    float * DataShifted;
    if(modalChoice[0]==modalChoice[1]){

     DataShifted=this->GetDataShiftedCorrected(dim,1);
    }
    else{
         DataShifted=this->GetDataShiftedCorrected(dim,0);
    }
    float * Data=this->GetDataBFCorrected();
//    float * Data=static_cast<float*>(this->GetDataImage()->data);
    float * DataShifted_PTR=DataShifted;
    float * Data_PTR=Data;
    int* L2S_PTR=this->GetL2S();


    // First calculate the mean and the standard deviation for shifted and not shifted
        DataShifted_PTR=&DataShifted[modalChoice[1]*numel];
        Data_PTR=&Data[modalChoice[0]*numelmasked];
        L2S_PTR=this->GetL2S();
        float tmpMeanShifted=0;
        float tmpMean=0;
        float tmpStd=0;
        float tmpStdShifted=0;
        for (int i=0; i<numel; i++,DataShifted_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
                tmpMeanShifted+=*DataShifted_PTR;
                tmpMean+=*Data_PTR;
                tmpStdShifted+=*DataShifted_PTR*(*DataShifted_PTR);
                tmpStd+=*Data_PTR*(*Data_PTR);
                Data_PTR++;
            }
        }
        MeanShifted=tmpMeanShifted/numelmasked;
        Mean=tmpMean/numelmasked;
//        Var[m]=powf(tmpStd/numelmasked-Mean[m]*Mean[m],0.5);
//        Var[m]=powf(tmpStd/(numelmasked-1)-(float)(numelmasked+1)/(float)(numelmasked-1)*Mean[m]*Mean[m], 0.5);
        Var=tmpStd/numelmasked-Mean*Mean;
//        VarShifted[m]=powf(tmpStdShifted/numelmasked-MeanShifted[m]*MeanShifted[m], 0.5);
//        VarShifted[m]=powf(tmpStdShifted/(numelmasked-1)-(float)(numelmasked+1)/(float)(numelmasked-1)*MeanShifted[m]*MeanShifted[m], 0.5);
        VarShifted=tmpStdShifted/numelmasked-MeanShifted*MeanShifted;


    // Then calculate the correlation


        L2S_PTR=this->GetL2S();
        DataShifted_PTR=&DataShifted[modalChoice[1]*numel];
        Data_PTR=&Data[modalChoice[0]*numelmasked];
        float tmpCorrelation=0;
        for (int i=0; i<numel; i++,DataShifted_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
//                tmpCorrelation+=((*Data_PTR)-Mean[m])*((*DataShifted_PTR)-MeanShifted[m]);
                tmpCorrelation+=*Data_PTR*(*DataShifted_PTR);
                Data_PTR++;
            }
        }
//        CorrelationResult[m]=tmpCorrelation/(StdShifted[m]*Std[m]);
        tmpCorrelation=tmpCorrelation-numelmasked*Mean*MeanShifted;
        CorrelationResult=tmpCorrelation/sqrt(VarShifted*Var);


    delete [] DataShifted;
    DataShifted=NULL;
    DataShifted_PTR=NULL;

//    float * DataShiftedMoins=this->GetDataShifted(dim, -1);
    float * DataShiftedMoins;
    if(modalChoice[0]==modalChoice[1]){
     DataShiftedMoins=this->GetDataShiftedCorrected(dim,-1);
    }
    else{
         DataShiftedMoins=this->GetDataShiftedCorrected(dim,0);
    }
    L2S_PTR=this->GetL2S();
    float * DataShiftedMoins_PTR=DataShiftedMoins;

    // Then calculate the mean and the standard deviation for shifted in other direction
        DataShiftedMoins_PTR=&DataShiftedMoins[modalChoice[1]*numel];
        L2S_PTR=this->GetL2S();
        float tmpMeanShiftedMoins=0;
        float tmpStdShiftedMoins=0;
        int CountActive=0;
        for (int i=0; i<numel; i++,DataShiftedMoins_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
                tmpMeanShiftedMoins+=*DataShiftedMoins_PTR;
                CountActive++;
                tmpStdShiftedMoins+=*DataShiftedMoins_PTR*(*DataShiftedMoins_PTR);
            }
        }
//        cout<< "Numelmasked is "<<CountActive<<endl;
        MeanShifted=tmpMeanShiftedMoins/(float)numelmasked;
//        StdShifted[m]=powf(tmpStdShifted/numelmasked-MeanShifted[m]*MeanShifted[m], 0.5);
//        VarShifted[m]=powf(tmpStdShifted/(numelmasked-1)-(float)(numelmasked+1)/(float)(numelmasked-1)*MeanShifted[m]*MeanShifted[m], 0.5);
        VarShifted=tmpStdShiftedMoins/(float)numelmasked-MeanShifted*MeanShifted;
//float test=0;
//float test2=0;
    // Then calculate the correlation
        DataShiftedMoins_PTR=&DataShiftedMoins[modalChoice[1]*numel];
        Data_PTR=&Data[modalChoice[0]*numelmasked];
        L2S_PTR=this->GetL2S();
        float tmpCorrelationMoins=0;
//        float tmpCorrelationtest=0;
        for (int i=0; i<numel; i++,DataShiftedMoins_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
//                tmpCorrelationMoins+=((*Data_PTR)-Mean[m])*((*DataShiftedMoins_PTR)-MeanShifted[m]);
                tmpCorrelationMoins+=*Data_PTR*(*DataShiftedMoins_PTR);
                Data_PTR++;
            }
        }

//        CorrelationResult[m]=(float)(CorrelationResult[m]+tmpCorrelation/(VarShifted[m]*Var[m]))/(2*(numelmasked-1));
//        CorrelationResult[m]=(float)(CorrelationResult[m]+tmpCorrelation/(VarShifted[m]*Var[m]))/(2*(numelmasked));
//        test=tmpCorrelationMoins/sqrt(Var[m]*VarShifted[m]);
        tmpCorrelationMoins=tmpCorrelationMoins-numelmasked*Mean*MeanShifted;
//        test2=test2/sqrt(VarShifted[m]*Var[m]);
        CorrelationResult=(float)(CorrelationResult+tmpCorrelationMoins/sqrt(VarShifted*Var))/(2.0*(float)numelmasked);

    delete [] DataShiftedMoins;
    DataShiftedMoins=NULL;
    DataShiftedMoins_PTR=NULL;
        if(fabs(CorrelationResult)>1){
            cout<<"Pb with the calculation of the correlation in dim "<<dim<<" "<<CorrelationResult<<endl;
        }
    return CorrelationResult;
}

// Returns the independence factor for dimension dim;
float TreeEM::MakeIndFactor(int dim){
    dim=dim>2?2:dim;
    dim=dim<0?0:dim;
    int numbmodal=this->GetNumberModalities();
    float * CorrelationDim=this->GetDataCorrelation(dim);
    float CorrMean=0;
    for (int m=0; m<numbmodal; m++) {
        CorrMean+=CorrelationDim[m];
    }
    CorrMean/=numbmodal;
//    cout<< CorrMean <<"corr mean for dim "<<dim<<endl;
    float IndFactorDim=0.9394/powf(-2*logf(2)/logf(CorrMean), 0.5);
    delete [] CorrelationDim;
    CorrelationDim=NULL;
    return IndFactorDim;
}

float TreeEM::MakeIndFactorMod(int dim){
    dim=dim>2?2:dim;
    dim=dim<0?0:dim;
    int numbModal=this->GetNumberModalities();
    // Initialisation of the corr matrix
    float * CorrMatrix=new float[numbModal*numbModal];
    for(int i=0;i<numbModal;i++){
        for(int j=0;j<numbModal;j++){
            if(i==j){
                CorrMatrix[i+j*numbModal]=1;
            }
            else{
                CorrMatrix[i+j*numbModal]=0;
            }
        }
    }
    int * modalChoice=new int[2];
    for(int m1=0;m1<numbModal;m1++){
        for(int m2=m1;m2<numbModal;m2++){
            modalChoice[0]=m1;
            modalChoice[1]=m2;
//            cout<< modalChoice[0]<<modalChoice[1];
            CorrMatrix[m1+m2*numbModal]=this->GetDataCorrelationMod(dim,modalChoice);
            CorrMatrix[m2+m1*numbModal]=CorrMatrix[m1+m2*numbModal];
        }
    }
    delete [] modalChoice;
    modalChoice=NULL;
    float DetCorr= determinant(CorrMatrix,numbModal);
    DetCorr=1;
    for(int m=0;m<numbModal;m++){
        DetCorr*=CorrMatrix[m+m*numbModal];
    }
    cout<< DetCorr << "det for dimension "<<dim<<endl;
    float IndFactorDim=0.9394/powf(-2*logf(2)/logf(fabs(DetCorr)), 0.5);
    return IndFactorDim;

}

float TreeEM::MakeIndFactorTotMod(){
    float IF0=this->MakeIndFactorMod(0);
    float IF1=this->MakeIndFactorMod(1);
    float IF2=this->MakeIndFactorMod(2);
    cout<<"the partial IF are "<< IF0<<" "<<IF1<<" "<<IF2<<endl;
    return IF0*IF1*IF2;
//    return this->MakeIndFactorMod(0)*this->MakeIndFactorMod(1)*MakeIndFactorMod(2);
}

float TreeEM::MakeIndFactorTot(){
    float IF0=this->MakeIndFactor(0);
    float IF1=this->MakeIndFactor(1);
    float IF2=this->MakeIndFactor(2);
    cout<<"the partial IF are "<< IF0<<" "<<IF1<<" "<<IF2<<endl;
    return IF0*IF1*IF2;
//    return this->MakeIndFactor(0)*this->MakeIndFactor(1)*MakeIndFactor(2);
}

float TreeEM::GetIndFactor(){
    if(this->GetParent()==NULL){
        if(this->IndFactor<=0){
            this->SetIndFactor(this->MakeIndFactorTot());
            return this->IndFactor;
        }
        else return this->IndFactor;
    }
    else{
        return this->GetParent()->GetIndFactor();
    }
}

float * TreeEM::GetDPChildren(){
    if(!this->IsRoot()){
        return this->GetParent()->GetDPChildren();
    }
    else{
        return this->DPChildren;
    }
}

float * TreeEM::GetDPChildrenDirect(){
    return this->DPChildren;
}

bool TreeEM::GetFlagDistClassInd(){
    if(!this->IsRoot()){
       return this->GetParent()->GetFlagDistClassInd();
    }
    else{
        return this->FlagDistClassInd;
    }
}

int TreeEM::GetFlagOutliers(){
    if(!this->IsRoot()){
        return this->GetParent()->GetFlagOutliers();
    }
    else{
        return this->FlagOutliers;
    }
}

// Returns a float array of size the number of elements in the image (in order to save it as a nifti image afterwards) with the needed values
float * TreeEM::GetPartialResult(PartialResultType ResultType, SEG_PARAMETERS * segment_param){
    if (!this->IsTreeValid()) {
        cout<<"Tree not valid : will not get any partial result"<<endl;
        return NULL;
    }
    if (ResultType>5 || ResultType<0) {
        cout<<"Not proper result type wanted"<<endl;
        return NULL;
    }
    int numel=this->GetNumberElements();
    int numelmasked=this->GetNumberMaskedElements();
    float * PartialResult=new float[numel];//{0};
    for (int i=0; i<numel; i++) {
        PartialResult[i]=0;
    }
    float * PartialResult_PTR=PartialResult;
    float * tmp=new float[numelmasked];
    float * tmp_PTR=tmp;
    for(int i=0;i<numelmasked;i++){
        tmp[i]=0;
    }
    int * L2S_PTR=this->GetL2S();
    float * tmp2;

    // Get the pointer to the data to put in the PartialResult array;
    switch (ResultType) {
    case DISTRIBUTION:
        this->MakeWeightedDist(tmp);
        break;
    case NONNORMRESP:
        this->MakeNonNormWeightedSum(tmp,segment_param);
        break;
    case NORMRESP:
        tmp2=this->GetNormResp();
        for(int i=0;i<numelmasked;i++,tmp2++){
            tmp[i]=*tmp2;
        }
        break;
    default:
        tmp2=static_cast<float *>(this->GetDataImage()->data);
        for(int i=0;i<numel;i++,tmp2++,L2S_PTR++){
            if(*L2S_PTR>=0){
                *tmp_PTR=*tmp2;
                tmp_PTR++;
            }
        }
        tmp_PTR=tmp;
        break;
    }
    L2S_PTR=this->GetL2S();
    // Copying the corresponding data at the correct position in the PartialResult array
    for (int i=0; i<numel; i++,PartialResult_PTR++,L2S_PTR++) {
        if (*L2S_PTR>=0) {
            *PartialResult_PTR=*tmp_PTR;
            tmp_PTR++;
        }
    }
    if(tmp!=NULL){
        delete [] tmp;
        tmp=NULL;
    }
    return PartialResult;
}

//float * TreeEM::GetBFCorrectionDirect(){
//    return this->BFCorrection;
//}

float * TreeEM::GetBFCoeffsDirect(){
    return this->BFCoeffs;
}

float * TreeEM::GetDataBFCorrectedDirect(){
    return this->DataBFCorrected;
}

//float * TreeEM::GetBFCorrection(){
//    if (this->WhatTypeOfTreeComponent()!=ROOT && this->GetBFCorrectionDirect()==NULL) {
//        return this->GetParent()->GetBFCorrection();
//    }
//    else{
//        // if not yet calculated when needed, then create the right Basis Functions
//        if (BFFlag && this->BFCorrection==NULL) {
//            this->BFCorrection=this->MakeBFCorrection();
//        }
//        return this->BFCorrection;
//    }
//}

float * TreeEM::GetDataBFCorrected(){
    if (this->WhatTypeOfTreeComponent()!=ROOT && this->GetDataBFCorrectedDirect()==NULL) {
        return this->GetParent()->GetDataBFCorrected();
    }
    else{
        // if not yet calculated when needed, then create the right Basis Functions
        if (this->DataBFCorrected==NULL) {
            float * DataCorrectedToSet=this->MakeDataBFCorrected();
            this->SetDataBFCorrected(DataCorrectedToSet);
            if(DataCorrectedToSet!=NULL){
                delete [] DataCorrectedToSet;
                DataCorrectedToSet=NULL;
            }
        }
        return this->DataBFCorrected;
    }
}

float * TreeEM::GetBFCoeffs(){
    if ((this->WhatTypeOfTreeComponent()!=ROOT && this->WhatTypeOfTreeComponent()!=INITIALNODE) && !this->AreBFCoeffsDirectUseable()) {
        return this->GetParent()->GetBFCoeffs();
    }
    else {
        //        if (BFFlag && !this->AreBFCoeffsDirectUseable()) {
        ////            if (this->GetBFCorrection()==NULL) {
        ////                this->BFCorrection=this->MakeBFCorrection();
        ////            }
        //            this->BFCoeffs=this->MakeFinalBFCoeffsChildren();
        //        }
        return this->BFCoeffs;
    }
}

/************** SET FUNCTIONS WITH NEEDED CHECKS *********************/

void TreeEM::SetParent(TreeEM *ParentInput){
    // First check if the ParentInput is a valid tree (Method to write, not done so far)
    // then change the pointer of parent
    this->Parent=ParentInput;
    // Reinitialise memory allocation and set everything back to zero
    this->ReinitialiseTree();
}

void TreeEM::ReinitialiseTree(){
    this->SetMask(this->GetMask());
    for (int c=0; c<this->GetNumberChildren(); c++) {
        this->GetChild(c)->ReinitialiseTree();
    }
}

void TreeEM::SetChildren(vector<TreeEM *> ChildrenInput){
    // First check if everything is right with each of the offered tree in the vector
    for (int c=0; c<ChildrenInput.size(); c++) { /*if valid, each in ChildrenInput is added to the Children vector of this
                                                  WARNING : here the previous possible existing children are not suppressed see */
        if (ChildrenInput[c]->IsTreeValid() && ChildrenInput[c]->GetParent()==this) {
            this->Children.push_back(ChildrenInput[c]);
        }
    }
    if (ChildrenInput.size()==0) { // if there is no children as input, be sure that this is still in mode simple distribution
        if (this->GetNumberChildren()==0) {
            this->CreateAllocateAndInitializeParameters(1);
        }
    }
    else{
        this->MakeParametersMixture();
    }
}

// Add Child if valid to current vector of Children
void TreeEM::AddChild(TreeEM * ChildAdded){
    if (ChildAdded->GetParent()==this) { // check that the Child added is set for the right parent
        this->Children.push_back(ChildAdded);
        if (!this->IsNumbChildOKWithDistType()) {
            this->MakeParametersMixture();
        }
    }
    else{
        cout<<"The added child does not correspond to the right parent"<<endl;
    }
}

void TreeEM::SetData(nifti_image * DataInput){
    // Check if the considered Tree on which applied is of type Root otherwise data cannot be changed nor modified
    if (this->WhatTypeOfTreeComponent()!=ROOT&&this->WhatTypeOfTreeComponent()!=INITIALNODE) {
        cout<< "Data image can only be modified at root level and pointer is only stored there"<<endl;
        return;
    }
    else{ // the data can be modified since we are at Root level;
        this->DataImage=DataInput;
        if(!IsDataFloat()){
            this->ConvertDataImageToFloat();
        }
//        cout<< "Re give data input"<<endl;
//        if (!this->IsDataImageNormalised()) {
//            cout<< "We have to normalise Image"<<endl;
//            this->NormaliseDataImage();
//        }
//        cout<< "Image is already normalised" << endl;
        // if there was previously a mask, check if we can keep it
        if(this->GetMask()!=NULL){ // check on the dimensions
            if (this->GetMask()->nx!=DataInput->nx || this->GetMask()->ny!=DataInput->ny || this->GetMask()->nz!=DataInput->nz){ // Mask reput to NULL as default and so everything is reallocated as needed
                this->SetMask(NULL);
                this->ReinitialiseTree();
            }
            else { // in case of existing Mask, need to reinitialise the tree : needed mostly for the parameters structure in case the number of modalities has changed. Allows for keeping a similar structure but reinitialising and reallocating memory
                this->ReinitialiseTree();
            }
        }
    }
}

void TreeEM::NormaliseDataImage(){
    if (this->GetDataImage()==NULL){
        cout<<"WARNING : there is no data to segment !"<<endl;
        return;
    }
    if (!this->IsDataFloat()) {
        this->ConvertDataImageToFloat();
    }
    int numbmodal=this->GetNumberModalities();
    int numel=this->GetNumberElements();
//    cout<<"Data type image"<<this->GetDataImage()->datatype<< endl;
    float * Data_PTR=static_cast<float *>(this->GetDataImage()->data);
    float * Data_PTRtmp=Data_PTR;
    for (int m=0; m<numbmodal; m++) {
        Data_PTRtmp=&Data_PTR[m*numel];
        if (this->GetMaxDataModal(m)==this->GetMinDataModal(m)) { // Case where all the values are equal in the image
            for (int i=0; i<numel; i++,Data_PTRtmp++) {
                (*Data_PTRtmp)=0;
            }
        }
        else{
            float MultiplicativeFactor=(expf(1.0f)-1.0f)/(this->GetMaxDataModal(m)-this->GetMinDataModal(m));
            float AdditiveFactor=1.0f-MultiplicativeFactor*this->GetMinDataModal(m);
            //cout<< MultiplicativeFactor << "and " <<AdditiveFactor << endl;
            int numel=this->GetNumberElements();
            for (int i=0; i<numel; i++,Data_PTRtmp++) {
                (*Data_PTRtmp)=MultiplicativeFactor*(*Data_PTRtmp)+AdditiveFactor;
                (*Data_PTRtmp)=logf((*Data_PTRtmp)); // Normally then between 0 and 1;
            }
        }
    }
    //cout<<this->GetMinDataModal(0)<<"and "<<this->GetMaxDataModal(0);
    return;
}

void TreeEM::QuantilizeDataImage(SEG_PARAMETERS * segment_param){
    // Check if there is Data to modify
    if(this->GetDataImage()==NULL){
        cout<<"No data to take care of"<<endl;
    }
    // Make image properly float
    if(!this->IsDataFloat()){
        this->ConvertDataImageToFloat();
    }
    // Look at values of the quantiles to choose as outliers limits
    float quantMin=segment_param->quantMin;
    float quantMax=segment_param->quantMax;
    quantMin=quantMin<=0?0:quantMin;
    quantMax=quantMax>=1?1:quantMax;
    if(quantMax<=quantMin){
        cout<<"quantiles not in right order therefore nothing done"<<endl;
        return;
    }
    // Copy Data image values to define quantiles values for each modality;
    int numbmodal=this->GetNumberModalities();
    int numel = this->GetNumberElements();
    int numelmasked=this->GetNumberMaskedElements();
    int IndQuantMin=(int)round(quantMin*numelmasked);
    IndQuantMin=IndQuantMin<=0?0:IndQuantMin;
    int IndQuantMax=(int)round(quantMax*numelmasked);
    IndQuantMax=IndQuantMax>=numelmasked-1?(numelmasked-1):IndQuantMax;

    for(int m=0;m<numbmodal;m++){
        float * DataToCopy=&static_cast<float*>(this->GetDataImage()->data)[m*numel];
        float * DataToCopy_PTR=DataToCopy;
        float * tmpCopyData=new float[numelmasked];
        float * tmpCopyData_PTR=tmpCopyData;
        int * L2S_PTR=this->GetL2S();
        for(int i=0;i<numel;i++,L2S_PTR++){
            if(*L2S_PTR>=0){
                *tmpCopyData_PTR=*DataToCopy_PTR;
                tmpCopyData_PTR++;
                DataToCopy_PTR++;
            }
        }
        // Sorting of the values in the image to consider
        if(IndQuantMin >0 || IndQuantMax < numelmasked -1){
        quickSort(tmpCopyData,numelmasked);
        float ValueQuantMin=tmpCopyData[IndQuantMin];
        float ValueQuantMax=tmpCopyData[IndQuantMax];
        // Putting all values according to the quantiles chosen.
        L2S_PTR=this->GetL2S();
        for(int i=0;i<numel;i++,L2S_PTR++,DataToCopy++){
            if(*L2S_PTR>=0){
                *DataToCopy=*DataToCopy<=ValueQuantMin?ValueQuantMin:*DataToCopy;
                *DataToCopy=*DataToCopy>=ValueQuantMax?ValueQuantMax:*DataToCopy;
            }
        }
        }
        delete [] tmpCopyData;
        tmpCopyData=NULL;

    }


}

void TreeEM::ConvertDataImageToFloat(){
    if (this->GetDataImage()==NULL) {
        cout<<"WARNING : there is no data to convert !"<<endl;
        return;
    }
    if (this->IsDataFloat()) {
        cout<<"Data already of type float no need to change"<<endl;
        return;
    }


    cout<<"Image has to be converted"<<endl;
//    // First changing the values to have only 0 and 1
//    float * Data_PTR=static_cast<float *>(this->GetDataImage()->data); // consider that then data is float
//    int numbvox=this->GetDataImage()->nvox;

//    // Then changing datatype
//    this->GetDataImage()->datatype=DT_FLOAT;

//    // the initial array is saved and freeed
//    float *initialValue = new float[this->GetDataImage()->nvox];
//    float * initialValue_PTR=initialValue;
//    for (int i=0; i<numbvox; i++,initialValue_PTR++,Data_PTR++) {
//        *initialValue_PTR=(float)(*Data_PTR);
//    }

//    float *initialValue_tmp=initialValue;
//    //cout<<initialValue<<endl;
//    // the new array is allocated and then filled

//    free(this->GetDataImage()->data);
//    this->GetDataImage()->nbyper = sizeof(float);
//    this->GetDataImage()->data = (void *)calloc(this->GetDataImage()->nvox,sizeof(float));
//    float *dataPtr = static_cast<float *>(this->GetDataImage()->data);
//    for (int i=0; i<numbvox; i++, dataPtr++,initialValue_tmp++) {
//        (*dataPtr)=(float)(*initialValue_tmp);
//    }
//    delete [] initialValue;
    seg_changeDatatype<float>(this->GetDataImage());
//    float * TestImage=static_cast<float*>(this->GetDataImage()->data);
//    SaveTmpResult(TestImage,"/Users/Carole/Documents/PhD/TestImage3.nii.gz");
//    nifti_set_filenames(this->GetDataImage(),"/Users/Carole/Documents/PhD/TestImage2.nii.gz",0,0);
//    nifti_image_write(this->GetDataImage());
    //cout<<this->GetDataImage()->datatype<<endl;
}



// After checking set the parameters to ParametersSet
void TreeEM::SetParameters(Parameters * ParametersToSet){

    // Check for compatibility between the Distribution Type and the Parameters and between the image and the Parameters
    if(!CheckForValidityOfParametersStructure(ParametersToSet)||!CheckForSizeParametersValidity( ParametersToSet)){
        cout<< "The structure of the parameters to set is not valid or the sizes are not compatible therefore creation and allocation of valid ones"<<endl;
        if(ParametersToSet !=NULL){
            delete ParametersToSet;
            ParametersToSet=NULL;
        }
        int TypeTreeComponent=this->WhatTypeOfTreeComponent();
        ParametersToSet=new Parameters();
        switch (TypeTreeComponent) {
        case 0: // Case of a root
            this->CreateAllocateAndInitializeParameters(0);
            break;
        case 1:// Case of a branch
            this->CreateAllocateAndInitializeParameters(0);
            break;
        default: // Case of a leaf;
            this->CreateAllocateAndInitializeParameters( 1);
            break;
        }
    }
    else{
        // All checks are ok : Delete previous one and set new one instead;
        if(this->ParametersDistribution!=NULL){
            delete this->ParametersDistribution;
        }
        this->ParametersDistribution=new Parameters;
        this->ParametersDistribution->DistributionType=ParametersToSet->DistributionType;
        this->ParametersDistribution->SizeParameters=ParametersToSet->SizeParameters;
        this->ParametersDistribution->ValueParameters=new float [this->ParametersDistribution->SizeParameters];//{0};
        for (int i=0; i<this->ParametersDistribution->SizeParameters; i++) {
            this->ParametersDistribution->ValueParameters[i]=0;
        }
        int SP=this->ParametersDistribution->SizeParameters;
        for (int sp=0; sp<SP; sp++) {
            this->ParametersDistribution->ValueParameters[sp]=ParametersToSet->ValueParameters[sp];
        }
    }
    return;
}

// Transforms the Parameters into a mixture (needed when creating child for a leaf becoming a branch)
void TreeEM::MakeParametersMixture(){
    if (CheckForValidityOfParametersStructure()&& this->GetDistributionType()==0) {
        return;
    }
    if(this->ParametersDistribution!=NULL){
        if (this->GetParametersValue()!=NULL) {
            delete [] this->GetParametersValue();
            this->ParametersDistribution->ValueParameters=NULL;
        }
        delete this->ParametersDistribution;
        this->ParametersDistribution=NULL;
    }
    this->CreateAllocateAndInitializeParameters(0);
}

void TreeEM::SetIndFactor(float IndFactorInput){
    if(this->GetParent()!=NULL){
        cout<<"Cannot change ind factor if not at root"<<endl;
        return;
    }
    if(IndFactorInput<=0){
        cout<<"IndFactor cannot be negative"<<endl;
        return;
    }
    else{
        this->IndFactor=IndFactorInput;
    }
}

// Set Mask to MaskImage after checking. For the moment : Mask can only be introduced at the root level and nowhere else
void TreeEM::SetMask(nifti_image *MaskImage){
    // First check if trying to set a mask at root or somewhere else
    if (this->GetParent()==NULL){
        cout<<"It is a root or an initial node"<<endl;
        // Then check for the dimension validity
        if (MaskImage!=NULL && (MaskImage->nx!=this->GetDataImage()->nx || MaskImage->ny!=this->GetDataImage()->ny ||MaskImage->nz!=this->GetDataImage()->nz) ) {
            cout<< "The dimensions between Mask and Image are not compatible"<<endl;
            this->SetMask(NULL);
            return;
        }

        // Set the pointer to Mask to MaskImage
        this->Mask=MaskImage;
        if (MaskImage!=NULL && this->GetMask()->datatype!=DT_BINARY) {
//            nifti_set_filenames(this->GetMask(),"/Users/Carole/Documents/PhD/Masktest.nii.gz",0,0);
//            nifti_image_write(this->GetMask());
            this->MakeMaskBinary();
//            nifti_set_filenames(this->GetMask(),"/Users/Carole/Documents/PhD/Masktest.nii.gz",0,0);
//            nifti_image_write(this->GetMask());
//            bool * MaskTest=static_cast<bool*>(this->GetMask()->data);
//            float * MaskToSave=new float[this->GetMask()->nvox];
//            int numel=this->GetMask()->nvox;
//            for(int i=0;i<numel;i++){
//                MaskToSave[i]=(float)MaskTest[i];
//            }
//            SaveTmpResult(MaskToSave,"/Users/Carole/Documents/PhD/Masktest2.nii.gz");

            cout<<"Binarisation performed "<<endl;
        }
        // Calculate S2L and L2S
        if (this->L2S!=NULL) {
            delete [] this->L2S;
            this->L2S=NULL;
        }
        this->MakeL2S();
        if(this->S2L!=NULL){
            delete [] this->S2L;
            this->S2L=NULL;
        }
        this->MakeS2L();
        this->SetNumberMaskedElements(this->MakeNumberMaskedElements());
        int numelmasked=this->GetNumberMaskedElements();
        //cout<< "the number of masked elements is "<<this->GetNumberMaskedElements()<<endl;
        // Allocate right amount of memory to Distribution, NonNormResp, NormResp
        //        if(this->GetDistribution()!=NULL){
        //            delete this->Distribution;
        //        }
        //        this->Distribution=new float[numelmasked];//{0};
        //        for (int i=0; i<numelmasked; i++) {
        //                    this->Distribution[i]=0;
        //                }
        //        if(this->NonNormResp!=NULL){
        //            delete this->NonNormResp;
        //        }
        //        this->NonNormResp=new float[numelmasked];//{0};
        //        for (int i=0; i<numelmasked; i++) {
        //                    this->NonNormResp[i]=0;
        //                }
        if(this->NormResp!=NULL){
            delete this->NormResp;
            this->NormResp=NULL;
        }
        this->NormResp=new float[numelmasked];//{0};
        for (int i=0; i<numelmasked; i++) {
            this->NormResp[i]=0;
        }
//        if(this->HardSeg!=NULL){
//            delete [] this->HardSeg;
//            this->HardSeg=NULL;
//        }
//        this->HardSeg=new int[numelmasked];
//        for(int i=0;i<numelmasked;i++){
//            this->HardSeg[i]=0;
//        }
        // Propagate Mask to children
        int numbchild=this->GetNumberChildren();
        for(int c=0;c<numbchild;c++){
            this->GetChild(c)->SetMask(MaskImage);
        }
        return;
    }
    else {
        //cout<<"Parent Not NULL and address Mask is"<< this->GetMask()<<endl;
        //        this->SetNumberMaskedElements(this->MakeNumberMaskedElements());
        int numelmasked=this->GetNumberMaskedElements();
        //cout<< "the number of masked elements is "<<this->GetNumberMaskedElements()<<endl;
        // Allocate right amount of memory to Distribution, NonNormResp, NormResp
        //        if(this->GetDistribution()!=NULL){
        //            delete this->Distribution;
        //        }
        //        this->Distribution=new float[numelmasked];//{0};
        //        for (int i=0; i<numelmasked; i++) {
        //                    this->Distribution[i]=0;
        //                }
        //        if(this->NonNormResp!=NULL){
        //            delete this->NonNormResp;
        //        }
        //        this->NonNormResp=new float[numelmasked];//{0};
        //        for (int i=0; i<numelmasked; i++) {
        //                    this->NonNormResp[i]=0;
        //                }
        if(this->NormResp!=NULL){
            delete [] this->NormResp;
        }
        this->NormResp=new float[numelmasked];//{0};
        for (int i=0; i<numelmasked; i++) {
            this->NormResp[i]=0;
        }

        // Ensuring that pointer to L2S and S2L and HardSeg are NULL if we are not at the node where the pointer to the Mask is stored
        if (this->L2S!=NULL) {
            delete [] this->L2S;
            this->L2S=NULL;
        }

        if(this->S2L!=NULL){
            delete [] this->S2L;
            this->S2L=NULL;
        }
//        if(this->HardSeg!=NULL){
//            delete [] this->HardSeg;
//            this->HardSeg=NULL;
//        }

        // Propagate memory allocation to children
        int numbchild=this->GetNumberChildren();
        for(int c=0;c<numbchild;c++){
            this->GetChild(c)->SetMask(MaskImage);
        }
    }
    return;
}

void TreeEM::SetNumberMaskedElements(int NumberMaskedElementsInput){
    if(this->GetParent()!=NULL){
        //        cout<<"Cannot set the number of masked elements at any other place than root"<<endl;
        return;
    }
    if(NumberMaskedElementsInput!=this->MakeNumberMaskedElements()){
        cout<<"Not appropriate Input for number masked elements"<<endl;
        return;
    }
    this->NumberMaskedElements=NumberMaskedElementsInput;
    return;
}

void TreeEM::MakeMaskBinary(){
    if (this->GetMask()==NULL){
        cout<< "No mask so nothing to binarise"<<endl;
        return;
    }
    if(this->GetMask()->datatype==DT_BINARY){
        cout<< "Already binarised"<<endl;
        return;
    }
    else{
        cout<<"Mask has to be binarised"<<endl;
//        // First changing the values to have only 0 and 1
//        float * Mask_PTR=static_cast<float *>(this->GetMask()->data); // consider that then data is float
//        float * Mask_PTRtmp=Mask_PTR;
//        int numel=this->GetNumberElements();
//        // if set as mask then it has been checked before that dimensions with dataimage are compatible
//        bool *initialValue = new bool[this->GetMask()->nvox];
//        bool *initialValue_PTR=initialValue;
//        for (int i=0; i<numel; i++,Mask_PTRtmp++,initialValue_PTR++) {
//            (*Mask_PTRtmp)=(*Mask_PTRtmp)>0?1:0;// Set to 1 if strictly positive, 0 otherwise
//            (*initialValue_PTR)=(bool)(*Mask_PTRtmp);

//        }
//        // Then changing datatype
//        this->GetMask()->datatype=DT_BINARY;

//        // the initial array is saved and freeed

//        initialValue_PTR=initialValue;
//        free(this->GetMask()->data);
//        this->GetMask()->nbyper = sizeof(bool);
//        this->GetMask()->data = (void *)calloc(this->GetMask()->nvox,sizeof(bool));
//        bool *dataPtr = static_cast<bool *>(this->GetMask()->data);
//        for (int i=0; i<numel; i++, dataPtr++,initialValue_PTR++) {
//            (*dataPtr)=(bool)(*initialValue_PTR);
//        }
//        delete [] initialValue;
//        cout<<this->GetMask()->datatype<<endl;
        seg_convert2binary(this->GetMask(),0.1);
    }

}

// Setting the Priors to PriorsInput first checking the validity of the dimension and datatype, deleting any previous set of the priors and making sure then they are of probability type
void TreeEM::SetPriors(nifti_image *PriorsInput){
    if (PriorsInput == NULL) {
        this->Priors=NULL;
        return;
    }
    else{
        bool ValidityPriors=(PriorsInput->nx==this->GetDataImage()->nx)&&(PriorsInput->ny==this->GetDataImage()->ny)&&(PriorsInput->nz==this->GetDataImage()->nz);
        if(ValidityPriors){
            if(this->Priors!=NULL){
//                delete [] this->Priors;
                nifti_image_free(this->Priors);
                this->Priors=NULL;
            }
            this->Priors=PriorsInput;
            if((PriorsInput->datatype!=DT_FLOAT32)){
                this->MakePriorsFloat();
            }
            this->MakePriorsProbabilityType();
        }
        return;
    }
}

void TreeEM::SetPriorsAdapted(float * PriorsAdaptedInput){
    if (PriorsAdaptedInput == NULL) {
        this->PriorsAdapted=NULL;
        return;
    }
    else{
        int numel =this->GetNumberElements();
        if(this->GetPriorsAdapted()!=NULL){
            delete [] this->PriorsAdapted;
            this->PriorsAdapted=NULL;
        }
        this->PriorsAdapted=new float[numel];
        for(int i=0;i<numel;i++){
            this->PriorsAdapted[i]=PriorsAdaptedInput[i];
        }
        return;
    }
}

// Change the priors referred so that we only handle float datatype
void TreeEM::MakePriorsFloat(){
    if (this->GetPriors()==NULL) {
        cout<<"WARNING : there is no data to convert !"<<endl;
        return;
    }
    if (this->IsPriorsFloat()) {
        cout<<"Priors already of type float no need to change"<<endl;
        return;
    }


    cout<<"Priors has to be converted"<<endl;
//    // First changing the values to have only 0 and 1
//    float * Priors_PTR=static_cast<float *>(this->GetPriors()->data); // consider that then data is float
//    int numbvox=this->GetPriors()->nvox;

//    // Then changing datatype
//    this->GetPriors()->datatype=DT_FLOAT;

//    // the initial array is saved and freeed
//    float *initialValue = new float[numbvox];
//    float * initialValue_PTR=initialValue;
//    float maxValue=0;
//    for (int i=0; i<numbvox; i++,initialValue_PTR++,Priors_PTR++) {

//        *initialValue_PTR=(float)(*Priors_PTR);
////        cout<<*initialValue_PTR<<"  ";
//        if (*initialValue_PTR>maxValue){
//            maxValue=*initialValue_PTR;
//        }
//    }
//    float *initialValue_tmp=initialValue;
//    if(maxValue>0){
//    for(int i=0;i<numbvox;i++,initialValue_tmp++)
//        *initialValue_tmp/=maxValue;
//    }

//    initialValue_tmp=initialValue;
//    //cout<<initialValue<<endl;
//    // the new array is allocated and then filled

//    free(this->GetPriors()->data);
//    this->GetPriors()->nbyper = sizeof(float);
//    this->GetPriors()->data = (void *)calloc(this->GetPriors()->nvox,sizeof(float));
//    float *priorsPtr = static_cast<float *>(this->GetPriors()->data);
//    for (int i=0; i<numbvox; i++, priorsPtr++,initialValue_tmp++) {
//        (*priorsPtr)=(float)(*initialValue_tmp);
//    }
//    delete [] initialValue;
    seg_changeDatatype<float>(this->GetPriors());
    return;
//    string filenameOut="/Users/Carole/Documents/PhD/Phantom/PriorsModified";
//    nifti_set_filenames(this->GetPriors(), filenameOut.c_str(), 0, 0);
//    nifti_image_write(this->GetPriors());
    //cout<<this->GetPriors()->datatype<<endl;
}

// Replaces the values in Priors so that they are between 0 and 1 and can be used as probabilities by rescaling the dispersion of the values
void TreeEM::MakePriorsAdaptedProbabilityType(){
    if (this->GetPriorsAdapted()==NULL) {
        return;
    }
    float * Priors_PTR=this->GetPriorsAdapted();
    float * Priors_PTRtmp=Priors_PTR;
    // Initialisation of the maximal and minimal values
    float maxPriors=-1E32;
    float minPriors=1E32;
    int numel=this->GetNumberElements();
    // Determination of the value of the maximum and the minimum
    for (int i=0; i<numel; i++,Priors_PTRtmp++) {
        if ((*Priors_PTRtmp)>maxPriors) { // update of the maximum value
            maxPriors=(*Priors_PTRtmp);
        }
        if((*Priors_PTRtmp)<minPriors){ // update of the minimum value
            minPriors=(*Priors_PTRtmp);
        }
    }

    // reinitialisation of the pointer to the beginning of the Priors
    Priors_PTRtmp=Priors_PTR;
    // modification of the values
    if (minPriors==maxPriors) { // case (strange) where all the values are equal, then we put them all to 1
        for (int i=0; i<numel; i++,Priors_PTRtmp++) {
            (*Priors_PTRtmp)=1;
        }
    }
    else {// rescaling of the values
        for (int i=0;i<numel;i++,Priors_PTRtmp++){ // Putting all the values between 0 and 1
            (*Priors_PTRtmp)=(*Priors_PTRtmp-minPriors)/(maxPriors-minPriors);
        }
    }
//    string filenameOut="/Users/Carole/Documents/PhD/Phantom/PriorsModified";
//    nifti_set_filenames(this->GetPriors(), filenameOut.c_str(), 0, 0);
//    nifti_image_write(this->GetPriors());
}


// Replaces the values in Priors so that they are between 0 and 1 and can be used as probabilities by rescaling the dispersion of the values
void TreeEM::MakePriorsProbabilityType(){
    if (this->GetPriors()==NULL) {
        return;
    }
    float * Priors_PTR=static_cast<float *>(this->GetPriors()->data);
    float * Priors_PTRtmp=Priors_PTR;
    // Initialisation of the maximal and minimal values
    float maxPriors=-1E32;
    float minPriors=1E32;
    int numel=this->GetNumberElements();
    // Determination of the value of the maximum and the minimum
    for (int i=0; i<numel; i++,Priors_PTRtmp++) {
        if ((*Priors_PTRtmp)>maxPriors) { // update of the maximum value
            maxPriors=(*Priors_PTRtmp);
        }
        if((*Priors_PTRtmp)<minPriors){ // update of the minimum value
            minPriors=(*Priors_PTRtmp);
        }
    }

    // reinitialisation of the pointer to the beginning of the Priors
    Priors_PTRtmp=Priors_PTR;
    // modification of the values
    if (minPriors==maxPriors) { // case (strange) where all the values are equal, then we put them all to 1 NO !!! Can correspond to weight spatially constant put to be afterwards adapted !!!
//        for (int i=0; i<numel; i++,Priors_PTRtmp++) {
//            (*Priors_PTRtmp)=1;
//        }
    }
    else {// rescaling of the values
        for (int i=0;i<numel;i++,Priors_PTRtmp++){ // Putting all the values between 0 and 1
            (*Priors_PTRtmp)=(*Priors_PTRtmp-minPriors)/(maxPriors-minPriors);
        }
    }
//    string filenameOut="/Users/Carole/Documents/PhD/Phantom/PriorsModified";
//    nifti_set_filenames(this->GetPriors(), filenameOut.c_str(), 0, 0);
//    nifti_image_write(this->GetPriors());
}

// Copies the content of NonNormRespInput into the NonNormResp of the tree on which it is applied
//void TreeEM::SetNonNormResp(float *NonNormRespInput){
//    int numelmasked=this->GetNumberMaskedElements();
//    if (NonNormRespInput==NULL) {
//        if (this->GetNonNormResp()!=NULL) {
//            delete []this->GetNonNormResp();
//        }
//        this->NonNormResp=NULL;
//        return;
//    }
//    if (this->GetNonNormResp()==NULL) {
//        this->NonNormResp=new float[this->GetNumberMaskedElements()];//{0};
//        for (int i=0; i<numelmasked; i++) {
//                    this->NonNormResp[i]=0;
//                }
//    }
//    float * NonNormResp_PTR=this->GetNonNormResp();
//    float * NonNormRespInput_PTR=NonNormRespInput;
//    for (int i=0; i<numelmasked; i++,NonNormResp_PTR++,NonNormRespInput_PTR++) {
//        *NonNormResp_PTR=*NonNormRespInput_PTR;
//    }
//    return;
//}

// Copies the content of NormRespInput into NormResp of the tree on which it is applied
void TreeEM::SetNormResp(float * NormRespInput){
    int numelmasked=this->GetNumberMaskedElements();
    if (NormRespInput==NULL) {
        if (this->GetNormResp()!=NULL) {
            delete []this->GetNormResp();
        }
        this->NormResp=NULL;
        return;
    }
    if (this->GetNormResp()==NULL) {
        this->NormResp=new float[this->GetNumberMaskedElements()];//{0};
        for (int i=0; i<numelmasked; i++) {
            this->NormResp[i]=0;
        }
    }
    float * NormResp_PTR=this->NormResp;
    float * NormRespInput_PTR=NormRespInput;

    for (int i=0; i<numelmasked; i++,NormResp_PTR++,NormRespInput_PTR++) {
        *NormResp_PTR=*NormRespInput_PTR;
    }
    return;
}

void TreeEM::SetHardSeg(int * HardSegInput){
    int numelmasked=this->GetNumberMaskedElements();
    if(!this->IsRoot()){
        cout<<"HardSeg only modified at root"<<endl;
        return;
    }
    if(HardSegInput==NULL){
        if(this->GetHardSegDirect()!=NULL){
            delete [] this->HardSeg;
        }
        this->HardSeg=NULL;
        return;
    }
    if(this->GetHardSegDirect()==NULL){
        this->HardSeg=new int[numelmasked];
        for(int i=0;i<numelmasked;i++){
            this->HardSeg[i]=0;
        }
    }
    int * HardSeg_PTR=this->HardSeg;
    int * HardSegInput_PTR=HardSegInput;
    for(int i=0;i<numelmasked;i++,HardSegInput_PTR++,HardSeg_PTR++){
        *HardSeg_PTR=*HardSegInput_PTR;
    }
    return;
}

// At this specific level of the tree reinitialisation of the NormResp with the values of the priors after making them usable in a probabilistic manner
void TreeEM::InitialiseNormRespWithPriors(){
    float * NormResp_PTR=this->GetNormResp();
    int numelmasked=this->GetNumberMaskedElements();
    if (NormResp_PTR==NULL) {
        float * NormRespToSet=new float[numelmasked];
        for (int i=0; i<numelmasked; i++) {
            NormRespToSet[i]=0;
        }
        this->SetNormResp( NormRespToSet);//{0});
    }
    if(!this->IsPriorsProbabilityType()){ // Conversion into probabilities if needed
        this->MakePriorsProbabilityType();
    }
    if(!this->ArePriorsNormalised()){
        this->NormalisePriors();
    }
    float * PriorsData_PTR=static_cast<float*>(this->GetPriors()->data);
    int * L2S_PTR=this->GetL2S();
    int numel=this->GetNumberElements();
    for (int i=0; i<numel; i++,PriorsData_PTR++,L2S_PTR++) {
        if ((*L2S_PTR)>=0) {
            (*NormResp_PTR)=(*PriorsData_PTR);
            NormResp_PTR++;
        }
    }
}
// Sets NonNormWeight of the tree on which it is applied to NonNormWeightInput
void TreeEM::SetNonNormWeight(float NonNormWeightInput){
    this->NonNormWeight=NonNormWeightInput;
}

// Sets NormWeight of the tree on which it is applied to NormWeightInput after checking that the value is between 0 and 1
void TreeEM::SetNormWeight(float NormWeightInput){
    if (NormWeightInput>=0&&NormWeightInput<=1) {
        this->NormWeight=NormWeightInput;
    }
    else{
        cout<<"the input value is not valid"<<endl;
    }
}

//void TreeEM::SetDistribution(float * DistributionInput){
//    int numelmasked=this->GetNumberMaskedElements();
//    if (DistributionInput==NULL) {
//        if (this->GetDistribution()!=NULL) {
//            delete []this->GetDistribution();
//        }
//        this->Distribution=NULL;
//        return;
//    }
//    if (this->GetDistribution()==NULL) {
//        this->Distribution=new float[this->GetNumberMaskedElements()];//{0};
//        for(int i=0;i<numelmasked;i++){
//            this->Distribution[i]=0;
//        }
//    }
//    float * Distribution_PTR=this->GetDistribution();
//    float * DistributionInput_PTR=DistributionInput;

//    //cout<<"Addresses distribution"<< Distribution_PTR << " "<< DistributionInput_PTR;
//    for (int i=0; i<numelmasked; i++,Distribution_PTR++,DistributionInput_PTR++) {
//        *Distribution_PTR=*DistributionInput_PTR;
//    }
//    return;

//}

// After checking the validity of such an update change the BF coeffs
void TreeEM::SetBFCoeffs(float * BFCoeffsInput){
    int numbmodal=this->GetNumberModalities();

    // Clearing current space for BFCoeffs if allocated
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    float * BFCoeffsToChange=this->GetBFCoeffs();
    if(this->GetParent()!=NULL){ // if it is not a root we do not change the BF coeffs
        cout<<"BF coeffs only set at the root"<<endl;
        return;
    }

    // Fill with BF CoeffsInput
    if(BFCoeffsInput!=NULL) { // it is only in the case we are at the root (checked before)
        if(BFCoeffsToChange==NULL){
            cout<<"Memory allocated for BF"<<endl;
            this->BFCoeffs=new float[numbmodal*numbBF];
        }
        BFCoeffsToChange=this->BFCoeffs;
        for (int l=0; l<numbmodal*numbBF;l++) {
            BFCoeffsToChange[l]=BFCoeffsInput[l];
        }
    }

}

//void TreeEM::SetBFCorrection(float * BFCorrectionInput){
//    int numbmodal=this->GetNumberModalities();
//    int numelmasked=this->GetNumberMaskedElements();
////    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
//    // Clearing current space for BFCoeffs if allocated

//    float * BFCorrectionToChange=this->BFCorrection;
//    if(this->GetParent()!=NULL){ // if it is not a root we do not change the BF coeffs
//        cout<<"BF correction only set at the root"<<endl;
//        return;
//    }

//    // Fill with BF CoeffsInput
//    if(BFCorrectionInput!=NULL) {
//        if(BFCorrectionToChange==NULL){
//        this->BFCorrection=new float[numelmasked*numbmodal];
//            BFCorrectionToChange=this->BFCorrection;
//        }
//        for (int i=0; i<numelmasked*numbmodal;i++) {
//            BFCorrectionToChange[i]=BFCorrectionInput[i];
//        }
//    }
//}


void TreeEM::SetDataBFCorrected(float * DataBFCorrectedInput){
    int numbmodal=this->GetNumberModalities();
    int numelmasked=this->GetNumberMaskedElements();

    // Clearing current space for BFCoeffs if allocated

    float * DataBFCorrectedToChange=this->DataBFCorrected;
    if(this->GetParent()!=NULL){ // if it is not a root we do not change the BF coeffs
        cout<<"Data correction only set at the root"<<endl;
        return;
    }
    // Fill with DataBFCorrectedInput
    if(DataBFCorrectedInput!=NULL) {
        if(DataBFCorrectedToChange==NULL){
            this->DataBFCorrected=new float[numelmasked*numbmodal];
            DataBFCorrectedToChange=this->DataBFCorrected;
        }
        for (int i=0; i<numelmasked*numbmodal;i++) {
            DataBFCorrectedToChange[i]=DataBFCorrectedInput[i];
        }
    }
}

void TreeEM::SetDPChildren(float* DPChildrenInput){
    if(!this->IsRoot()){
        cout<<"DP children only set at root"<<endl;
        return;
    }
    int numbchild=this->GetNumberChildren();
//    if(numbchild!=DPChildrenInput.size()){
//        cout<<"Incompatibility in size for DP"<<endl;
//        if(numbchild>DPChildrenInput.size()){
//            for(int i=DPChildrenInput.size();i<numbchild;i++){
//                DPChildrenInput.push_back(NULL);
//            }
//        }

//    }
    float * DPChildrenInitial=this->GetDPChildrenDirect();
    if(DPChildrenInitial!=NULL){
        delete [] DPChildrenInitial;
        DPChildrenInitial=NULL;
    }
    if(this->GetFlagDistClassInd()){
    float * DPChildrenToInput=new float[MaxSupport*numbchild];
    for(int c=0;c<numbchild;c++){
        for(int i=0;i<MaxSupport;i++){
            DPChildrenToInput[i+c*MaxSupport]=DPChildrenInput[i+c*MaxSupport];
        }
    }
    this->DPChildren=DPChildrenToInput;
    }
    else{
        int SizeHistogram=(int)pow_int(MaxSupport,numbchild);
        float *DPChildrenToInput=new float[SizeHistogram];
        for(int i=0;i<SizeHistogram;i++){
            DPChildrenToInput[i]=DPChildrenInput[i];
        }
        this->DPChildren=DPChildrenToInput;
    }



//        int DPsize=this->GetDPChildrenDirect().size();
//        vector<float *> DPChildrenInitial=this->GetDPChildrenDirect();
//        if(DPsize!=0){
//            for(int c=0;c<DPsize;c++){
//                if(DPChildrenInitial[c]!=NULL){
//                    delete [] DPChildrenInitial[c];
//                    DPChildrenInitial[c]=NULL;
//                }
//            }
//            DPChildrenInitial.clear();
//        }
//        this->DPChildren=DPChildrenInput;
}

void TreeEM::SetFlagDistClassInd(bool flag_DistClassInd){
    if(!this->IsRoot()){
        this->GetParent()->SetFlagDistClassInd(flag_DistClassInd);
    }
    else{
        this->FlagDistClassInd=flag_DistClassInd;
    }
}

void TreeEM::SetFlagOutliers(int flag_Outliers){
    if(!this->IsRoot()){
        this->GetParent()->SetFlagOutliers(flag_Outliers);
    }
    else{
        this->FlagOutliers=flag_Outliers;
    }
}

/********************* EM RELATED METHODS ***************************/

// Update all the parameters of the model in a recursive way
void TreeEM::UpdateParameters(){
    int TypeTreeComponent=this->WhatTypeOfTreeComponent();
    if(TypeTreeComponent==LEAF){// Means it is a leaf
        /* For the moment only Gaussian distribution as simple distribution so no other case than the default one in the switch*/
        switch (this->GetDistributionType()) {
        case 2: // case of uniform distribution with outliers : no updating
            break;
        default:
            this->UpdateGaussianParameters();
            //cout<<"Gaussian parameters updated";
            break;
        }
    }
    else{ // Recursive part
        for (int c=0; c<this->GetChildren().size(); c++) {
            //cout<<"Updating parameters of child "<<c<<endl;
            this->GetChild(c)->UpdateParameters();
        }
    }
    return;
}

void TreeEM::UpdateGaussianParameters(){
    if (this->GetDistributionType()!=1){
        cout<<"It is not a Gaussian distribution"<<endl;
        return;
    }
    this->UpdateGaussianMean();
    this->UpdateGaussianVariance();
    return;
}

// Update the Adapted Priors : only needed on nodes that have direct priors and if AtlasWeight is non zero
void TreeEM::UpdatePriorsAdapted(SEG_PARAMETERS *segment_param){
    float Weight=segment_param->AtlasWeight;
    if(Weight==0){
        return;
    }
//    int numbchild=this->GetNumberChildren();
    vector<TreeEM*> NodePriorsVector=this->GetPriorsNodeVector();
    int numbclasses=NodePriorsVector.size();
    if(numbclasses >0){
        for(int c=0;c<numbclasses;c++){
            float * PriorsAdaptedToUpdate=NodePriorsVector[c]->CopyPriorsAdapted();
            float * UpdatedPriors = NodePriorsVector[c]->AdaptPriors(segment_param);
            NodePriorsVector[c]->SetPriorsAdapted(UpdatedPriors);
            delete [] PriorsAdaptedToUpdate;
            PriorsAdaptedToUpdate=NULL;
            delete [] UpdatedPriors;
            UpdatedPriors=NULL;

//                // Check for number of zero values in PriorsAdapted :
//                float * PriorsAdapted_PTR=this->GetChild(c)->GetPriorsAdapted();
//                int * L2S_PTR=this->GetL2S();
//                int numel = this->GetNumberElements();
//                int CountPriorsAdaptedZero=0;
//                for(int i=0;i<numel;i++,L2S_PTR++,PriorsAdapted_PTR++){
//                    if(*L2S_PTR>=0){
//                        if(*PriorsAdapted_PTR==0){
//                            CountPriorsAdaptedZero++;
//                        }
//                    }
//                }
//                cout<<"Number of zero values in new priors masked is "<<CountPriorsAdaptedZero<<" for child "<<c<<endl;
        }
        this->NormalisePriorsAdapted();

    }
}
//    }
//    if(numbchild >0){
//        if(this->GetChild(0)->GetPriorsAdaptedDirect()!=NULL){
//            for(int c=0;c<numbchild;c++){
//                float * PriorsAdaptedToUpdate=this->GetChild(c)->CopyPriorsAdapted();
//                float * UpdatedPriors = this->GetChild(c)->AdaptPriors(segment_param);
//                this->GetChild(c)->SetPriorsAdapted(UpdatedPriors);
//                delete [] PriorsAdaptedToUpdate;
//                PriorsAdaptedToUpdate=NULL;
//                delete [] UpdatedPriors;
//                UpdatedPriors=NULL;

////                // Check for number of zero values in PriorsAdapted :
////                float * PriorsAdapted_PTR=this->GetChild(c)->GetPriorsAdapted();
////                int * L2S_PTR=this->GetL2S();
////                int numel = this->GetNumberElements();
////                int CountPriorsAdaptedZero=0;
////                for(int i=0;i<numel;i++,L2S_PTR++,PriorsAdapted_PTR++){
////                    if(*L2S_PTR>=0){
////                        if(*PriorsAdapted_PTR==0){
////                            CountPriorsAdaptedZero++;
////                        }
////                    }
////                }
////                cout<<"Number of zero values in new priors masked is "<<CountPriorsAdaptedZero<<" for child "<<c<<endl;
//            }
//            this->NormalisePriorsAdapted();

//        }
//    }

//    // RecursivePart
//    for(int c=0;c<numbchild;c++){
//        this->GetChild(c)->UpdatePriorsAdapted(segment_param);
//    }
//}


void TreeEM::UpdateGaussianMean(){
//    int numel=this->GetNumberElements(); // number of elements in the image
    int numelmasked=this->GetNumberMaskedElements();
    int numbmodal=this->GetNumberModalities(); // number of modalities to consider
    //float * PointerToDataBegin=static_cast<float *>(this->GetDataImage()->data);
    //    float * PointerToDataBegin=this->MakeDataBFCorrected();
    float * PointerToDataBegin=this->GetDataBFCorrected();
    for (int m=0;m<numbmodal;m++){
        float * PointerToImageBegin_PTR=&PointerToDataBegin[m*numelmasked]; // Points the data to the beginning of each modality image
        float * PointerToNormRespBegin_PTR=this->GetNormResp(); // Points to the beginning of the NormalizedResponsabilities
        double meantmp=0; // Initialisation of the numerator in the update of the mean
        double sumResp=0; // Initialisation of the denominator in the update of the mean
        float meantmpMarc=0.f; // Initialisation of the numerator in the update of the mean
        float sumRespMarc=0.f; // Initialisation of the denominator in the update of the mean
        //        int * L2S_PTR=this->GetL2S();
        int CountMoreThanOne=0;
        for(int i=0;i<numelmasked;i++,PointerToImageBegin_PTR++,PointerToNormRespBegin_PTR++){
            /* Each time the considered voxel is active (<=> L2S contains an index), update the numerator and the denominator*/
            meantmp+=(double)(*PointerToImageBegin_PTR)*(*PointerToNormRespBegin_PTR);
            meantmpMarc+=(*PointerToImageBegin_PTR)*(*PointerToNormRespBegin_PTR);
            if (*PointerToImageBegin_PTR>1) {
                CountMoreThanOne++;
                //                    cout<<"Pb in normalisation image";
            }
            sumResp+=(double)(*PointerToNormRespBegin_PTR);
            sumRespMarc+=(*PointerToNormRespBegin_PTR);
        }
        if(CountMoreThanOne>0){
            //            cout<<"Transformed values more than one is "<<CountMoreThanOne<<endl;
        }
        if (sumResp!=0) {
            this->GetMean()[m]=meantmp/sumResp;
        }
        else{
            this->GetMean()[m]=meantmp/this->GetNumberMaskedElements();
        }
//        printf("HERE We are float %f / %f = %f\n", meantmpMarc, sumRespMarc,meantmpMarc/sumRespMarc);
//        printf("HERE We are double %g / %g = %g\n", meantmp, sumResp,meantmp/sumResp);
//        exit(1);

        //cout<<"Mean is now "<< this->GetMean()[m];
    }
    // Clearing memory
    //    if (PointerToDataBegin!=NULL) {
    //        delete [] PointerToDataBegin;
    //        PointerToDataBegin=NULL;
    //    }
    return;
}

void TreeEM::UpdateGaussianVariance(){

//    int numel=this->GetNumberElements();
    int numbmodal=this->GetNumberModalities();
    int numelmasked=this->GetNumberMaskedElements();
    PrecisionTYPE sumResp=0; // Initialisation of the denominator
    //float * PointerToDataBegin = static_cast<float *>(this->GetDataImage()->data);
    //    float * PointerToDataBegin=this->MakeDataBFCorrected();
    float * PointerToDataBegin=this->GetDataBFCorrected();
    float * PointerToDataBegin_PTR1=PointerToDataBegin;
    float * PointerToDataBegin_PTR2=PointerToDataBegin;
    float * NormalisedResponsabilities_PTR=this->GetNormResp();
    //    int * L2S_PTR=this->GetL2S();

    // Calculation of the denominator (sum over the active voxels of the normalised responsabilities)
    for (int i=0; i<numelmasked; i++,NormalisedResponsabilities_PTR++) {
        sumResp+=(PrecisionTYPE)(*NormalisedResponsabilities_PTR);
    }

    float * VarianceToUpdateB=this->GetVariance();
    float VarianceToUpdate[MaxNumbModal*MaxNumbModal];
    for(int m1=0;m1<MaxNumbModal;m1++){
        for(int m2=0;m2<MaxNumbModal;m2++){
            if(m1<numbmodal&&m2<numbmodal){
                VarianceToUpdate[m1+m2*MaxNumbModal]=VarianceToUpdateB[m1+m2*numbmodal];
            }
            else{
                VarianceToUpdate[m1+m2*MaxNumbModal]=0;
            }
        }
    }
    float * MeanToUseB=this->GetMean();
    float MeanToUse[MaxNumbModal];
    for(int m=0;m<MaxNumbModal;m++){
        if(m<numbmodal){
            MeanToUse[m]=MeanToUseB[m];
        }
    }
    for(int m1=0;m1<numbmodal;m1++){
        // First data pointer to the beginning of the modality m1 considered
        for(int m2=0;m2<numbmodal;m2++){
            PointerToDataBegin_PTR1=&PointerToDataBegin[m1*numelmasked];
            PointerToDataBegin_PTR2=&PointerToDataBegin[m2*numelmasked]; // Second Data pointer to the beginning of the modality m2 considered
            NormalisedResponsabilities_PTR=this->GetNormResp(); // Reinitialisation of the responsabilities pointer to the beginning
            VarianceToUpdate[m1+m2*MaxNumbModal]=0;
            PrecisionTYPE VarianceToUpdate_tmp=0;
            for(int i=0;i<numelmasked;i++,PointerToDataBegin_PTR1++,PointerToDataBegin_PTR2++,NormalisedResponsabilities_PTR++){
                // Update of the numerator of the Variance Calculation only if in the case of an active voxel
                VarianceToUpdate_tmp+=(PrecisionTYPE)(*NormalisedResponsabilities_PTR)*((*PointerToDataBegin_PTR1)-MeanToUse[m1])*((*PointerToDataBegin_PTR2)-MeanToUse[m2]);
            }
            if (sumResp !=0) {
                VarianceToUpdate[m1+m2*MaxNumbModal]=(float)VarianceToUpdate_tmp/sumResp;
                if (m1==m2) {
                    VarianceToUpdate[m1+m2*MaxNumbModal]=VarianceToUpdate[m1+m2*MaxNumbModal]<=1E-6?1E-6:VarianceToUpdate[m1+m2*MaxNumbModal]; // in order to avoid going to 0 if too sharp distribution but not changing non diagonal of variance
                }


            }
            else{
                VarianceToUpdate[m1+m2*MaxNumbModal]=VarianceToUpdate_tmp/numelmasked;// this->GetNumberMaskedElements();
            }
            // Use of the symmetry property of the Variance matrix
            VarianceToUpdate[m2+m1*MaxNumbModal]=VarianceToUpdate[m1+m2*MaxNumbModal];
        }
    }
    for(int m1=0;m1<numbmodal;m1++){
        for(int m2=0;m2<numbmodal;m2++){
            VarianceToUpdateB[m1+m2*numbmodal]=VarianceToUpdate[m1+m2*MaxNumbModal];
        }
    }
    // Clearing memory
    //    if (PointerToDataBegin!=NULL) {
    //        delete []  PointerToDataBegin;
    //        PointerToDataBegin=NULL;
    //    }
    return;
}


// Update all the parameters of the model in a recursive way
void TreeEM::UpdateParametersCEM(int Child){
    int TypeTreeComponent=this->WhatTypeOfTreeComponent();
    if(TypeTreeComponent==LEAF){// Means it is a leaf
        /* For the moment only Gaussian distribution as simple distribution so no other case than the default one in the switch*/
        switch (this->GetDistributionType()) {
        default:
            this->UpdateGaussianParametersCEM(Child);
            //cout<<"Gaussian parameters updated";
            break;
        }
    }
    else{ // Recursive part
        for (int c=0; c<this->GetChildren().size(); c++) {
            //cout<<"Updating parameters of child "<<c<<endl;
            this->GetChild(c)->UpdateParametersCEM(Child);
        }
    }
    return;
}

void TreeEM::UpdateGaussianParametersCEM(int c){
    if (this->GetDistributionType()!=1){
        cout<<"It is not a Gaussian distribution"<<endl;
        return;
    }
    this->UpdateGaussianMeanCEM(c);
    this->UpdateGaussianVarianceCEM(c);
    return;
}


void TreeEM::UpdateGaussianMeanCEM(int c){
//    int numel=this->GetNumberElements(); // number of elements in the image
    int numelmasked=this->GetNumberMaskedElements();
    int numbmodal=this->GetNumberModalities(); // number of modalities to consider
    //float * PointerToDataBegin=static_cast<float *>(this->GetDataImage()->data);
    //    float * PointerToDataBegin=this->MakeDataBFCorrected();
    float * PointerToDataBegin=this->GetDataBFCorrected();
    int * HardSegUsed=this->GetHardSeg();
    for (int m=0;m<numbmodal;m++){
        float * PointerToImageBegin_PTR=&PointerToDataBegin[m*numelmasked]; // Points the data to the beginning of each modality image
        float * PointerToNormRespBegin_PTR=this->GetNormResp(); // Points to the beginning of the NormalizedResponsabilities
        int * HardSeg_PTR=this->GetHardSeg();
        double meantmp=0; // Initialisation of the numerator in the update of the mean
        double sumResp=0; // Initialisation of the denominator in the update of the mean
        float meantmpMarc=0.f; // Initialisation of the numerator in the update of the mean
        float sumRespMarc=0.f; // Initialisation of the denominator in the update of the mean
        //        int * L2S_PTR=this->GetL2S();
        int CountMoreThanOne=0;
        for(int i=0;i<numelmasked;i++,PointerToImageBegin_PTR++,PointerToNormRespBegin_PTR++,HardSeg_PTR++){
            /* Each time the considered voxel is active (<=> L2S contains an index), update the numerator and the denominator*/
            if(*HardSeg_PTR==c){
            meantmp+=(double)(*PointerToImageBegin_PTR)*(*PointerToNormRespBegin_PTR);
            meantmpMarc+=(*PointerToImageBegin_PTR)*(*PointerToNormRespBegin_PTR);
            if (*PointerToImageBegin_PTR>1) {
                CountMoreThanOne++;
                //                    cout<<"Pb in normalisation image";
            }
            sumResp+=(double)(*PointerToNormRespBegin_PTR);
            sumRespMarc+=(*PointerToNormRespBegin_PTR);
            }
        }
        if(CountMoreThanOne>0){
            //            cout<<"Transformed values more than one is "<<CountMoreThanOne<<endl;
        }
        if (sumResp!=0) {
            this->GetMean()[m]=meantmp/sumResp;
        }
        else{
            this->GetMean()[m]=meantmp/this->GetNumberMaskedElements();
        }
//        printf("HERE We are float %f / %f = %f\n", meantmpMarc, sumRespMarc,meantmpMarc/sumRespMarc);
//        printf("HERE We are double %g / %g = %g\n", meantmp, sumResp,meantmp/sumResp);
//        exit(1);

        //cout<<"Mean is now "<< this->GetMean()[m];
    }
    // Clearing memory
    //    if (PointerToDataBegin!=NULL) {
    //        delete [] PointerToDataBegin;
    //        PointerToDataBegin=NULL;
    //    }
    return;
}

void TreeEM::UpdateGaussianVarianceCEM(int c){

//    int numel=this->GetNumberElements();
    int numbmodal=this->GetNumberModalities();
    int numelmasked=this->GetNumberMaskedElements();
    PrecisionTYPE sumResp=0; // Initialisation of the denominator
    //float * PointerToDataBegin = static_cast<float *>(this->GetDataImage()->data);
    //    float * PointerToDataBegin=this->MakeDataBFCorrected();
    float * PointerToDataBegin=this->GetDataBFCorrected();
    float * PointerToDataBegin_PTR1=PointerToDataBegin;
    float * PointerToDataBegin_PTR2=PointerToDataBegin;
    float * NormalisedResponsabilities_PTR=this->GetNormResp();
    int * HardSeg_PTR=this->GetHardSeg();
    //    int * L2S_PTR=this->GetL2S();

    // Calculation of the denominator (sum over the active voxels of the normalised responsabilities)
    for (int i=0; i<numelmasked; i++,NormalisedResponsabilities_PTR++,HardSeg_PTR++) {
        if(*HardSeg_PTR==c){
        sumResp+=(PrecisionTYPE)(*NormalisedResponsabilities_PTR);
        }
    }

    float * VarianceToUpdateB=this->GetVariance();
    float VarianceToUpdate[MaxNumbModal*MaxNumbModal];
    for(int m1=0;m1<MaxNumbModal;m1++){
        for(int m2=0;m2<MaxNumbModal;m2++){
            if(m1<numbmodal&&m2<numbmodal){
                VarianceToUpdate[m1+m2*MaxNumbModal]=VarianceToUpdateB[m1+m2*numbmodal];
            }
            else{
                VarianceToUpdate[m1+m2*MaxNumbModal]=0;
            }
        }
    }
    float * MeanToUseB=this->GetMean();
    float MeanToUse[MaxNumbModal];
    for(int m=0;m<MaxNumbModal;m++){
        if(m<numbmodal){
            MeanToUse[m]=MeanToUseB[m];
        }
    }
    for(int m1=0;m1<numbmodal;m1++){
        // First data pointer to the beginning of the modality m1 considered
        for(int m2=0;m2<numbmodal;m2++){
            PointerToDataBegin_PTR1=&PointerToDataBegin[m1*numelmasked];
            PointerToDataBegin_PTR2=&PointerToDataBegin[m2*numelmasked]; // Second Data pointer to the beginning of the modality m2 considered
            NormalisedResponsabilities_PTR=this->GetNormResp(); // Reinitialisation of the responsabilities pointer to the beginning
            HardSeg_PTR=this->GetHardSeg();
            VarianceToUpdate[m1+m2*MaxNumbModal]=0;
            PrecisionTYPE VarianceToUpdate_tmp=0;
            for(int i=0;i<numelmasked;i++,PointerToDataBegin_PTR1++,PointerToDataBegin_PTR2++,NormalisedResponsabilities_PTR++,HardSeg_PTR++){
                // Update of the numerator of the Variance Calculation only if in the case of an active voxel
                if(*HardSeg_PTR==c){
                VarianceToUpdate_tmp+=(PrecisionTYPE)(*NormalisedResponsabilities_PTR)*((*PointerToDataBegin_PTR1)-MeanToUse[m1])*((*PointerToDataBegin_PTR2)-MeanToUse[m2]);
                }
            }
            if (sumResp !=0) {
                VarianceToUpdate[m1+m2*MaxNumbModal]=(float)VarianceToUpdate_tmp/sumResp;
                if (m1==m2) {
                    VarianceToUpdate[m1+m2*MaxNumbModal]=VarianceToUpdate[m1+m2*MaxNumbModal]<=1E-6?1E-6:VarianceToUpdate[m1+m2*MaxNumbModal]; // in order to avoid going to 0 if too sharp distribution but not changing non diagonal of variance
                }


            }
            else{
                VarianceToUpdate[m1+m2*MaxNumbModal]=VarianceToUpdate_tmp/numelmasked;// this->GetNumberMaskedElements();
            }
            // Use of the symmetry property of the Variance matrix
            VarianceToUpdate[m2+m1*MaxNumbModal]=VarianceToUpdate[m1+m2*MaxNumbModal];
        }
    }
    for(int m1=0;m1<numbmodal;m1++){
        for(int m2=0;m2<numbmodal;m2++){
            VarianceToUpdateB[m1+m2*numbmodal]=VarianceToUpdate[m1+m2*MaxNumbModal];
        }
    }
    // Clearing memory
    //    if (PointerToDataBegin!=NULL) {
    //        delete []  PointerToDataBegin;
    //        PointerToDataBegin=NULL;
    //    }
    return;
}




// Update the obtention of the bias field coefficients
void TreeEM::UpdateBFCoeffs(){
    //    float * BFCoeffsToUpdate=this->MakeFinalBFCoeffsChildren();
    float * BFCoeffsToUpdate=this->MakeBFCoeffsDirect();
    this->SetBFCoeffs(BFCoeffsToUpdate);
    // Clearing memory not needed anymore.
    if (BFCoeffsToUpdate!=NULL) {
        delete [] BFCoeffsToUpdate;
        BFCoeffsToUpdate=NULL;
    }
}

//// Update the Bias Field Correction
//void TreeEM::UpdateBFCorrection(){
//    float * BFCorrectionToUpdate=this->MakeBFCorrection();
//    this->SetBFCorrection(BFCorrectionToUpdate);
//    // Clearing memory not needed anymore
//    if(BFCorrectionToUpdate!=NULL){
//        delete [] BFCorrectionToUpdate;
//        BFCorrectionToUpdate=NULL;
//    }
//}

// Update the corrected Data
void TreeEM::UpdateDataBFCorrected(){
    float * DataBFCorrectedToUpdate=this->MakeDataBFCorrected();
    this->SetDataBFCorrected(DataBFCorrectedToUpdate);
    // Clearing memory not needed anymore
    if(DataBFCorrectedToUpdate!=NULL){
        delete [] DataBFCorrectedToUpdate;
        DataBFCorrectedToUpdate=NULL;
    }
}


// Recursively obtain the NonNormWeight (sum over the voxels for one class of the normalised responsabilities)
void TreeEM::UpdateNonNormWeights(){
    int numelmasked=this->GetNumberMaskedElements();
    int numbchild=this->GetNumberChildren();
    this->SetNonNormWeight(0);
    float * NormResp_PTR=this->GetNormResp();
    for(int i=0;i<numelmasked;i++,NormResp_PTR++){
        this->SetNonNormWeight((PrecisionTYPE)this->GetNonNormWeight()+(*NormResp_PTR));
    }
    for(int c=0;c<numbchild;c++){
        this->GetChild(c)->UpdateNonNormWeights();
    }
    return;
}

// Recursively obtain the NonNormWeight (sum over the voxels for one class of the normalised responsabilities)
void TreeEM::UpdateNonNormWeightsCEM(int c){
    int numelmasked=this->GetNumberMaskedElements();
    int numbchild=this->GetNumberChildren();
    this->SetNonNormWeight(0);
    float * NormResp_PTR=this->GetNormResp();
    int * HardSeg_PTR=this->GetHardSeg();
    for(int i=0;i<numelmasked;i++,NormResp_PTR++,HardSeg_PTR++){
        if(*HardSeg_PTR==c){
        this->SetNonNormWeight((PrecisionTYPE)this->GetNonNormWeight()+(*NormResp_PTR));
        }
    }
    for(int c=0;c<numbchild;c++){
        this->GetChild(c)->UpdateNonNormWeights();
    }
    return;
}


// Recursively Update the normalised weights when the non normalised ones are obtained
void TreeEM::UpdateNormWeights(){
    // First determination of the normalisation factor for the children
    PrecisionTYPE SumNonNormWeights=0;
    int numbchild=this->GetNumberChildren();
    for (int c=0;c<numbchild; c++) {
        SumNonNormWeights+=(PrecisionTYPE)this->GetChild(c)->GetNonNormWeight();
    }
    // Then normalisation of the not normalised weights for the children
    for (int c=0; c<numbchild; c++) {
        if (SumNonNormWeights>0) {
            this->GetChild(c)->SetNormWeight((PrecisionTYPE)this->GetChild(c)->GetNonNormWeight()/SumNonNormWeights);
        }
        else{// in case sum of the weights is 0
            this->GetChild(c)->SetNormWeight(1.0/this->GetNumberChildren());
        }
        //cout<< "Norm weight of child "<<c<<" is of value "<<this->GetChild(c)->GetNormWeight()<<endl;
        this->GetChild(c)->UpdateNormWeights(); // Recursive part
    }

    return;
}


void TreeEM::MakeNonNormWeightedSum(float * NonNormSum, SEG_PARAMETERS * segment_param){
    int numbchild=this->GetNumberChildren();
    int numelmasked=this->GetNumberMaskedElements();
    int numel=this->GetNumberElements();
    if(numbchild==0){// case where it is a leaf
        float* DistSum_tmp=NULL;
        float * DistSum_tmp_PTR=NULL;
        switch (this->GetDistributionType()){
        case 2:{
             DistSum_tmp=this->MakeUniformDistribution();
             DistSum_tmp_PTR=DistSum_tmp;
            // Supposedly no priors
            float NormWeightUsed=this->GetPartNormWeight();
            for(int i=0;i<numelmasked;i++,DistSum_tmp_PTR++){
                (*DistSum_tmp_PTR)*=NormWeightUsed;
            }
            break;
        }

        default : {
            int CountZeroValuesNonNorm=0;
             DistSum_tmp=this->MakeGaussianDistribution();
             DistSum_tmp_PTR=DistSum_tmp;
            if(this->GetPriorsDirect()==NULL){// case we have to multiply by the normalised weight until
                float NormWeightUsed=this->GetPartNormWeight();
                TreeEM* NodePriors=this->FindGeneralClassPriors();
                if(NodePriors!=NULL){ // case where there are some priors in the tree but not at the first level (outliers model used there)
                float ComplementNormWeightUsed=NodePriors->GetPartNormWeight();
                NormWeightUsed*=ComplementNormWeightUsed;
                }
                for(int i=0;i<numelmasked;i++,DistSum_tmp_PTR++){
                    (*DistSum_tmp_PTR)*=NormWeightUsed;
                    if(*DistSum_tmp_PTR==0){
                        CountZeroValuesNonNorm++;
                    }
                }
//                cout<< "Zero values after weighting is "<< CountZeroValuesNonNorm<<endl;
            }
        break;
        }
        }

            if(this->GetPriors()!=NULL){// case where we are using atlases
                int CountZeroValuesNonNorm=0;
                DistSum_tmp_PTR=DistSum_tmp;
                int *L2S_PTR=this->GetL2S();
//                float *Priors_PTR=new float[numel];
//                this->AdaptPriors(segment_param,Priors_PTR);
                float * Priors_PTR=this->GetPriorsAdapted();
                CountZeroValuesNonNorm=0;
//                float * Priors_PTR=static_cast<float*>(this->GetPriors()->data);
                for(int i=0;i<numel;i++,L2S_PTR++){
                    if(*L2S_PTR>=0){
                        (*DistSum_tmp_PTR)*=Priors_PTR[i];
                        if(*DistSum_tmp_PTR==0){
                            CountZeroValuesNonNorm++;
                        }
                        DistSum_tmp_PTR++;
                    }
                }
//                cout<< "Zero Values after priors is "<< CountZeroValuesNonNorm<<endl;
//                delete[] Priors_PTR;
//                Priors_PTR=NULL;
            }

            if(this->GetMRF()!=NULL && !this->IsMRFZero()){ // Case where we have an MRF to apply
                DistSum_tmp_PTR=DistSum_tmp;
                float * MRFToUse=this->GetMRF();
                int CountZeroValuesNonNorm=0;
                for(int i=0;i<numelmasked;i++,DistSum_tmp_PTR++){
                    (*DistSum_tmp_PTR)*=MRFToUse[i];
                    if(*DistSum_tmp_PTR==0){
                        CountZeroValuesNonNorm++;
                    }
                }
//                cout<< "Zero value after MRF is "<<CountZeroValuesNonNorm<<endl;
            }
//            SaveTmpResultMasked(DistSum_tmp,"/Users/Carole/Documents/PhD/NormRespMRFTest.nii.gz");
            // Copying the final result in the destined float array
            float * NonNormSum_PTR=NonNormSum;
            DistSum_tmp_PTR=DistSum_tmp;
            int CountZeroValuesNonNorm=0;
            for(int i=0;i<numelmasked;i++,NonNormSum_PTR++,DistSum_tmp_PTR++){
                *NonNormSum_PTR=*DistSum_tmp_PTR;
                if(*NonNormSum_PTR==0){
                    CountZeroValuesNonNorm++;
                }
            }
//            cout<< "NonNorm zero values are : "<<CountZeroValuesNonNorm<<endl;
//            SaveTmpResultMasked(NonNormSum,"/Users/Carole/Documents/PhD/NormRespMRFTest.nii.gz");
            // Clearing memory not needed anymore
            delete[] DistSum_tmp;

            DistSum_tmp=NULL;

    }

    else {
        for(int c=0;c<numbchild;c++){
            // Initialisation of temporary result
            float * DistSum_tmp=new float[numelmasked];
            for(int i=0;i<numelmasked;i++){
                DistSum_tmp[i]=0;
            }
            this->GetChild(c)->MakeNonNormWeightedSum(DistSum_tmp,segment_param); // recursive part
            float * NonNormSum_PTR=NonNormSum;
            float * DistSum_tmp_PTR=DistSum_tmp;
            float NormWeightUsed=1;
            //            if(!this->GetChild(c)->IsLeaf()&&this->GetChild(c)->GetPriorsDirect()==NULL){
            //                NormWeightUsed=this->GetChild(c)->GetNormWeight();
            //            }
            for(int i=0;i<numelmasked;i++,NonNormSum_PTR++,DistSum_tmp_PTR++){
                (*NonNormSum_PTR)+=*DistSum_tmp_PTR*NormWeightUsed;
            }
            delete [] DistSum_tmp;
            DistSum_tmp=NULL;
        }

    }

    return;
}

void TreeEM::InitialiseNormResp(){
    int numelmasked=this->GetNumberMaskedElements();
    int numel=this->GetNumberElements();
    float * NormRespInit=NULL;
    float * NormRespInit_PTR=NULL;
    if(this->IsLeaf()){
//        NormRespInit=this->GetNormResp();
//        if(NormRespInit!=NULL){
//            delete [] NormRespInit;
//            NormRespInit=NULL;
//        }
            NormRespInit=new float[numelmasked];
            for(int i=0;i<numelmasked;i++){
                NormRespInit[i]=0;
            }
            NormRespInit_PTR=NormRespInit;
            if(this->GetPriorsAdaptedDirect()!=NULL){//case where there are some priors to put as initial values for the normalised responsabilities
                float * Priors_PTR=this->GetPriorsAdapted();
                int * L2S_PTR=this->GetL2S();
                for(int i=0;i<numel;i++,Priors_PTR++,L2S_PTR++){
                    if(*L2S_PTR>=0){
                    *NormRespInit_PTR=*Priors_PTR;
                        NormRespInit_PTR++;
                    }
                }
        }
            else{
                float NormWeightUsed=this->GetNormWeight();
                for(int i=0;i<numelmasked;i++){
                    NormRespInit[i]=NormWeightUsed;
                }
            }
            float NormWeightNormalisation=this->GetPartNormWeightAbovePriors();
            for(int i=0;i<numelmasked;i++){
                NormRespInit[i]*=NormWeightNormalisation;
            }
            this->SetNormResp(NormRespInit);
            delete[] NormRespInit;
            NormRespInit=NULL;
    }
    else{
        int numbchild=this->GetNumberChildren();
//        NormRespInit=this->GetNormResp();
//        if(NormRespInit!=NULL){
//            delete [] NormRespInit;
//            NormRespInit=NULL;
//        }
            NormRespInit=new float[numelmasked];
            for(int i=0;i<numelmasked;i++){
                NormRespInit[i]=0;
            }
        for(int c=0;c<numbchild;c++){
            this->GetChild(c)->InitialiseNormResp(); // Recursive part
            float * NormRespToAdd_PTR=this->GetChild(c)->GetNormResp();
            for(int i=0;i<numelmasked;i++,NormRespToAdd_PTR++){
                NormRespInit[i]+=*NormRespToAdd_PTR;
            }
        }
//        float NormWeightToUse=this->GetNormWeight();
//        for(int i=0;i<numelmasked;i++){
//            NormRespInit[i]*=NormWeightToUse;
//        }
        this->SetNormResp(NormRespInit);
        delete [] NormRespInit;
        NormRespInit=NULL;
    }
    return;
}

void TreeEM::MakeNonNormWeightedSumCEM(float * NonNormSum){
    int numbchild=this->GetNumberChildren();
    int numelmasked=this->GetNumberMaskedElements();
//    int numel=this->GetNumberElements();
    if(numbchild==0){// case where it is a leaf
        switch (this->GetDistributionType()){
        default : {
            float * DistSum_tmp=this->MakeGaussianDistribution();
            float * DistSum_tmp_PTR=DistSum_tmp;
            if(this->GetPriorsDirect()==NULL){// case we have to multiply by the normalised weight
                float NormWeightUsed=this->GetPartNormWeight();
                for(int i=0;i<numelmasked;i++,DistSum_tmp_PTR++){
                    (*DistSum_tmp_PTR)*=NormWeightUsed;
                }
            }
//            if(this->GetPriors()!=NULL){// case where we are using atlases
//                DistSum_tmp_PTR=DistSum_tmp;
//                int *L2S_PTR=this->GetL2S();
//                float * Priors_PTR=static_cast<float*>(this->GetPriors()->data);
//                for(int i=0;i<numel;i++,L2S_PTR++,Priors_PTR++){
//                    if(*L2S_PTR>=0){
//                        (*DistSum_tmp_PTR)*=*Priors_PTR;
//                        DistSum_tmp_PTR++;
//                    }
//                }
//            }
            // Copying the final result in the destined float array
            float * NonNormSum_PTR=NonNormSum;
            DistSum_tmp_PTR=DistSum_tmp;
            for(int i=0;i<numelmasked;i++,NonNormSum_PTR++,DistSum_tmp_PTR++){
                *NonNormSum_PTR=*DistSum_tmp_PTR;
            }
            // Clearing memory not needed anymore
            delete[] DistSum_tmp;
            DistSum_tmp=NULL;

        }
        }
    }

    else {
        for(int c=0;c<numbchild;c++){
            // Initialisation of temporary result
            float * DistSum_tmp=new float[numelmasked];
            for(int i=0;i<numelmasked;i++){
                DistSum_tmp[i]=0;
            }
            this->GetChild(c)->MakeNonNormWeightedSumCEM(DistSum_tmp); // recursive part
            float * NonNormSum_PTR=NonNormSum;
            float * DistSum_tmp_PTR=DistSum_tmp;
            float NormWeightUsed=1;
            //            if(!this->GetChild(c)->IsLeaf()&&this->GetChild(c)->GetPriorsDirect()==NULL){
            //                NormWeightUsed=this->GetChild(c)->GetNormWeight();
            //            }
            for(int i=0;i<numelmasked;i++,NonNormSum_PTR++,DistSum_tmp_PTR++){
                (*NonNormSum_PTR)+=*DistSum_tmp_PTR*NormWeightUsed;
            }
            delete [] DistSum_tmp;
            DistSum_tmp=NULL;
        }

    }

    return;
}


void TreeEM::MakeWeightedDist(float * WeightedDist){
    int numbchild=this->GetNumberChildren();
    int numelmasked=this->GetNumberMaskedElements();
    //    int numel=this->GetNumberElements();
    if(numbchild==0){// case where it is a leaf
        switch (this->GetDistributionType()){
        case 2: {// case of uniform distribution : appears when outlier model is used
            float * DistSum_tmp=this->MakeUniformDistribution();
            float * WeightedDist_PTR=WeightedDist;
            float * DistSum_tmp_PTR=DistSum_tmp;
            for(int i=0;i<numelmasked;i++,WeightedDist_PTR++,DistSum_tmp_PTR++){
                *WeightedDist_PTR=*DistSum_tmp_PTR;
            }
            // Clearing memory not needed anymore
            delete[] DistSum_tmp;
            DistSum_tmp=NULL;
            break;
            }

        default : {
            float * DistSum_tmp=this->MakeGaussianDistribution();
            float * DistSum_tmp_PTR=DistSum_tmp;

            // Copying the final result in the destined float array
            float * WeightedDist_PTR=WeightedDist;
            DistSum_tmp_PTR=DistSum_tmp;
            for(int i=0;i<numelmasked;i++,WeightedDist_PTR++,DistSum_tmp_PTR++){
                *WeightedDist_PTR=*DistSum_tmp_PTR;
            }
            // Clearing memory not needed anymore
            delete[] DistSum_tmp;
            DistSum_tmp=NULL;
            break;
        }

        }
    }
    else {
        for(int c=0;c<numbchild;c++){
            // Initialisation of temporary result
            float * DistSum_tmp=new float[numelmasked];
            for(int i=0;i<numelmasked;i++){
                DistSum_tmp[i]=0;
            }
            this->GetChild(c)->MakeWeightedDist(DistSum_tmp); // recursive part

            float * WeightedDist_PTR=WeightedDist;
            float * DistSum_tmp_PTR=DistSum_tmp;
            float NormWeightUsed=this->GetChild(c)->GetNormWeight();
            for(int i=0;i<numelmasked;i++,WeightedDist_PTR++,DistSum_tmp_PTR++){
                (*WeightedDist_PTR)+=NormWeightUsed*(*DistSum_tmp_PTR);
            }
            delete [] DistSum_tmp;
            DistSum_tmp=NULL;
        }
    }
    return;
}

void TreeEM::UpdateNonNormResp(SEG_PARAMETERS * segment_param){
    int numbchild = this->GetNumberChildren();
    int numelmasked=this->GetNumberMaskedElements();
    // Allocation and initialisation of the needed float array
    float * NonNormRespUpdate=new float[numelmasked];
    for(int i=0;i<numelmasked;i++){
        NonNormRespUpdate[i]=0;
    }
    this->MakeNonNormWeightedSum(NonNormRespUpdate, segment_param);
    // SetNormResp
    this->SetNormResp(NonNormRespUpdate);
    // Clearing Memory
    if(NonNormRespUpdate!=NULL){
        delete [] NonNormRespUpdate;
        NonNormRespUpdate=NULL;
    }
    // Recursivity
    for(int c=0;c<numbchild;c++){
        this->GetChild(c)->UpdateNonNormResp(segment_param);
    }
    return;
}

void TreeEM::UpdateNonNormRespCEM(){
    int numbchild = this->GetNumberChildren();
    int numelmasked=this->GetNumberMaskedElements();
    // Allocation and initialisation of the needed float array
    float * NonNormRespUpdate=new float[numelmasked];
    for(int i=0;i<numelmasked;i++){
        NonNormRespUpdate[i]=0;
    }
    this->MakeNonNormWeightedSumCEM(NonNormRespUpdate);
    // SetNormResp
    this->SetNormResp(NonNormRespUpdate);
    // Clearing Memory
    if(NonNormRespUpdate!=NULL){
        delete [] NonNormRespUpdate;
        NonNormRespUpdate=NULL;
    }
    // Recursivity
    for(int c=0;c<numbchild;c++){
        this->GetChild(c)->UpdateNonNormRespCEM();
    }
    return;
}

void TreeEM::UpdateNormResp(){
    int CountNormRespUpdatePb=0;
    int numbchild = this->GetNumberChildren();
    int numelmasked=this->GetNumberMaskedElements();
    if(this->GetParent()!=NULL){
        float * RootNonNormResp=this->FindRoot()->GetNormResp();
        float * NonNormRespToNormalise=this->GetNormResp();
        // Initialisation of the needed pointers for the normalisation
        float * RootNonNormResp_PTR=RootNonNormResp;
        float * NonNormRespToNormalise_PTR=NonNormRespToNormalise;
        // Normalisation if not root
        float MaxDifferenceToOne=0;
        int CountPbZero=0;
        int CountNeg=0;
        for(int i=0;i<numelmasked;i++,RootNonNormResp_PTR++,NonNormRespToNormalise_PTR++){
            if(*RootNonNormResp_PTR>-1E-6){
                if(*RootNonNormResp_PTR<*NonNormRespToNormalise_PTR){
                    //                cout<<"Possible problem at index "<<i<<endl;
                }
                if(*RootNonNormResp_PTR==0){
                    CountPbZero++;
//                    cout<<"RootNonNormResp is equal to 0 at index "<<i<<endl;
                    *NonNormRespToNormalise_PTR=0;
                    *NonNormRespToNormalise_PTR=1E-6;
                }
                else{
                *NonNormRespToNormalise_PTR=(*NonNormRespToNormalise_PTR)/(*RootNonNormResp_PTR);
                }
                if(*NonNormRespToNormalise_PTR>1+1E-6){
                    CountNormRespUpdatePb++;
                    if(*NonNormRespToNormalise_PTR-1>MaxDifferenceToOne){
                        MaxDifferenceToOne=*NonNormRespToNormalise_PTR-1;
                    }
                }
            }
            else{
                CountNeg++;
                *NonNormRespToNormalise_PTR=1.0;
//                cout<<"RootNonNormResp is negative at index "<<i<<endl;
            }
        }
        if(CountPbZero>0){
//        cout<< "CountZeroPb at root is "<<CountPbZero<<endl;
        }
        if(CountNeg>0){
//            cout<<"CountNeg is "<<CountNeg<<endl;
        }
        if(CountNormRespUpdatePb>0){
//            cout<<"CountNormRespPb is "<<CountNormRespUpdatePb<<" and the max difference is "<<MaxDifferenceToOne<<endl;
        }
    }
    else{
        //        cout<<"We do not change the norm resp of the root at first"<<endl;
    }

    // Recursivity
    for(int c=0;c<numbchild;c++){
        this->GetChild(c)->UpdateNormResp();
    }
    return;
}

void TreeEM::UpdateHardSeg(){
    if(!this->IsRoot()){
        return this->GetParent()->UpdateHardSeg();
    }
    else{
        this->MakeHardSeg();
    }
}

void TreeEM::UpdateNormRespRoot(){
    int numelmasked=this->GetNumberMaskedElements();
    float * RootNonNormResp=this->FindRoot()->GetNormResp();
    float * RootNonNormResp_PTR=RootNonNormResp;
    for(int i=0;i<numelmasked;i++,RootNonNormResp_PTR++){
        if(*RootNonNormResp_PTR>0){
            *RootNonNormResp_PTR/=*RootNonNormResp_PTR;
        }
        else{
            *RootNonNormResp_PTR=1;
        }
    }
}

//void TreeEM::UpdateNormResp(){
//    int CountNormRespUpdatePb=0;
//    int numbchild = this->GetNumberChildren();
//    int numelmasked=this->GetNumberMaskedElements();
//    // Allocation and initialisation of the needed float array
//    float * RootNonNormResp=new float [numelmasked];
//    float * NonNormRespToNormalise=new float[numelmasked];
//    float * NormRespUpdate=new float[numelmasked];
//    for(int i=0;i<numelmasked;i++){
//        RootNonNormResp[i]=0;
//        NonNormRespToNormalise[i]=0;
//        NormRespUpdate[i]=0;
//    }
//    this->FindRoot()->MakeNonNormWeightedSum(RootNonNormResp);
//    this->MakeNonNormWeightedSum(NonNormRespToNormalise);
//    // Initialisation of the needed pointers for the normalisation
//    float * RootNonNormResp_PTR=RootNonNormResp;
//    float * NonNormRespToNormalise_PTR=NonNormRespToNormalise;
//    float * NormRespUpdate_PTR=NormRespUpdate;
//    // Normalisation
//    for(int i=0;i<numelmasked;i++,NormRespUpdate_PTR++,RootNonNormResp_PTR++,NonNormRespToNormalise_PTR++){
//        if(*RootNonNormResp_PTR>0){
//            *NormRespUpdate_PTR=(*NonNormRespToNormalise_PTR)/(*RootNonNormResp_PTR);
//            if(*NormRespUpdate_PTR>1+1E-6){
//                CountNormRespUpdatePb++;
//            }
//        }
//        else{
//            *NormRespUpdate_PTR=1.0;
//        }
//    }
////    cout<<"CountNormRespUpdatePb is "<<CountNormRespUpdatePb<<endl;
//    // SetNormResp
//    this->SetNormResp(NormRespUpdate);
//    // Clearing Memory
//    if(RootNonNormResp!=NULL){
//        delete [] RootNonNormResp;
//        RootNonNormResp=NULL;
//    }
//    if(NonNormRespToNormalise!=NULL){
//        delete [] NonNormRespToNormalise;
//        NonNormRespToNormalise=NULL;
//    }
//    if(NormRespUpdate!=NULL){
//        delete [] NormRespUpdate;
//        NormRespUpdate=NULL;
//    }
//    // Recursivity
//    for(int c=0;c<numbchild;c++){
//        this->GetChild(c)->UpdateNormResp();
//    }
//return;
//}

float * TreeEM::MakeUniformDistribution(){
    if(this->GetDistributionType()!=2){
        cout<<"This is not a leaf with uniform distribution : cannot be made uniform"<<endl;
        return NULL;
    }
    int numelmasked=this->GetNumberMaskedElements();
    float * UniformDistOutput=new float[numelmasked];
    for(int i=0;i<numelmasked;i++){
        UniformDistOutput[i]=1;
    }
    return UniformDistOutput;
}

float* TreeEM::MakeGaussianDistribution(){
    //cout<<"Updating Gaussian Distribution"<<endl;
    if (this->GetDistributionType()!=1) {
        cout<<"This is not a leaf with Gaussian distribution : cannot be updated as wanted"<<endl;
        return NULL;
    }
    // Distribution put back to zero every where
    //    this->ReputDistributionToZero();
//    int numel=this->GetNumberElements();
    int numelmasked=this->GetNumberMaskedElements();
    int numbmodal=this->GetNumberModalities();
    float * GaussianDistOutput=new float[numelmasked];
    for(int i=0;i<numelmasked;i++){
        GaussianDistOutput[i]=0;
    }
    //    float * DataPointer=this->MakeDataBFCorrected();
    float * DataPointer=this->GetDataBFCorrected();
//    matrix<float> VarianceMatrix=matrix<float>(numbmodal);
    float * VarianceBis=this->GetVariance();
    float VarianceToUse[MaxNumbModal*MaxNumbModal];
    for(int m1=0;m1<MaxNumbModal;m1++){
        for(int m2=0;m2<MaxNumbModal;m2++){
            if(m1<numbmodal&&m2<numbmodal){
                VarianceToUse[m1+m2*MaxNumbModal]=VarianceBis[m1+m2*numbmodal];
            }
            else{
                VarianceToUse[m1+m2*MaxNumbModal]=0;
            }
        }
    }
    float * VarianceToUseTest=new float[numbmodal*numbmodal];
    for(int i=0;i<numbmodal*numbmodal;i++){
        VarianceToUseTest[i]=VarianceBis[i];
    }
    float * MeanBis=this->GetMean();
    float MeanToUse[MaxNumbModal];
    for(int i=0;i<MaxNumbModal;i++){
        if(i<numbmodal){
            MeanToUse[i]=MeanBis[i];
        }
        else{
            MeanToUse[i]=0;
        }

    }
    //    for(int m1=0;m1<MaxNumbModal;m1++){
    //        for(int m2=0;m2<MaxNumbModal;m2++){
    //            VarianceMatrix.setvalue(m1,m2,m1==m2?1:0);
    //        }
    //    }
//    for (int m1=0; m1<numbmodal; m1++) {
//        for (int m2=0; m2<numbmodal; m2++) {
//            VarianceMatrix.setvalue(m1, m2, VarianceToUse[m1+m2*MaxNumbModal]);

//        }
//    }

    // Calculation of the factor in front of the exponential
//    float DeterminantVariance=VarianceMatrix.determinant();
    float DeterminantVariance=determinant(VarianceToUseTest,numbmodal);
    //cout<<"the Variance determinant is "<<DeterminantVariance<<endl;
    float NormalisationFactor=1.0/(float)(powf(2*M_PI , (float)((float)numbmodal/2.0))*powf(DeterminantVariance, 0.5));
    //cout <<"The normalisation factor is "<<NormalisationFactor<<endl;
    //Initialisation of the needed element to calculate the inside of the exponential
    float * Distribution_PTR=GaussianDistOutput;
    float InvertedVariance[MaxNumbModal*MaxNumbModal];
    for(int i=0;i<MaxNumbModal*MaxNumbModal;i++){
        InvertedVariance[i]=0;
    }
    if (numbmodal>1) {
//        VarianceMatrix.invert();

//        float DetTest=determinant(VarianceToUseTest,numbmodal);
        invertMatrix(VarianceToUseTest,numbmodal);

//        bool success;
        for(int m1=0;m1<numbmodal;m1++){
            for (int m2=0; m2<numbmodal; m2++) {
                InvertedVariance[m1+m2*MaxNumbModal]=VarianceToUseTest[m1+m2*numbmodal];
//                VarianceMatrix.getvalue(m1, m2, InvertedVariance[m1+m2*MaxNumbModal], success);
            }
        }
        delete[] VarianceToUseTest;
        VarianceToUseTest=NULL;

    }
    else{
//        bool success;
//        VarianceMatrix.getvalue(0, 0, InvertedVariance[0], success);
        *InvertedVariance=1.0/(*VarianceToUseTest);
        delete[] VarianceToUseTest;
        VarianceToUseTest=NULL;
    }
    PrecisionTYPE temp;

    // Calculation of the inside of the exponential
    //cout<<InvertedVariance[0]<<endl;
    //        int * L2S_PTR=this->GetL2S();
    for (int i=0; i<numelmasked; i++,Distribution_PTR++) {

        //        if (L2S_PTR[i]>=0){
        for (int m1=0; m1<numbmodal; m1++) {
            temp=0;
            for (int m2=0; m2<numbmodal; m2++) {
                temp+=(PrecisionTYPE)(DataPointer[i+m2*numelmasked]-MeanToUse[m2])*InvertedVariance[m2+m1*MaxNumbModal];
            }
            *Distribution_PTR=(*Distribution_PTR)+temp*(DataPointer[i+m1*numelmasked]-MeanToUse[m1]);
            //            }
            //            Distribution_PTR++;
        }
    }

    // Filling the distribution array with the value
    Distribution_PTR=GaussianDistOutput;


    for (int i=0;i<numelmasked;i++,Distribution_PTR++){
        *Distribution_PTR=NormalisationFactor*expf(-0.5*(*Distribution_PTR));
    }
    // Clearing space needing for inverse of variance

    //    if (DataPointer!=NULL) {
    //        delete [] DataPointer;
    //        DataPointer=NULL;
    //    }
    //    string Filename="/Users/Carole/Documents/PhD/TryDistributionDataT1Nifty";
    //    nifti_image * PartialRes=SavePartialResult(this->GetPartialResult(DISTRIBUTION),this->GetDataImage(),(char *)Filename.c_str());
    //    nifti_image_write(PartialRes);
    //    nifti_image_free(PartialRes);    return;
    return GaussianDistOutput;
}

float * TreeEM::NormalisedPriorsForCCL(){
    vector<TreeEM*> LeavesVector=this->GetAllLeaves();

    int numbLeaves=LeavesVector.size();
    int numelmasked=this->GetNumberMaskedElements();
//    int numel=this->GetNumberElements();

    // Initialisation of the result
    float * NormalisedPriorsResult=new float[numbLeaves*numelmasked];
    for(int i = 0;i<numbLeaves*numelmasked;i++){
        NormalisedPriorsResult[i]=1;
    }
    // Copy of the adapted priors into the results
    for(int l=0;l<numbLeaves;l++){
        int * S2L_PTR=this->GetS2L();
        float * PriorsLeaves=LeavesVector[l]->GetPriorsAdapted();
        float * MRFValue = LeavesVector[l]->GetMRF();
        float PartNormWeight=LeavesVector[l]->GetPartNormWeight();
        for(int i=0;i<numelmasked;i++){
            if(PriorsLeaves!=NULL){
                NormalisedPriorsResult[i]*=PriorsLeaves[S2L_PTR[i]];
            }
            if(MRFValue!=NULL){
                NormalisedPriorsResult[i]*=MRFValue[i];
            }
            if(LeavesVector[l]->GetPriorsDirect()==NULL){
                NormalisedPriorsResult[i]*=PartNormWeight;
            }
        }
    }

    // Normalisation step
    int CountNormZero=0;
    for(int i=0;i<numelmasked;i++){
        float NormFactor=0;
        // Determination of normalisation factor so that sum at each voxel is to 1 after normalisation
        for(int l=0;l<numbLeaves;l++){
            NormFactor+=NormalisedPriorsResult[i+l*numelmasked];
        }
        //Normalisation at each voxel
        if(NormFactor==0){
            CountNormZero++;
        }
        else{
            for(int l=0;l<numbLeaves;l++){
                NormalisedPriorsResult[i+l*numelmasked]/=NormFactor;
            }
        }
    }
    if(CountNormZero>0){
        cout<<"In normalisation of adapted atlas + MRF number of zero pb is "<<CountNormZero<<endl;
    }
    return NormalisedPriorsResult;

}

float TreeEM::CompleteLogLikelihoodFinBis(){
    PrecisionTYPE CompleteLogLikelihood=0;
    int numelmasked=this->GetNumberMaskedElements();
//    int numel=this->GetNumberElements();
//    int numbchild=this->GetNumberChildren();
    vector<TreeEM*> Leaves=this->GetAllLeaves();
    float * NormalisedPriors=this->NormalisedPriorsForCCL();
    int numbLeaves=Leaves.size();

    for (int c=0; c<numbLeaves; c++) {

        float * NormResp_PTR=Leaves[c]->GetNormResp();
        float * Priors_PTR=&NormalisedPriors[c*numelmasked];
        float * Distribution= new float[numelmasked];
        float * Distribution_PTR=Distribution;
        for(int i=0;i<numelmasked;i++){
            Distribution[i]=0;
        }
        Leaves[c]->MakeWeightedDist(Distribution);

        for(int i=0;i<numelmasked;i++,NormResp_PTR++,Distribution_PTR++,Priors_PTR++){
            if(*Distribution_PTR>0.00001){
                CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*log(*Distribution_PTR);
            }
            if(*Priors_PTR>0.00001){
                CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*log(*Priors_PTR);
            }
        }
        if(Distribution!=NULL){
            delete [] Distribution;
            Distribution = NULL;
        }
    }
    if(NormalisedPriors!=NULL){
        delete [] NormalisedPriors;
        NormalisedPriors=NULL;
    }
    return (float)CompleteLogLikelihood;
}

// REMARK : Be really careful when Mask is changed from one layer to another. Always check that the right amount of memory is then allocated
float TreeEM::CompleteLogLikelihoodFin(SEG_PARAMETERS * segment_param){
    PrecisionTYPE CompleteLogLikelihood=0;
    int numelmasked=this->GetNumberMaskedElements();
    int numel=this->GetNumberElements();
    int numbchild=this->GetNumberChildren();
    vector<TreeEM*> Leaves=this->GetAllLeaves();
    int numbLeaves=Leaves.size();
    for (int c=0; c<numbLeaves; c++) {
        int * L2S_PTR=this->GetL2S();
        float * NormResp_PTR=Leaves[c]->GetNormResp();
        float * Priors_PTR=NULL;
        float * Priors_PTR2=NULL;
        if(Leaves[c]->GetPriors()!=NULL){
//            Priors_PTR=new float[numel];
//            Priors_PTR2=Priors_PTR;
//            Priors_PTR2=Leaves[c]->AdaptPriors(segment_param);
//            Priors_PTR=static_cast<float *>(Leaves[c]->GetPriors()->data);
            Priors_PTR=Leaves[c]->GetPriorsAdapted();
        }
        float * Distribution= new float[numelmasked];
        float * Distribution_PTR=Distribution;
        for(int i=0;i<numelmasked;i++){
            Distribution[i]=0;
        }
        Leaves[c]->MakeWeightedDist(Distribution);

        if (Priors_PTR!=NULL) {
            for(int i=0;i<numel;i++){
                if ((*L2S_PTR)>=0){
                    if((*Distribution_PTR)>0.00001){
                        CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(*Distribution_PTR);
                    }
                    else{
                        CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(0.00001);
                    }

                    if((*Priors_PTR)>0.00001){
                        CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(*Priors_PTR);
                    }
                    else{
                        CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(0.00001);
                    }
                    if((Leaves[c]->GetPriorsDirect()==NULL)){
                        CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(Leaves[c]->GetPartNormWeight());
                    }


                    NormResp_PTR++;
                    Distribution_PTR++;
                }
                Priors_PTR++;
                L2S_PTR++;
            }
//            if(Priors_PTR2!=NULL){
//                delete [] Priors_PTR2;
//                Priors_PTR2=NULL;
//            }
        }

        else {
            for(int i=0;i<numelmasked;i++,Distribution_PTR++,NormResp_PTR++){

                if((*Distribution_PTR)>0){
                    CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(*Distribution_PTR);
                }
                else{
                    CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(0.00001);
                }

                if(this->GetChild(Leaves[c]->FindGeneralClass())->GetNormWeight()>0){
                    CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(this->GetChild(Leaves[c]->FindGeneralClass())->GetNormWeight());
                }
                else{
                    CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(0.00001);
                }
                if(!Leaves[c]->GetParent()->IsRoot()){
                    CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(Leaves[c]->GetNormWeight());
                }

            }

        }
        if(Distribution!=NULL){
            delete [] Distribution;
            Distribution=NULL;
        }
    }

    return (float)CompleteLogLikelihood;
}

// REMARK : Be really careful when Mask is changed from one layer to another. Always check that the right amount of memory is then allocated
float TreeEM::CompleteLogLikelihood(){
    PrecisionTYPE CompleteLogLikelihood=0;
    int numelmasked=this->GetNumberMaskedElements();
    int numel=this->GetNumberElements();
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        int * L2S_PTR=this->GetL2S();
        float * NormResp_PTR=this->GetChild(c)->GetNormResp();
        float * Priors_PTR=NULL;
        if(this->GetChild(c)->GetPriors()!=NULL){
            Priors_PTR=static_cast<float *>(this->GetChild(c)->GetPriors()->data);
        }
        float * Distribution= new float[numelmasked];
        float * Distribution_PTR=Distribution;
        for(int i=0;i<numelmasked;i++){
            Distribution[i]=0;
        }
        this->GetChild(c)->MakeWeightedDist(Distribution);

        if (Priors_PTR!=NULL) {
            for(int i=0;i<numel;i++){
                if ((*L2S_PTR)>=0){
                    if((*Distribution_PTR)>0.00001){
                        CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(*Distribution_PTR);
                    }
                    else{
                        CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(0.00001);
                    }

                    if((*Priors_PTR)>0.00001){
                        CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(*Priors_PTR);
                    }
                    else{
                        CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(0.00001);
                    }

                    NormResp_PTR++;
                    Distribution_PTR++;
                }
                Priors_PTR++;
                L2S_PTR++;
            }
        }

        else {
            for(int i=0;i<numelmasked;i++,Distribution_PTR++,NormResp_PTR++){

                if((*Distribution_PTR)>0){
                    CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(*Distribution_PTR);
                }
                else{
                    CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(0.00001);
                }

                if(this->GetChild(c)->GetNormWeight()>0){
                    CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(this->GetChild(c)->GetNormWeight());
                }
                else{
                    CompleteLogLikelihood+=(PrecisionTYPE)(*NormResp_PTR)*logf(0.00001);
                }

            }

        }
        if(Distribution!=NULL){
            delete [] Distribution;
            Distribution=NULL;
        }
    }

    return (float)CompleteLogLikelihood;
}

void TreeEM::EMCompleteIteration(float & LogLikelihood, float & OldLogLikelihood, int & Iteration,SEG_PARAMETERS * segment_param){

    OldLogLikelihood=LogLikelihood;


    this->UpdateNonNormResp(segment_param);
    //cout<<"Update NonNormResp done"<<endl;

    this->UpdateNormResp();
//    this->IsNormRespValidGeneral();
//    this->IsNormRespValidLeaves();
    //cout<<"Update NormResp done"<<endl;
    this->UpdateNormRespRoot();

    LogLikelihood=this->CompleteLogLikelihoodFinBis();

    if(segment_param->flag_MRF){
        if(!segment_param->flag_GMatrix){
            float * GMatrixToUpdate=this->MRFOptSolveLS();
            this->SetGMatrix(GMatrixToUpdate);
            delete [] GMatrixToUpdate;
            GMatrixToUpdate=NULL;
        }
        this->UpdateMRF();
    }


//    LogLikelihood=this->GetLogLikelihood();

    this->UpdateNonNormWeights();
    //cout<<"Update NonNormWeights done"<<endl;
    this->UpdateNormWeights();
    this->UpdatePriorsAdapted(segment_param);

//    CompleteLogLikelihood=this->CompleteLogLikelihoodFin();
    if(LogLikelihood>this->CompleteLogLikelihoodFinBis()){
//        cout<<"Pb in maximization"<<endl;
    }
    LogLikelihood=this->CompleteLogLikelihoodFinBis();

    this->UpdateParameters();

    if(LogLikelihood>this->CompleteLogLikelihoodFinBis()){
//        cout<<"Pb in maximization parameters"<<endl;
    }

//    CompleteLogLikelihood=this->CompleteLogLikelihoodFin();
//    LL=this->GetLogLikelihood();
    //cout<<"Update Parameters done"<<endl;

    //    this->UpdateDistribution();
    //cout<<"Update distribution done"<<endl;
//    this->UpdateNonNormResp();
//    //cout<<"Update NonNormResp done"<<endl;

//    this->UpdateNormResp();
//    //cout<<"Update NormResp done"<<endl;
//    this->UpdateNormRespRoot();
//    this->UpdateNonNormWeights();
//    //cout<<"Update NonNormWeights done"<<endl;
//    this->UpdateNormWeights();
//    bool TestWeights=this->AreWeightsValid();
    if (BFFlag) {
        this->UpdateBFCoeffs();
        //        this->UpdateBFCorrection();
        this->UpdateDataBFCorrected();
    }
//    CompleteLogLikelihood=this->CompleteLogLikelihood();
//    CompleteLogLikelihood=this->CompleteLogLikelihoodFin();

    LogLikelihood=this->GetLogLikelihood();
    cout<< LogLikelihood<<" ";
    Iteration++;
    if(this->GetParent()==NULL){
//            this->SaveAllClasses("/Users/Carole/Documents/PhD/TestResult.nii.gz",segment_param);
    }
}

// Performs a complete EM optimisation
void TreeEM::RunFullEM(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration,SEG_PARAMETERS * segment_param){
    bool ValidityInitialisationTree = 1;// this->IsTreeWellInitialised();
    if (!ValidityInitialisationTree){
        cout<< "Tree not well initialised, EM cannot be performed" << endl;
        return;
    }
    cout<<"Tree well initialised we can run the EM"<<endl;
    while (((CompleteLogLikelihood-OldCompleteLogLikelihood)/fabs(OldCompleteLogLikelihood)>Threshold || Iteration<6) && Iteration<MaxIteration) {
        //        cout<<"we are at iteration "<<Iteration<<endl;
        //        cout<< CompleteLogLikelihood <<" and the old CLL "<< OldCompleteLogLikelihood<<endl;
        this->EMCompleteIteration(CompleteLogLikelihood, OldCompleteLogLikelihood, Iteration,segment_param);
//        this->SaveGeneralClasses("/Users/Carole/Documents/PhD/MRFResultWOG.nii.gz",segment_param);
    }
    cout << Iteration<<" iterations before convergence and the LogLikelihood at convergence is "<<CompleteLogLikelihood<<endl;
}

void TreeEM::CEMPartialIteration(float & LogLikelihood, float & OldLogLikelihood, int & Iteration,int c){
    OldLogLikelihood=LogLikelihood;


    this->UpdateNonNormRespCEM();
    //cout<<"Update NonNormResp done"<<endl;

    this->UpdateNormResp();
//    this->IsNormRespValidGeneral();
//    this->IsNormRespValidLeaves();
    //cout<<"Update NormResp done"<<endl;
    this->UpdateNormRespRoot();

//    LogLikelihood=this->GetLogLikelihood();

    this->UpdateNonNormWeightsCEM(c);
    //cout<<"Update NonNormWeights done"<<endl;
    this->UpdateNormWeights();

//    CompleteLogLikelihood=this->CompleteLogLikelihoodFin();
    if(LogLikelihood>this->GetLogLikelihoodCEM(c)){
//        cout<<"Pb in maximization"<<endl;
    }
    LogLikelihood=this->GetLogLikelihoodCEM(c);

    this->UpdateParametersCEM( c);

    if(LogLikelihood>this->GetLogLikelihoodCEM(c)){
//        cout<<"Pb in maximization parameters"<<endl;
    }


    if (BFFlag) {
        this->UpdateBFCoeffs();
        //        this->UpdateBFCorrection();
        this->UpdateDataBFCorrected();
    }

    LogLikelihood=this->GetLogLikelihoodCEM(c);
//    cout<< LogLikelihood<<" ";
    Iteration++;
    if(this->GetParent()==NULL){
        //    this->SaveAllClasses("/Users/Carole/Documents/PhD/TestResult");
    }
}

// Performs a complete CEM optimisation for child c
void TreeEM::RunPartialCEM(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration, int c){
    bool ValidityInitialisationTree = 1;// this->IsTreeWellInitialised();
    if (!ValidityInitialisationTree){
        cout<< "Tree not well initialised, EM cannot be performed" << endl;
        return;
    }
    cout<<"Tree well initialised we can run the EM"<<endl;
    if(this->GetNumberChildren()==0){ // If nothing to really optimise since not split
        cout<<"No embedded EM to run"<<endl;
        this->UpdateParametersCEM(c);
    }
    else{
    while (((CompleteLogLikelihood-OldCompleteLogLikelihood)/fabs(OldCompleteLogLikelihood)>Threshold || Iteration<6) && Iteration<MaxIteration) {
        //        cout<<"we are at iteration "<<Iteration<<endl;
        //        cout<< CompleteLogLikelihood <<" and the old CLL "<< OldCompleteLogLikelihood<<endl;
        this->CEMPartialIteration(CompleteLogLikelihood, OldCompleteLogLikelihood, Iteration,c);
    }
    }
    cout << Iteration<<" iterations before convergence and final CompleteLikelihood is "<<CompleteLogLikelihood<<endl;

}


void TreeEM::CEMCompleteIteration(float & LogLikelihood, float & OldLogLikelihood, int & Iteration,SEG_PARAMETERS * segment_param){
    OldLogLikelihood=LogLikelihood;
    int numbchild=this->GetNumberChildren();

    this->UpdateNonNormResp(segment_param);
    //cout<<"Update NonNormResp done"<<endl;

    this->UpdateNormResp();
//    this->IsNormRespValidGeneral();
//    this->IsNormRespValidLeaves();
    //cout<<"Update NormResp done"<<endl;
    this->UpdateNormRespRoot();

//    LogLikelihood=this->GetLogLikelihood();

    this->UpdateNonNormWeights();
    //cout<<"Update NonNormWeights done"<<endl;
    this->UpdateNormWeights();

//    CompleteLogLikelihood=this->CompleteLogLikelihoodFin();
    if(LogLikelihood>this->GetLogLikelihood()){
//        cout<<"Pb in maximization"<<endl;
    }
    LogLikelihood=this->GetLogLikelihood();

    if(numbchild>0){ // Applied Embedded EM if it is not a leaf
    for(int c=0;c<numbchild;c++){
        float LLc=0;
        float OldLLc=0;
        int IterationC=0;
        TreeEM * TmpCEMTree=this->GetChild(c)->CopyTree(NULL);
        TmpCEMTree->SetHardSeg(this->GetHardSeg());
        TmpCEMTree->RunPartialCEM(LLc,OldLLc,IterationC,c);
        this->GetChild(c)->ReplugParameters(TmpCEMTree);
        delete TmpCEMTree;
        TmpCEMTree=NULL;
    }
    }


    if(LogLikelihood>this->GetLogLikelihood()){
//        cout<<"Pb in maximization parameters"<<endl;
    }


    if (BFFlag) {
        this->UpdateBFCoeffs();
        //        this->UpdateBFCorrection();
        this->UpdateDataBFCorrected();
    }

    LogLikelihood=this->GetLogLikelihood();
//    cout<< LogLikelihood<<" ";
    Iteration++;
    if(this->GetParent()==NULL){
        //    this->SaveAllClasses("/Users/Carole/Documents/PhD/TestResult");
    }
}

void TreeEM::ReplugParameters(TreeEM* TreeToReplug){
    // Check if trees are of similar structures
    if(!this->IsStructureSimilar(TreeToReplug)){
        cout<<"The considered trees do not share a similar structure"<<endl;
        return;
    }
    else{
        int numbchild=this->GetNumberChildren();
        this->SetParameters(TreeToReplug->CopyParameters());
        for(int c=0;c<numbchild;c++){
            this->GetChild(c)->ReplugParameters(TreeToReplug->GetChild(c));
        }
    }
}

void TreeEM::RunFullCEM(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration,SEG_PARAMETERS * segment_param){
    while (((CompleteLogLikelihood-OldCompleteLogLikelihood)/fabs(OldCompleteLogLikelihood)>Threshold || Iteration<6) && Iteration<MaxIteration) {
        //        cout<<"we are at iteration "<<Iteration<<endl;
        //        cout<< CompleteLogLikelihood <<" and the old CLL "<< OldCompleteLogLikelihood<<endl;
        this->CEMCompleteIteration(CompleteLogLikelihood, OldCompleteLogLikelihood, Iteration,segment_param);
    }
    cout << Iteration<<" iterations before convergence and the final loglikelihood is "<<CompleteLogLikelihood<<endl;
}


TreeEM * TreeEM::RunFullBiASM(SEG_PARAMETERS * segment_param){
    float CompleteLogLikelihood=0;
    float OldCompleteLogLikelihood=0;
    int Iteration=0;
    cout<< "The weights of the general classes are : ";
    float * DPChildren=this->GetDPChildrenDirect();

    for (int c=0;c<this->GetNumberChildren();c++){
        cout<<this->GetChild(c)->GetNormWeight()<< "    ******      ";
    }
    BFFlag=segment_param->flag_Bias;
    BForder=segment_param->bias_order;
    Threshold=segment_param->Bias_threshold;
    MaxIteration=segment_param->maxIterationBF;
    MinIteration=segment_param->minIteration;
    NumbMaxLeaves=segment_param->maxNumbLeaves;
    KernelSize=segment_param->KernelSize;
    int numbmodal=this->GetNumberModalities();
//    segment_param->choiceInitSplit=segment_param->choiceInitSplit<1?1:segment_param->choiceInitSplit;
//    segment_param->choiceInitSplit=segment_param->choiceInitSplit>numbmodal?numbmodal:segment_param->choiceInitSplit;
    //    BFFlag=0;
//    this->SaveBFCorrectedData("/Users/Carole/Documents/PhD/BRATS_DataCorrected");
//    this->MakeIndFactorTot();
//    cout<<this->GetNumberChildren();
    cout<<this->MakeIndFactorTotMod()<<" "<<this->MakeIndFactorTot()<<endl;
    this->IsNormRespValidGeneral();
    this->IsNormRespValidLeaves();
    float AtlasWeight=segment_param->AtlasWeight;
//    segment_param->AtlasWeight=0;
    if(!segment_param->flag_intxt){
            segment_param->AtlasWeight=0;
    this->RunFullEM(CompleteLogLikelihood, OldCompleteLogLikelihood, Iteration,segment_param);
    if(segment_param->flag_EMfirst_out){
        this->SaveAllClasses(segment_param->filename_EMfirst,segment_param);
        string FilenameEMFirst=segment_param->filename_EMfirst;
        FilenameEMFirst+=".txt";
        char FilenameToUse[200];
        strcpy(FilenameToUse,FilenameEMFirst.c_str());
        this->SaveTreeInTextFile(FilenameToUse,segment_param);
    }
    this->SavePriorsAdapted(segment_param);
    if(BFFlag){
//        this->SaveBFCorrectedData("/Users/Carole/Documents/PhD/BRATS_DataCorrected");
//        this->SaveBFCorrection("/Users/Carole/Documents/PhD/BRATS_BFCorrection");
        this->SetIndFactor(this->MakeIndFactorTot());
    }
    }
    this->SetIndFactor(this->MakeIndFactorTot());
    int CountRunEM=1;
    cout<< "The weights of the general classes are : ";
    for (int c=0;c<this->GetNumberChildren();c++){
        cout<<this->GetChild(c)->GetNormWeight()<< "    ******      ";
    }
    MaxIteration=segment_param->maxIteration;
    Threshold=segment_param->ConvThreshold;
    BFFlag=0;
    segment_param->AtlasWeight=AtlasWeight;
//    if(this->GetFlagOutliers()>0){
//        segment_param->AtlasWeight=1;
//        this->UpdatePriorsAdapted(segment_param);
//        int LastChild=this->GetNumberChildren()-1;
//        float * OutlierClass= this->GetChild(LastChild)->GetPriorsAdaptedDirect();
//        nifti_image * AdaptedOutlierAtlas= nifti_copy_nim_info(this->GetDataImage());
//        AdaptedOutlierAtlas->dim[0]=3;
//        AdaptedOutlierAtlas->dim[4]=1;
//        nifti_update_dims_from_array(AdaptedOutlierAtlas);
//        AdaptedOutlierAtlas->data=(void *)calloc(AdaptedOutlierAtlas->nvox,sizeof(float));
//        float * AdaptedOutlierAtlas_PTR=static_cast<float*>(AdaptedOutlierAtlas->data);
//        int numel=this->GetNumberElements();
//        for(int i=0;i<numel;i++,AdaptedOutlierAtlas_PTR++,OutlierClass++){
//            *AdaptedOutlierAtlas_PTR=*OutlierClass;
//        }
//            nifti_set_filenames(AdaptedOutlierAtlas, "/Users/Carole/Documents/PhD/OutlierClass.nii.gz", 0, 0);
//            nifti_image_write(AdaptedOutlierAtlas);
//        this->GetChild(LastChild)->SetPriors(AdaptedOutlierAtlas);
//        if(this->GetChild(LastChild)->GetPriorsDirect()!=NULL){
//            //delete [] this->GetChild(LastChild)->GetPriorsDirect();
//            this->GetChild(LastChild)->SetPriors(AdaptedOutlierAtlas);
//        }
//    }
    TreeEM * TreeSplit=this;
    TreeEM * TreeMerge=this;
    if(segment_param->flag_BiASM){
        cout<<"The BiASM part will be conducted"<<endl;
        MergeKLD * MergeTest=TreeMerge->GetToMerge();
        SplitKLD * SplitTest=TreeSplit->GetToSplit(segment_param);
        int numberAllLeaves=this->GetNumberAllLeaves();
        bool TestNormResp=TreeSplit->AreNormRespValid();
        bool TestNormWeights=TreeSplit->AreWeightsValid();
        if(!TestNormResp){
            cout << "Norm resp non valid for TreeMerge"<<endl;
        }
        if(!TestNormWeights){
            cout << "Norm weights non valid for TreeMerge"<<endl;
        }
        bool AcceptanceDecision=0;
        while ((SplitTest!=NULL&& numberAllLeaves<NumbMaxLeaves) || ((MergeTest!=NULL)&&numberAllLeaves<=NumbMaxLeaves)) {
            AcceptanceDecision=0;
            if(SplitTest!=NULL && numberAllLeaves<NumbMaxLeaves){
                CountRunEM+=1;
            }
            if(MergeTest!=NULL){
                CountRunEM+=1;
            }

            int numbchild=TreeSplit->GetNumberChildren();
            if(numberAllLeaves<=NumbMaxLeaves){
                TreeSplit=TreeSplit->RunSplitOperation(SplitTest,AcceptanceDecision,segment_param);

//                for (int c=0; c<numbchild; c++) {
//                    (TreeSplit->GetChild(c))->PutAllLeavesToChildrenLevel();
//                }
                TreeSplit->PutAllLeavesToMainNodesChildrenLevel();
            }

            TestNormResp=TreeSplit->AreNormRespValid();
            TestNormWeights=TreeSplit->AreWeightsValid();
            if(!TestNormResp){
                cout << "Norm resp non valid for TreeSplit"<<endl;
            }
            if(!TestNormWeights){
                cout << "Norm weights non valid for TreeSplit"<<endl;
            }
            if(AcceptanceDecision){
                TreeSplit->ClearSMChecks();
//                MergeTest=NULL;
//                MergeTest=TreeSplit->GetToMerge();
                AcceptanceDecision=0;
//                if(MergeTest!=NULL){
//                    delete MergeTest;
//                    MergeTest=NULL;
//                }
//                MergeTest=TreeSplit->GetToMerge();
            }
            MergeTest=NULL;
            MergeTest=TreeSplit->GetToMerge();
            TreeMerge=TreeSplit->RunMergeOperation(MergeTest,AcceptanceDecision,segment_param);
            TestNormResp=TreeMerge->AreNormRespValid();
            TestNormWeights=TreeMerge->AreWeightsValid();

            if(!TestNormResp){
                cout << "Norm resp non valid for TreeMerge"<<endl;
            }
            if(!TestNormWeights){
                cout << "Norm weights non valid for TreeMerge"<<endl;
            }
            TreeSplit=TreeMerge;
//            TreeSplit->SaveTreeInTextFile("/Users/Carole/Documents/PhD/DataTree.txt");

            SplitTest=TreeSplit->GetToSplit(segment_param);
            MergeTest=TreeMerge->GetToMerge();
            numberAllLeaves=TreeSplit->GetNumberAllLeaves();
            numbchild=TreeMerge->GetNumberChildren();
            cout<< "The weights of the general classes are : ";
            for (int c=0;c<numbchild;c++){
                cout<<TreeMerge->GetChild(c)->GetNormWeight()<< "    ******      ";
            }
            cout<<endl;
//            TreeSplit->SaveMRFImage(segment_param);
//            TreeSplit->SaveAllClasses("/Users/Carole/Documents/PhD/TmpResult.nii.gz",segment_param);
//            TreeSplit->SaveTreeInTextFile("/Users/Carole/Documents/PhD/DataTree.txt",segment_param);
        }
        if(MergeTest!=NULL){
            delete MergeTest;
            MergeTest=NULL;
        }
        if(SplitTest!=NULL){
            delete SplitTest;
            SplitTest=NULL;
        }

    }
    TreeEM * TreeResult=TreeSplit->CopyTree(NULL);
    //    if(this!=NULL){
    //        delete this;
    //    }
    if(TreeSplit!=NULL){
        delete TreeSplit;
        TreeSplit=NULL;
    }
    cout<<CountRunEM<<" have been performed "<<endl;
    return TreeResult;
}

TreeEM * TreeEM::RunFullBiASM_bis(SEG_PARAMETERS * segment_param){
    float CompleteLogLikelihood=0;
    float OldCompleteLogLikelihood=0;
    int Iteration=0;
    cout<< "The weights of the general classes are : ";
    float * DPChildren=this->GetDPChildrenDirect();

    for (int c=0;c<this->GetNumberChildren();c++){
        cout<<this->GetChild(c)->GetNormWeight()<< "    ******      ";
    }
    BFFlag=segment_param->flag_Bias;
    BForder=segment_param->bias_order;
    Threshold=segment_param->Bias_threshold;
    MaxIteration=segment_param->maxIterationBF;
    MinIteration=segment_param->minIteration;
    NumbMaxLeaves=segment_param->maxNumbLeaves;
    KernelSize=segment_param->KernelSize;
    int numbmodal=this->GetNumberModalities();
    float AtlasWeight=segment_param->AtlasWeight;
    this->SaveBFCorrectedData(segment_param->filename_datacorrected);
    if(!segment_param->flag_intxt){
//        segment_param->AtlasWeight=0;
    this->RunFullEM(CompleteLogLikelihood, OldCompleteLogLikelihood, Iteration,segment_param);
    if(segment_param->flag_EMfirst_out){
        this->SaveAllClasses(segment_param->filename_EMfirst,segment_param);
        string FilenameEMFirst=nifti_makebasename(segment_param->filename_EMfirst.c_str());
        FilenameEMFirst+=".txt";
        char FilenameToUse[200];
        strcpy(FilenameToUse,FilenameEMFirst.c_str());
        this->SaveTreeInTextFile(FilenameToUse,segment_param);
    }
    this->SavePriorsAdapted(segment_param);
    if(BFFlag){
        this->SaveBFCorrectedData(segment_param->filename_datacorrected);

//        this->SaveBFCorrectedData("/Users/Carole/Documents/PhD/BRATS_DataCorrected");
//        this->SaveBFCorrection("/Users/Carole/Documents/PhD/BRATS_BFCorrection");
        this->SetIndFactor(this->MakeIndFactorTot());
    }
    }
    this->SetIndFactor(this->MakeIndFactorTot());
    int CountRunEM=1;
    cout<< "The weights of the general classes are : ";
    for (int c=0;c<this->GetNumberChildren();c++){
        cout<<this->GetChild(c)->GetNormWeight()<< "    ******      ";
    }
    MaxIteration=segment_param->maxIteration;
    Threshold=segment_param->ConvThreshold;
    BFFlag=0;
//    segment_param->AtlasWeight=AtlasWeight;
//    if(this->GetFlagOutliers()>0){
//        this->UpdatePriorsAdapted(segment_param);
//        int LastChild=this->GetNumberChildren()-1;
//        float * OutlierClass= this->GetChild(LastChild)->GetPriorsAdaptedDirect();
//        nifti_image * AdaptedOutlierAtlas= nifti_copy_nim_info(this->GetDataImage());
//        AdaptedOutlierAtlas->dim[0]=3;
//        AdaptedOutlierAtlas->dim[4]=1;
//        nifti_update_dims_from_array(AdaptedOutlierAtlas);
//        AdaptedOutlierAtlas->data=(void *)calloc(AdaptedOutlierAtlas->nvox,sizeof(float));
//        float * AdaptedOutlierAtlas_PTR=static_cast<float*>(AdaptedOutlierAtlas->data);
//        int numel=this->GetNumberElements();
//        for(int i=0;i<numel;i++,AdaptedOutlierAtlas_PTR++,OutlierClass++){
//            *AdaptedOutlierAtlas_PTR=*OutlierClass;
//        }
////            nifti_set_filenames(AdaptedOutlierAtlas, "/Users/Carole/Documents/PhD/OutlierClass.nii.gz", 0, 0);
////            nifti_image_write(AdaptedOutlierAtlas);
//        this->GetChild(LastChild)->SetPriors(AdaptedOutlierAtlas);
//        if(this->GetChild(LastChild)->GetPriorsDirect()!=NULL){
//            //delete [] this->GetChild(LastChild)->GetPriorsDirect();
//            this->GetChild(LastChild)->SetPriors(AdaptedOutlierAtlas);
//        }
//        this->NormalisePriors();
//        this->NormalisePriorsAdapted();
//    }
    int numbchild=this->GetNumberChildren();
    if (segment_param->PriorsKept==0){
        for (int c=0;c<numbchild;c++){
            if(this->GetChild(c)->GetPriorsDirect()!=NULL){
                nifti_image_free(this->GetChild(c)->GetPriorsDirect());
                this->GetChild(c)->SetPriors(NULL);
            }
        }
    }
    if(segment_param->PriorsKept==2){
        for(int c=0;c<numbchild;c++){

            nifti_image * NewPriors=this->GetChild(c)->TransformNormRespIntoPriors(segment_param);
            this->GetChild(c)->SetPriors(NewPriors);
            float * NewPriorsAdapted=static_cast<float *>(NewPriors->data);
            this->GetChild(c)->SetPriorsAdapted(NewPriorsAdapted);
        }
    }
    this->MakeHardSeg();
    TreeEM * TreeSplit=this;
    TreeEM * TreeMerge=this;
    if(segment_param->flag_BiASM){
        cout<<"The BiASM part will be conducted"<<endl;

        int numberAllLeaves=this->GetNumberAllLeaves();
        bool TestNormResp=TreeSplit->AreNormRespValid();
        bool TestNormWeights=TreeSplit->AreWeightsValid();
        if(!TestNormResp){
            cout << "Norm resp non valid for TreeMerge"<<endl;
        }
        if(!TestNormWeights){
            cout << "Norm weights non valid for TreeMerge"<<endl;
        }
        bool AcceptanceDecisionSplit=0;
        bool AcceptanceDecisionMerge=0;
        int SeenSplit=0;
        int SeenMerge=0;
//        MergeKLD * MergeTest=TreeMerge->GetToMerge();
        int numbchild=this->GetNumberChildren();

//        int * CombTest=Combination(TreeSplit->GetNumberChildren(),2);
        int NumbCommonChanges=1;
        if(segment_param->flag_CommonChanges && segment_param->SMOrder){
            NumbCommonChanges=numbchild;
        }
        for(int c=NumbCommonChanges;c>0;c--){
            SeenSplit=0;
            SeenMerge=0;
            vector<SplitKLD *> SplitTest;

            vector<MergeKLD*> MergeTest;
            if(segment_param->SMOrder){
                SplitTest  =TreeSplit->GetSplitMoreVertical(c,segment_param);
              MergeTest=TreeMerge->GetMergeMoreVertical(c);
            }
            else{
                SplitTest=TreeSplit->OrderingSplittingLeaves(segment_param);
                MergeTest=TreeMerge->OrderingMergingLeaves();
            }

            cout<<"Now looking at "<<c<<" classes among "<<numbchild<<endl;
        while ((SplitTest.size()>0 && numberAllLeaves<NumbMaxLeaves) || ((MergeTest.size()>0)&&numberAllLeaves<=NumbMaxLeaves)) {
            TreeSplit->IsNormRespValidGeneral();
            TreeSplit->IsNormRespValidLeaves();
            AcceptanceDecisionSplit=0;
            AcceptanceDecisionMerge=0;
            if(SplitTest.size()>0 && numberAllLeaves<NumbMaxLeaves){
                CountRunEM+=1;
            }
            if(MergeTest.size()!=0){
                CountRunEM+=1;
            }

            int numbchild=TreeSplit->GetNumberChildren();
            if(numberAllLeaves<=NumbMaxLeaves && SplitTest.size()>0){
                cout<<SplitTest.size()<< "split operations to test";
                vector<SplitKLD*> SplitTested (SplitTest.begin(),SplitTest.begin()+c);
                TreeSplit=TreeSplit->RunSplitOperation(SplitTested,AcceptanceDecisionSplit,segment_param);

//                for (int c=0; c<numbchild; c++) {
//                    (TreeSplit->GetChild(c))->PutAllLeavesToChildrenLevel();
//                }
                TreeSplit->PutAllLeavesToMainNodesChildrenLevel();
            }

            TestNormResp=TreeSplit->AreNormRespValid();
            TestNormWeights=TreeSplit->AreWeightsValid();
            if(!TestNormResp){
                cout << "Norm resp non valid for TreeSplit"<<endl;
            }
            if(!TestNormWeights){
                cout << "Norm weights non valid for TreeSplit"<<endl;
            }


            if (AcceptanceDecisionSplit){
                SeenSplit=0;
                SeenMerge=0;
            }
            else{
                SeenSplit++;
                cout<< SeenSplit<< " split operations observed not accepted "<<endl;
            }


            int MergeSize=MergeTest.size();
//            for(int i=0;i<MergeSize;i++){
//                delete MergeTest[i];
//                MergeTest[i]=NULL;
//            }
            MergeTest.clear();
            if(segment_param->SMOrder)
            MergeTest=TreeSplit->GetMergeMoreVertical(c);
            else{
                MergeTest=TreeSplit->OrderingMergingLeaves();
            }
            MergeSize=MergeTest.size();
            vector<MergeKLD*> Merge_tmp;
            for(int i=c*SeenMerge;i<MergeSize;i++){
                Merge_tmp.push_back(MergeTest[i]->CopyMergeKLD());
            }
//            for(int i=0;i<MergeSize;i++){
//                if(MergeTest[i]!=NULL){
//                delete MergeTest[i];
//                MergeTest[i]=NULL;
//                }
//            }
            MergeTest.clear();
            MergeTest=Merge_tmp;
            Merge_tmp.clear();
            MergeSize=MergeTest.size();
            cout<<MergeSize<<" merge operations to test "<<endl;
            if(MergeSize>0){
                vector<MergeKLD*>MergeTested(MergeTest.begin(),MergeTest.begin()+c);
                TreeMerge=TreeSplit->RunMergeOperation(MergeTested,AcceptanceDecisionMerge,segment_param);
            }
            else{
                TreeMerge=TreeSplit;
            }

            if(AcceptanceDecisionMerge){
                SeenMerge=0;
                SeenSplit=0;
            }
            else{
                SeenMerge++;
                cout<<SeenMerge<<" operations not accepted"<<endl;
            }
            TestNormResp=TreeMerge->AreNormRespValid();
            TestNormWeights=TreeMerge->AreWeightsValid();

            if(!TestNormResp){
                cout << "Norm resp non valid for TreeMerge"<<endl;
            }
            if(!TestNormWeights){
                cout << "Norm weights non valid for TreeMerge"<<endl;
            }
            TreeSplit=TreeMerge;

            // Recalculating the lists of Split and Merge
            int SplitSize=SplitTest.size();
//            for(int i=0;i<SplitSize;i++){
//                if(SplitTest[i]!=NULL){
//                delete SplitTest[i];
//                SplitTest[i]=NULL;
//                }
//            }
            SplitTest.clear();
            if(segment_param->SMOrder){
            SplitTest=TreeSplit->GetSplitMoreVertical(c,segment_param);
            }
            else{
                SplitTest=TreeSplit->OrderingSplittingLeaves(segment_param);
            }
            SplitSize=SplitTest.size();
            vector<SplitKLD*> Split_tmp;
            for(int i=c*SeenSplit;i<SplitSize;i++){
//                cout<<i<<" ";
                SplitKLD * SplitCopied=SplitTest[i]->CopySplitKLD();
                Split_tmp.push_back(SplitCopied);
            }
//            for(int i=0;i<SplitSize;i++){
//                cout<< i <<endl;
//                if(SplitTest[i]!=NULL){
//                    if(SplitTest[i]->KLD!=NULL){
//                delete SplitTest[i];
//                    }
//                SplitTest[i]=NULL;
//                }
//            }
            SplitTest.clear();
            SplitTest=Split_tmp;

            MergeSize=MergeTest.size();
//            for(int i=0;i<MergeSize;i++){
//                if(MergeTest[i]!=NULL){
//                delete MergeTest[i];
//                MergeTest[i]=NULL;
//                }
//            }
            MergeTest.clear();
            if(segment_param->SMOrder){
            MergeTest=TreeSplit->GetMergeMoreVertical(c);
            }
            else{
                MergeTest=TreeSplit->OrderingMergingLeaves();
            }
            MergeSize=MergeTest.size();
//            vector<MergeKLD*> Merge_tmp;
            for(int i=c*SeenMerge;i<MergeSize;i++){
                Merge_tmp.push_back(MergeTest[i]->CopyMergeKLD());
            }

            MergeTest.clear();
            MergeTest=Merge_tmp;
            MergeSize=MergeTest.size();
            SplitSize=SplitTest.size();




//TreeSplit->SaveAllClasses("/Users/Carole/Documents/PhD/TmpResultPK2.nii.gz",segment_param);
//TreeSplit->SavePriorsAdapted(segment_param);
//           TreeSplit->SaveTreeInTextFile("/Users/Carole/Documents/PhD/DataTreePK2.txt",segment_param);



//            SplitTest=TreeSplit->GetSplitMoreVertical(3);
//            MergeTest=TreeMerge->GetToMerge();
            numberAllLeaves=TreeSplit->GetNumberAllLeaves();
            numbchild=TreeMerge->GetNumberChildren();
            cout<< "The weights of the general classes are : ";
            for (int c=0;c<numbchild;c++){
                cout<<TreeMerge->GetChild(c)->GetNormWeight()<< "    ******      ";
            }
            cout<<endl;
//            AcceptanceDecision=0;
        }
//        if(MergeTest!=NULL){
//            delete MergeTest;
//            MergeTest=NULL;
//        }


    }
    }
    TreeEM * TreeResult=TreeSplit->CopyTree(NULL);
    //    if(this!=NULL){
    //        delete this;
    //    }
    if(TreeSplit!=NULL){
        delete TreeSplit;
        TreeSplit=NULL;
    }
    cout<<CountRunEM<<" have been performed "<<endl;
    return TreeResult;
}


/**************  TREE RELATED METHODS ******************/


/*A root element is defined by the fact that the pointer to Parent is NULL. FindRoot recursively finds the element whose Parent is NULL. Should be the only one to use links upwards*/
TreeEM * TreeEM::FindRoot(){
    if(this->GetParent()==NULL){
        return this;
    }
    else{
        return this->GetParent()->FindRoot();
    }
}

// If the current node is considered as a mixture but has only one child, the node is replaced by the child itself and transformed as a simple distribution.
void TreeEM::CollapseOnlyChild(){
    if (this->GetNumberChildren()!=1) {
        cout<<"Collapsing only when number of children is 1"<<endl;
        return;
    }
    this->SetParameters(this->GetChild(0)->GetParameters());
    delete this->GetChild(0);
    this->Children.erase(this->Children.begin());

}

// Recursively collapsing all only children in the tree to avoid to have only one child
void TreeEM::CollapseOnlyChildTot(){
    if (this->GetNumberChildren()==1) {
        this->CollapseOnlyChild();
        return;
    }
    else{
        int numbchild=this->GetNumberChildren();
        for (int c=0; c<numbchild; c++) {
            this->GetChild(c)->CollapseOnlyChildTot();
        }
    }
    return;
}

void TreeEM::CreateAndAddChildPriors(SEG_PARAMETERS * segment_param,nifti_image * PriorsInput=NULL,int DistributionType=1 ){

    // First taking care of all needed changes on the created child to make a coherent tree
    TreeEM * ChildCreated=new TreeEM();
    //cout<<ChildCreated->GetDataDirect()<<endl;
    ChildCreated->SetParent(this); // set the parent and allocate the proper amount of memory
    //cout<<ChildCreated->GetDataDirect()<<endl;
    //    cout<<"Allocate memory when creating child"<<endl;
    // Allocation of good amount of memory for Distribution, NonNormResp, NormResp
    int numelmasked=ChildCreated->GetNumberMaskedElements();
    int numel = ChildCreated->GetNumberElements();
    if (ChildCreated->NormResp!=NULL) {
        delete [] ChildCreated->NormResp;
        ChildCreated->NormResp=NULL;
    }
    float * NormRespToSet=new float[numelmasked];
    for(int i=0;i<numelmasked;i++){
        NormRespToSet[i]=0;
    }
    ChildCreated->SetNormResp(NormRespToSet);
    if(NormRespToSet!=NULL){
        delete [] NormRespToSet;
        NormRespToSet=NULL;
    }
    //cout<< "NormResp set";

    // Handling the Priors input
    ChildCreated->SetPriors(PriorsInput); // Takes care of the validity checking of the priors and of making usable as probabilities

    // If priors used and not null natural initialisation of the adapted priors with the priors themselves
    if(ChildCreated->GetPriorsDirect()!=NULL){
        float * PriorsAdaptedToSet=new float[numel];
        float * Priors_PTR=static_cast<float*>(ChildCreated->GetPriors()->data);
        int CountPriorsZero=0;
        for(int i=0;i<numel;i++,Priors_PTR++){
            PriorsAdaptedToSet[i]=*Priors_PTR;
            if(PriorsAdaptedToSet[i]==0){
                CountPriorsZero++;
            }
        }
        Priors_PTR=static_cast<float*>(ChildCreated->GetPriors()->data);
//        SaveTmpResult(PriorsAdaptedToSet,"/Users/Carole/Documents/PhD/PriorsInit.nii.gz");
        ChildCreated->SetPriorsAdapted(PriorsAdaptedToSet);
        if(PriorsAdaptedToSet!=NULL){
            delete [] PriorsAdaptedToSet;
            PriorsAdaptedToSet=NULL;
        }
    }

    // If leaf, and flag_MRF ON, allocate memory for MRF
    if(ChildCreated->IsLeaf()&&segment_param->flag_MRF){
        float * MRFToSet=new float[numelmasked];
        for(int i=0;i<numelmasked;i++){
            MRFToSet[i]=0;
        }
        ChildCreated->SetMRF(MRFToSet);
        if(MRFToSet !=NULL){
            delete [] MRFToSet;
            MRFToSet=NULL;
        }
    }
    // Handling the distribution type input. Must necessary be a simple distribution. By default Gaussian distribution chosen
    if (DistributionType==0){ // Cannot be a mixture if leaf
        DistributionType=1;
    }
    //cout<<"Normally the distribution type is now "<<DistributionType;
    ChildCreated->CreateAllocateAndInitializeParameters(DistributionType);
    //    cout<<"Parameters allocated for child"<<endl;
    // Then taking care of changes needed on parent on which ChildCreated is added
    this->MakeParametersMixture(); // Modify parameters structure : check if it is already a mixture otherwise make it a mixture parameters
    this->Children.push_back(ChildCreated);
}

void TreeEM::CreateAndAddChildWeight(SEG_PARAMETERS * segment_param,float Weight=0,int DistributionType=1 ){

    // First taking care of all needed changes on the created child to make a coherent tree
    TreeEM * ChildCreated=new TreeEM();
    //cout<<ChildCreated->GetDataDirect()<<endl;
    ChildCreated->SetParent(this); // set the parent and allocate the proper amount of memory
    //cout<<ChildCreated->GetDataDirect()<<endl;
    //    cout<<"Allocate memory when creating child"<<endl;
    // Allocation of good amount of memory for Distribution, NonNormResp, NormResp
    int numelmasked=ChildCreated->GetNumberMaskedElements();
    int numel = ChildCreated->GetNumberElements();
    if (ChildCreated->NormResp!=NULL) {
        delete [] ChildCreated->NormResp;
        ChildCreated->NormResp=NULL;
    }
    float * NormRespToSet=new float[numelmasked];
    for(int i=0;i<numelmasked;i++){
        NormRespToSet[i]=0;
    }
    ChildCreated->SetNormResp(NormRespToSet);
    if(NormRespToSet!=NULL){
        delete [] NormRespToSet;
        NormRespToSet=NULL;
    }
    //cout<< "NormResp set";

    // Handling the Priors input
    ChildCreated->SetPriors(NULL); // Takes care of the validity checking of the priors and of making usable as probabilities
    // Handling the weight input
    ChildCreated->SetNormWeight(Weight);

    // If priors used and not null natural initialisation of the adapted priors with the priors themselves
    if(ChildCreated->GetPriorsDirect()!=NULL){
        float * PriorsAdaptedToSet=new float[numel];
        float * Priors_PTR=static_cast<float*>(ChildCreated->GetPriors()->data);
        for(int i=0;i<numel;i++,Priors_PTR++){
            PriorsAdaptedToSet[i]=*Priors_PTR;
        }
        Priors_PTR=static_cast<float*>(ChildCreated->GetPriors()->data);
        ChildCreated->SetPriorsAdapted(PriorsAdaptedToSet);
        if(PriorsAdaptedToSet!=NULL){
            delete [] PriorsAdaptedToSet;
            PriorsAdaptedToSet=NULL;
        }
    }

    // If leaf, and flag_MRF ON, allocate memory for MRF
    if(ChildCreated->IsLeaf()&&segment_param->flag_MRF){
        float * MRFToSet=new float[numelmasked];
        for(int i=0;i<numelmasked;i++){
            MRFToSet[i]=0;
        }
        ChildCreated->SetMRF(MRFToSet);
        if(MRFToSet !=NULL){
            delete [] MRFToSet;
            MRFToSet=NULL;
        }
    }
    // Handling the distribution type input. Must necessary be a simple distribution. By default Gaussian distribution chosen
    if (DistributionType==0){ // Cannot be a mixture if leaf
        DistributionType=1;
    }
    //cout<<"Normally the distribution type is now "<<DistributionType;
    ChildCreated->CreateAllocateAndInitializeParameters(DistributionType);
    //    cout<<"Parameters allocated for child"<<endl;
    // Then taking care of changes needed on parent on which ChildCreated is added
    this->MakeParametersMixture(); // Modify parameters structure : check if it is already a mixture otherwise make it a mixture parameters
    this->Children.push_back(ChildCreated);
}

//void TreeEM::InitialiseBasicTree(SEG_PARAMETERS * segment_param){
//    if (this->IsBasicTree()) {
//        srand (time(NULL));
//        cout<<"This is a basic Tree !!";
//        //vector<nifti_image*> PriorsVector=this->GetPriorsVector();
//        //cout<<this->GetPriorsVector().size();
//        if(ArePriorsNormalised()){
//            if (this->GetPriorsVector()[this->GetNumberChildren()-1]==NULL) {// <=> no priors considered
//                /* If there is no prior considered, the same value is put everywhere, equal 1/number of initial children*/
//                cout<<"No priors set beforehand"<<endl;
//                float * NormResp_PTR=NULL;
//                int numbchild=this->GetNumberChildren();
//                int numbmodal=this->GetNumberModalities();
//                int numelmasked=this->GetNumberMaskedElements();
//                for (int c=0;c<numbchild; c++) {
//                    NormResp_PTR=this->GetChild(c)->GetNormResp();
//                    for (int i=0; i<numelmasked; i++,NormResp_PTR++) {
//                        (*NormResp_PTR)=(float)1.0/numbchild;
//                    }
//                    this->GetChild(c)->SetNormWeight((float)(1.0/(float)this->GetNumberChildren()));

//                    // Random initialisation of the mean, variance set to the identity matrix
//                    if(!this->GetChild(c)->CheckForValidityOfParametersStructure()){
//                        this->GetChild(c)->CreateAllocateAndInitializeParameters(1);
//                    }
//                    for (int m=0; m<numbmodal; m++) {

//                        this->GetChild(c)->GetMean()[m]=(float)((float)rand()/(float)RAND_MAX);
//                        //cout<<"Address mean "<<this->GetChild(c)->GetMean()<<" and "<<this->GetChild(c)->ParametersDistribution->ValueParameters;
//                        //cout<<"Value mean = "<<this->GetChild(c)->GetMean()[m]<<endl;
//                        this->GetChild(c)->GetVariance()[m+numbmodal*m]=0.1;
//                    }
//                }
//                //               this->UpdateDistribution();

//                this->UpdateNonNormResp(segment_param);
//                this->UpdateNormResp();
//                this->UpdateNormRespRoot();
//                if(segment_param->flag_MRF){
//                    float * GMatrixToSet;
//                    if(segment_param->flag_GMatrixIn){
//                        GMatrixToSet =this->PrepareGMatrixFromFile(segment_param->filename_GMatrix);
//                        if(GMatrixToSet==NULL){
//                            segment_param->flag_GMatrix=0;
//                        }
//                    }
//                    else {
//                        GMatrixToSet=this->MRFOptSolveLS();
//                    }

//                    this->SetGMatrix(GMatrixToSet);
//                    if(GMatrixToSet!=NULL){
//                        delete [] GMatrixToSet;
//                        GMatrixToSet=NULL;
//                    }
//                }
//                if(segment_param->flag_MRF){
//                    this->UpdateMRF();
//                }
//                this->UpdateNonNormWeights();
//                //
//                this->UpdateNormWeights();
//                if (BFFlag) {
//                    this->UpdateBFCoeffs();
//                    //                    this->UpdateBFCorrection();
//                    this->UpdateDataBFCorrected();
//                }
//            }
//            else {// Priors are normalised as well as priors adapted and point to valid nifti_images
//                int numel=this->GetNumberElements();
//                int numbchild=this->GetNumberChildren();
//                int numelmasked=this->GetNumberMaskedElements();


////                float * Priors_PTR=NULL;
//                PrecisionTYPE tmpNormWeight;


//                for (int c=0; c<numbchild; c++) {
//                    //cout<<"Initialisation of the pointers"<<endl;
//                    this->GetChild(c)->SetNormWeight(0);
//                    float * RootNormResp_PTR=this->FindRoot()->GetNormResp();
////                    float * NormResp_PTR=NULL;
//                    tmpNormWeight=0;
//                    float * NormResp_PTR=this->GetChild(c)->GetNormResp();
//                    float * NormRespToSet=new float[numelmasked];
//                    for(int i=0;i<numelmasked;i++){
//                        NormRespToSet[i]=0;
//                    }
//                    float * NormRespToSet_PTR=NormRespToSet;
//                    //cout<<"the first value of NormResp is "<<(*NormResp_PTR)<<endl;
//                    float * Priors_PTR=static_cast<float *>(this->GetPriorsVector()[c]->data);
//                    float * Priorsb=Priors_PTR;
//                    int* L2S_PTR=this->GetL2S();
//                    for (int i=0; i<numel; i++, L2S_PTR++,Priors_PTR++) {
//                        if ((*L2S_PTR)>=0) {
//                            //cout<< "the index is "<<(*L2S_PTR)<< " and Priors are "<<(*Priors_PTR)<<endl;
//                            (*NormResp_PTR)=(*Priors_PTR);
//                            *NormRespToSet_PTR=*Priors_PTR;
//                            tmpNormWeight+=(PrecisionTYPE)(*Priors_PTR);
//                            *RootNormResp_PTR+=*Priors_PTR;
//                            NormResp_PTR++;
//                            NormRespToSet_PTR++;
//                            RootNormResp_PTR++;
//                        }
//                    }
//                    this->GetChild(c)->SetNormResp(NormRespToSet);
//                    delete[] NormRespToSet;
//                    NormRespToSet=NULL;
//                    this->GetChild(c)->SetNormWeight(tmpNormWeight/numelmasked);
////                    this->GetChild(c)->WithPriors();

//                    // check if value of Priors is the same as norm resp at each numelmasked value
//                    NormResp_PTR=this->GetChild(c)->GetNormResp();
//                    //cout<<"the first value of NormResp is "<<(*NormResp_PTR)<<endl;
//                    Priors_PTR=static_cast<float *>(this->GetPriorsVector()[c]->data);
////                    float * Priors_PTRb=static_cast<float *>(this->GetChild(c)->GetPriors()->data);
//                    L2S_PTR=this->GetL2S();
//                    int CountDiff=0;
//                    int CountDiffb=0;
//                    int j=0;
//                    for(int i=0;i<numel;i++,Priors_PTR++,L2S_PTR++){
//                        if(*L2S_PTR>=0){
//                            if(*Priors_PTR!=*NormResp_PTR){
//                                CountDiff++;
////                                cout<<*Priors_PTR-*NormResp_PTR<<" ";
//                            }
//                            j++;
////                            if(*Priors_PTRb!=*NormResp_PTR){
////                                CountDiffb++;
//////                                cout<<*Priors_PTRb-*NormResp_PTR<<" ";
////                            }
//                            NormResp_PTR++;
//                        }
//                    }
//                    cout << "CountDiff is "<<CountDiff<<endl;


//                    // check if value of Priors is the same as norm resp at each numelmasked value
//                    NormResp_PTR=this->GetChild(0)->GetNormResp();
//                    //cout<<"the first value of NormResp is "<<(*NormResp_PTR)<<endl;
//                    Priors_PTR=static_cast<float *>(this->GetPriorsVector()[0]->data);
////                    float * Priors_PTRb=static_cast<float *>(this->GetChild(c)->GetPriors()->data);
//                    L2S_PTR=this->GetL2S();
//                     CountDiff=0;
//                     CountDiffb=0;
//                     j=0;
//                    for(int i=0;i<numel;i++,Priors_PTR++,L2S_PTR++){
//                        if(*L2S_PTR>=0){
//                            if(*Priors_PTR!=*NormResp_PTR){
//                                CountDiff++;
////                                cout<<*Priors_PTR-*NormResp_PTR<<" ";
//                            }
//                            j++;
////                            if(*Priors_PTRb!=*NormResp_PTR){
////                                CountDiffb++;
//////                                cout<<*Priors_PTRb-*NormResp_PTR<<" ";
////                            }
//                            NormResp_PTR++;
//                        }
//                    }
//                    cout << "CountDiff is "<<CountDiff<<endl;

//                    this->UpdateParameters();

//                    //cout<<"Initialisation done for this prior"<<endl;
//                    // In basic tree there is always more than one child so no division by 0;
//                }

//                for(int c=0;c<numbchild;c++){
//                                    // check if value of Priors is the same as norm resp at each numelmasked value
//                                    float * NormResp_PTR=this->GetChild(c)->GetNormResp();
//                                    //cout<<"the first value of NormResp is "<<(*NormResp_PTR)<<endl;
//                                    float * Priors_PTR=static_cast<float *>(this->GetChild(c)->GetPriors()->data);
//                                    int * L2S_PTR=this->GetL2S();
//                                    int CountDiff=0;
//                                    for(int i=0;i<numel;i++,Priors_PTR++,L2S_PTR++){
//                                        if(*L2S_PTR>=0){
//                                            if(*Priors_PTR!=*NormResp_PTR){
//                                                CountDiff++;
////                                                cout<<*Priors_PTR-*NormResp_PTR<<" ";
//                                            }
//                                            NormResp_PTR++;
//                                        }
//                                    }
//                                    cout<< "Count Diff is "<<CountDiff<<"for child "<<c<<endl;
//                }


//                bool test=this->IsNormRespValidLeaves();
//                bool testb=this->IsNormRespValidGeneral();
//                bool test2=this->ArePriorsNormalised();
//                if(segment_param->flag_MRF){
//                    float * GMatrixToSet;
//                    if(segment_param->flag_GMatrixIn){
//                        GMatrixToSet =this->PrepareGMatrixFromFile(segment_param->filename_GMatrix);
//                        if(GMatrixToSet==NULL){
//                            segment_param->flag_GMatrix=0;
//                        }
//                    }
//                    else {
//                        GMatrixToSet=this->MRFOptSolveLS();
//                    }

//                    this->SetGMatrix(GMatrixToSet);
//                    if(GMatrixToSet!=NULL){
//                        delete [] GMatrixToSet;
//                        GMatrixToSet=NULL;
//                    }
//                    this->UpdateMRF();
//                }
//            }
//        }
//        else{ // In case the priors are not normalised, then need to normalise them
//            this->NormalisePriors(); // We normalise and call back InitialiseBasicTree
//            this->NormalisePriorsAdapted();
//            cout<<this->ArePriorsNormalised();
//            cout<<this->ArePriorsAdaptedNormalised();
//            cout<<"Normalised Priors"<<endl;
//            this->InitialiseBasicTree(segment_param);
//        }
//    }
//    else{
//        cout<<"This tree is not to initialise"<<endl;
//    }
//    return;
//}


void TreeEM::InitialiseBasicTree(SEG_PARAMETERS * segment_param){
    BFFlag=segment_param->flag_Bias;
    if (this->IsBasicTree()) {
        srand (time(NULL));
        cout<<"This is a basic Tree !!";
        //vector<nifti_image*> PriorsVector=this->GetPriorsVector();
        //cout<<this->GetPriorsVector().size();
        if(ArePriorsNormalised()){
//            bool TestWeight=this->TestInitWeights();
            this->InitialiseNormResp();
//            this->SaveTmpResultMasked(this->GetChild(0)->NormResp,"/Users/Carole/Documents/PhD/NormRespInit.nii.gz");
            this->UpdateNormResp();
            this->UpdateNormRespRoot();
            this->UpdateNonNormWeights();
            this->UpdateNormWeights();
            if(this->GetPriorsVector().size()==0){ // Need for parameters random creation if no priors are provided
                vector<TreeEM*> GeneralNodesVector;
                if(this->FlagOutliers==1){
                    GeneralNodesVector=this->GetChild(0)->GetChildren();
                }
                else{
                    GeneralNodesVector=this->GetChildren();
                }
                int numbmodal=this->GetNumberModalities();
                int numbclasses=this->GetNumberGeneralClasses();
                for(int c=0;c<numbclasses;c++){
                    // Random initialisation of the mean, variance set to the identity matrix
                    if(!GeneralNodesVector[c]->CheckForValidityOfParametersStructure()){
                        GeneralNodesVector[c]->CreateAllocateAndInitializeParameters(1);
                    }
                    for (int m=0; m<numbmodal; m++) {
                        this->GetChild(c)->GetMean()[m]=(float)((float)rand()/(float)RAND_MAX);
                        this->GetChild(c)->GetVariance()[m+numbmodal*m]=0.1;
                    }
                }

                this->UpdateNonNormResp(segment_param);
                this->UpdateNormResp();
                this->UpdateNormRespRoot();
                this->UpdateNonNormWeights();
                this->UpdateNormWeights();
            }
            else { // case we have initialised everything with the priors
                this->UpdateParameters();
            }
                if(segment_param->flag_MRF){
                    float * GMatrixToSet;
                    if(segment_param->flag_GMatrixIn){
                        GMatrixToSet =this->PrepareGMatrixFromFile(segment_param->filename_GMatrix);
                        if(GMatrixToSet==NULL){
                            segment_param->flag_GMatrix=0;
                        }
                    }
                    else {
                        GMatrixToSet=this->MRFOptSolveLS();
                    }

                    this->SetGMatrix(GMatrixToSet);
                    if(GMatrixToSet!=NULL){
                        delete [] GMatrixToSet;
                        GMatrixToSet=NULL;
                    }
                }
                if(segment_param->flag_MRF){
                    this->UpdateMRF();
                }
                this->UpdateNonNormWeights();
                //
                this->UpdateNormWeights();
                if (BFFlag) {
                    this->UpdateBFCoeffs();
                    //                    this->UpdateBFCorrection();
                    this->UpdateDataBFCorrected();
                }
            }

        else{ // In case the priors are not normalised, then need to normalise them
            this->NormalisePriors(); // We normalise and call back InitialiseBasicTree
            this->NormalisePriorsAdapted();
            cout<<this->ArePriorsNormalised();
            cout<<this->ArePriorsAdaptedNormalised();
            cout<<"Normalised Priors"<<endl;
            this->InitialiseBasicTree(segment_param);
        }
    }
    else{
        cout<<"This tree is not to initialise"<<endl;
    }
    return;
}







/* Reinitialise the checks for the SM operations.*/
void TreeEM::ClearSMChecks(){
    //    int numbDirectLeaves=this->GetNumberDirectLeaves();
    int numbchild=this->GetNumberChildren();

    if (numbchild!=0) {
        if (this->GetSplitCheck()!=NULL) {
            delete [] this->SplitCheck;
        }
        this->SplitCheck=NULL;
        this->SplitCheck=new bool[numbchild];//{0};
        for(int c=0;c<numbchild;c++){
            this->SplitCheck[c]=0;
            if(!this->GetChild(c)->IsLeaf()){ // If the considered child is not a leaf it won't be split therefore set as already checked
                this->SplitCheck[c]=1;
            }
        }
        if (this->GetMergeCheck()!=NULL) {
            delete [] this->MergeCheck;
        }
        this->MergeCheck=NULL;
        this->MergeCheck=new bool[numbchild*numbchild];//{0};
        for(int i=0;i<numbchild*numbchild;i++){
            this->MergeCheck[i]=0;
        }
        for (int dl=0; dl<numbchild; dl++) {
            this->MergeCheck[dl+dl*numbchild]=1; // We cannot authorize a merging of the same class with itself !!!
        }
        for(int c1=0;c1<numbchild;c1++){
            if(!this->GetChild(c1)->IsLeaf()){
                for(int c2=0;c2<numbchild;c2++){
                    this->MergeCheck[c1+c2*numbchild]=1; // If one of the child considered is not a leaf, merging cannot be tested
                }

            }
            for(int c2=0;c2<numbchild;c2++){
                if(this->GetChild(c1)->GetDistributionType()!=this->GetChild(c2)->GetDistributionType()){ // cannot accept merging between classes without same distribution
                    this->MergeCheck[c1+c2*numbchild]=1;
                }
            }
        }
    }
    else{
        if (this->GetMergeCheck()!=NULL) {
            delete [] this->MergeCheck;
        }
        this->MergeCheck=NULL;
        if (this->GetSplitCheck()!=NULL) {
            delete [] this->SplitCheck;
        }
        this->SplitCheck=NULL;
    }
    for (int c=0; c<numbchild; c++) {
        this->GetChild(c)->ClearSMChecks();
    }
    return;
}

// Modify NormWeights in the tree in order to just afterwards collapse everything to the children level
void TreeEM::ModifyNormWeightsForChildrenLevel(){
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        if (this->GetPriorsDirect()==NULL&& !this->GetParent()->IsRoot()) {
            this->GetChild(c)->SetNormWeight(this->GetNormWeight()*this->GetChild(c)->GetNormWeight());
        }
        //        else{
        //            this->GetChild(c)->SetPriors(this->GetPriorsDirect());
        //        }
        this->GetChild(c)->ModifyNormWeightsForChildrenLevel(); // Recursive part
    }
}

void TreeEM::ModifyNormWeightsForMainNodesChildrenLevel(){
    TreeEM * MainNode=this->FindMainNode();
    int numbchild=this->GetNumberChildren();
    for(int c=0;c<numbchild;c++){
        if(this->GetPriorsDirect()==NULL && this!=MainNode && MainNode!=NULL){
            this->GetChild(c)->SetNormWeight(this->GetNormWeight()*this->GetChild(c)->GetNormWeight());
        }
        this->GetChild(c)->ModifyNormWeightsForMainNodesChildrenLevel();
    }
}

// Using the preliminary function ModifyNormWeights for ChildrenLevel then collect all the leaves, copy them, delete the children and put the copied leaves instead. Allows for changing the value of parent for the copied leaves right away before adding them as children
void TreeEM::PutAllLeavesToChildrenLevel(){
    this->ModifyNormWeightsForChildrenLevel();
    vector<TreeEM *> Leaves=this->GetAllLeaves();
    int numbLeaves=this->GetNumberAllLeaves();
    vector<TreeEM *> CopiedLeaves;
    // Copy the leaves with the right weights
    for (int l=0; l<numbLeaves; l++) {
        CopiedLeaves.push_back(Leaves[l]->CopyTree(this));
    }
    // Delete the children (not useful anymore)
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        delete this->GetChild(0);
        this->Children.erase(this->Children.begin());
    }
    // Set the copied leaves as the children
    for (int l=0; l<numbLeaves; l++) {
        this->AddChild(CopiedLeaves[l]);
    }
}

void TreeEM::PutAllLeavesToMainNodesChildrenLevel(){
    this->FindRoot()->ModifyNormWeightsForMainNodesChildrenLevel();
    vector<TreeEM *> MainNodes=this->FindRoot()->GetMainNodesVector();
    int numbNodes=MainNodes.size();
    for(int n=0;n<numbNodes;n++){
        int numbLevels=MainNodes[n]->GetNumberLevels();
        if(numbLevels>1){
        vector<TreeEM*> LeavesVectorNode=MainNodes[n]->GetAllLeaves();
        int numbLeaves=LeavesVectorNode.size();
        vector<TreeEM *> CopiedLeaves;
        // Copy the leaves with the right weights
        for (int l=0; l<numbLeaves; l++) {
            CopiedLeaves.push_back(LeavesVectorNode[l]->CopyTree(MainNodes[n]));
        }
        // Delete the children (not useful anymore)
        int numbchild=MainNodes[n]->GetNumberChildren();
        for (int c=0; c<numbchild; c++) {
            delete MainNodes[n]->GetChild(0);
            MainNodes[n]->Children.erase(MainNodes[n]->Children.begin());
        }
        // Set the copied leaves as the children
        for (int l=0; l<numbLeaves; l++) {
            MainNodes[n]->AddChild(CopiedLeaves[l]);
        }
    }
    }
}

void TreeEM::RecPutAllLeavesToChildrenLevel(){
    if (this->GetNumberChildren()==this->GetNumberDirectLeaves()) {
        return;
    }
    else{
        for (int c=0; c<this->GetNumberChildren(); c++) {
            if (this->GetChild(c)->GetNumberChildren()!=0) {
                int numbchild2=this->GetChild(c)->GetNumberChildren();
                for (int c2=0; c2<numbchild2; c2++) {
                    this->GetChild(c)->GetChild(c2)->SetNormWeight(this->GetChild(c)->GetNormWeight()*this->GetChild(c)->GetChild(c2)->GetNormWeight());
                    TreeEM * CopiedChild=this->GetChild(c)->GetChild(c2)->CopyTree(this);
                    this->AddChild(CopiedChild);
                }
                delete this->GetChild(c);
                this->Children.erase(this->Children.begin()+c);
            }
        }
    }
}


/**************** TESTS ON TREE *************************/

/* Test if the tree on which applied is a leaf. A leaf is defined as a tree that does not have any children*/
bool TreeEM::IsLeaf(){
    if(this->GetNumberChildren()==0){
        return 1;
    }
    return 0;
}

/* Test if the tree on which applied is a branch. A branch is defined as a tree that has both a defined Parent and Children*/
bool TreeEM::IsBranch(){
    if(this->GetParent()!=NULL && this->GetNumberChildren()>0){
        return 1;
    }
    return 0;
}

/* Test if the tree on which applied is a root. A root is defined as a tree that has no Parent*/
bool TreeEM::IsRoot(){
    if(this->GetParent()==NULL){
        return 1;
    }
    return 0;
}

/*Returns as an integer the type of component (root,branch,leaf,initialNode) considered. Could be coded afterwards as enum*/
int TreeEM::WhatTypeOfTreeComponent(){
    if (this->IsBranch()) {
        return BRANCH;
    }
    if (this->IsRoot()&&this->IsLeaf()) {
        return INITIALNODE; // is an initial node with no parent and no children
    }
    if (this->IsLeaf()) {
        return LEAF;
    }

    return ROOT; // Is a root
}



/* Tests if the tree considered is a basic tree or not. A basic tree is defined as a tree with only a root and one level of children*/
bool TreeEM::IsBasicTree(){
    int numbchild=this->GetNumberChildren();
    if (this->WhatTypeOfTreeComponent()!=ROOT) {
        return 0;
    }
    if(this->FlagOutliers==1){
        int numbclasses=this->GetChild(0)->GetNumberChildren();
        for (int c=0; c<numbclasses; c++) {
            if (this->GetChild(0)->Children[c]->GetNumberChildren()!=0) {
                return 0;
            }
        }
    }
    else{
    for (int c=0; c<numbchild; c++) {
        if (this->Children[c]->GetNumberChildren()!=0) {
            return 0;
        }
    }
    }
    return 1;
}

bool TreeEM::IsTreeValid(){
    if (this->GetDataImage()==NULL){ // Force the tree to make reference to an image to segment
        cout<< "Not valid because no data image referenced"<<endl;
        return 0;
    }
    if (!this->IsDataImageNormalised()) {
        cout<<"Not valid because data not normalised"<<endl;
        return 0;
    }
    if (this->WhatTypeOfTreeComponent()==LEAF||this->WhatTypeOfTreeComponent()==INITIALNODE){ // Mask is only set at the root level for allocation problems and sizes incompatibilities if changed through the tree
        if(this->GetMaskDirect()!=NULL){
            cout<<"Not valid because Mask defined at wrong place"<<endl;
            return 0;
        }
    }
    if (!this->IsPriorsProbabilityType()) { // Not valid if priors not of probability type
        cout<<"Not valid because Priors not set as probabilities"<<endl;
        return 0;
    }
    if(!this->IsNumbChildOKWithDistType()){ // not valid if number of children not valid with parameters
        cout<<"Not valid of incompatibility between number of children and distribution type"<<endl;
        return 0;
    }
    if(!this->CheckForValidityOfParametersStructure()){ // not valid if parameters structure not valid
        cout<<"Not valid because of unvalidity of Parameters structure"<<endl;
        return 0;
    }
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        if (!this->GetChild(c)->IsTreeValid()) {
            return 0;
        }
    }
    return 1;
}

bool TreeEM::AreWeightsValid(){
    int numbchild=this->GetNumberChildren();
    bool weightsValidity=1;
    float sumWeight=0;
    if(numbchild==0){
        sumWeight=1;
    }
    for(int c=0;c<numbchild;c++){
        sumWeight+=this->GetChild(c)->GetNormWeight();
    }
    if (sumWeight > 1+1E-6 || sumWeight<1-1E-6){
        weightsValidity*= 0;
        return weightsValidity;
    }
    else{
        weightsValidity*=1;
        for(int c=0;c<numbchild;c++){
            weightsValidity*=this->GetChild(c)->AreWeightsValid();
        }
        return weightsValidity;
    }
}

bool TreeEM::IsNormRespValid(){
    float * NormResp_PTR=this->GetNormResp();
    int numelmasked=this->GetNumberMaskedElements();
    for(int i=0;i<numelmasked;i++,NormResp_PTR++){
        if(*NormResp_PTR>1+10E-6){
            cout<< "High Pb"<<endl;
            return 0;
        }
            else if( *NormResp_PTR<=-10E-6){
            cout<<"Negative Pb"<< *NormResp_PTR<< "at index "<<i<<endl;
            return 0;
        }
    }
    return 1;
}

bool TreeEM::IsNormRespValidGeneral(){
    bool CheckResult=1;
//    Check if it is applied to the Root
    if(!this->IsRoot()){
        cout<<"Validity test not applied on proper level"<<endl;
        return 0;
    }
    int numbchild=this->GetNumberChildren();
    int numelmasked=this->GetNumberMaskedElements();
    float * SumCheckNormResp=new float[numelmasked];
    for(int i=0;i<numelmasked;i++){
        SumCheckNormResp[i]=0;
    }
    // Filling of the checking array and checking validity of NormResp values
    float * SumCheckNormResp_PTR=SumCheckNormResp;
    int NbNegativePb=0;
    int NbIndHighPb=0;
    for(int c=0;c<numbchild;c++){
        float * ChildrenNormResp_PTR=this->GetChild(c)->GetNormResp();
        SumCheckNormResp_PTR=SumCheckNormResp;
        for(int i=0;i<numelmasked;i++,SumCheckNormResp_PTR++,ChildrenNormResp_PTR++){
            *SumCheckNormResp_PTR+=*ChildrenNormResp_PTR;
            if(*ChildrenNormResp_PTR<0){
                NbNegativePb++;
            }
            else if(*ChildrenNormResp_PTR>1+10E-6){
                NbIndHighPb++;
            }
        }
    }
    // Checking operation
    SumCheckNormResp_PTR=SumCheckNormResp;
    int NbLowPb=0;
    int NbHighPb=0;
    for(int i=0;i<numelmasked;i++,SumCheckNormResp_PTR++){
        if(*SumCheckNormResp_PTR>1+10E-6){
            NbHighPb++;
            CheckResult=0;
        }
        else if(*SumCheckNormResp_PTR<1-10E-6){
            NbLowPb++;
            CheckResult=0;
        }
    }
    if(NbIndHighPb!=0 || NbNegativePb!=0){
        CheckResult=0;
    }
    delete [] SumCheckNormResp;
    SumCheckNormResp=NULL;
    if(!CheckResult){
    cout<<"HighPb is "<<NbHighPb<<" and LowPb is "<<NbLowPb<< " NegativePb "<<NbNegativePb<<" IndHighPb "<< NbIndHighPb<<endl;
    }
    return CheckResult;
}

bool TreeEM::IsNormRespValidLeaves(){
    bool CheckResult=1;

//    Check if it is applied to the Root
    if(!this->IsRoot()){
        cout<<"Validity test not applied on proper level"<<endl;
        return 0;
    }

    vector<TreeEM*> LeavesVector=this->GetAllLeaves();
    int numbLeaves=LeavesVector.size();
    int numelmasked=this->GetNumberMaskedElements();
    float * SumCheckNormResp=new float[numelmasked];
    for(int i=0;i<numelmasked;i++){
        SumCheckNormResp[i]=0;
    }
    // Filling of the checking array
    float * SumCheckNormResp_PTR=SumCheckNormResp;
    for(int c=0;c<numbLeaves;c++){
        float * ChildrenNormResp_PTR=LeavesVector[c]->GetNormResp();
        SumCheckNormResp_PTR=SumCheckNormResp;
        for(int i=0;i<numelmasked;i++,SumCheckNormResp_PTR++,ChildrenNormResp_PTR++){
            *SumCheckNormResp_PTR+=*ChildrenNormResp_PTR;
        }
    }
    // Checking operation
    int NbLowPb=0;
    int NbHighPb=0;
    SumCheckNormResp_PTR=SumCheckNormResp;
    for(int i=0;i<numelmasked;i++,SumCheckNormResp_PTR++){
        if(*SumCheckNormResp_PTR>1+10E-6){
            NbHighPb++;
            CheckResult=0;
        }
        else if(*SumCheckNormResp_PTR<1-10E-6){
            NbLowPb++;
            CheckResult=0;
        }
    }
    delete [] SumCheckNormResp;
    SumCheckNormResp=NULL;
    if(!CheckResult){
    cout<<"HighPb is "<<NbHighPb<<" and LowPb is "<<NbLowPb<<endl;
    }
    return CheckResult;
}



bool TreeEM::AreNormRespValid(){
    bool NormRespValidity=1;
    int numbchild=this->GetNumberChildren();
    NormRespValidity*=this->IsNormRespValid();
    if(!NormRespValidity){
        return NormRespValidity;
    }
    else{
        for(int c=0;c<numbchild;c++){
            NormRespValidity*=this->GetChild(c)->AreNormRespValid();
        }
        return NormRespValidity;
    }
}

bool TreeEM::IsNumbChildOKWithDistType(){
    if (this->GetNumberChildren()==0&&this->GetDistributionType()==0){ // case considered as mixture whereas does not have any child
        return 0;
    }
    if (this->GetNumberChildren()>0 && this->GetDistributionType()>0){ // case where considered as simple distribution whereas has children
        return 0;
    }
    return 1;
}

bool TreeEM::IsDataImageNormalised(){
    if (this->GetDataImage()==NULL) {// no data image in the tree to look at
        cout<<"No data image as reference"<<endl;
        return 0;
    }
    int numbmodal=this->GetNumberModalities();
    for (int m=0; m<numbmodal; m++) {
        //cout << "The number of modalities is "<<this->GetNumberModalities()<<endl;
        float minRes=this->GetMinDataModal(m);
        float maxRes=this->GetMaxDataModal(m);
        //cout<<(minRes-0.0)<<"and "<<(maxRes-1.0)<<endl;
        if (minRes>1E-6||minRes<-1E-6){

            return 0;
        }
        if (maxRes>1+1E-6||maxRes<1-1E-6){
            return 0;
        }
    }
    return 1;
}

bool TreeEM::IsDataFloat(){
    if (this->GetDataImage()==NULL) {
        cout<<"No data to check"<<endl;
        return 0;
    }
    if (this->GetDataImage()->datatype==DT_FLOAT) {
        return 1;
    }
    return 0;
}

bool TreeEM::IsPriorsFloat(){
    if (this->GetPriors()==NULL) {
        cout<<"No priors to check"<<endl;
        return 0;
    }
    if (this->GetPriors()->datatype==DT_FLOAT) {
        return 1;
    }
    return 0;
}

bool TreeEM::IsStructureSimilar(TreeEM *TreeTest){
    bool Result=1;
    if(this->GetNumberChildren()==0 && TreeTest->GetNumberChildren()==0){
        return Result*1;
    }
    else if( this->GetNumberChildren()!=TreeTest->GetNumberChildren()){
        return Result*0;
    }
    else{
        int numbchild=this->GetNumberChildren();
        for(int c=0;c<numbchild;c++){
            Result*=this->GetChild(c)->IsStructureSimilar(TreeTest->GetChild(c));
        }
    }
    return Result;
}

// Check if the direct BF coefficients obtained are directly useable or not
bool TreeEM::AreBFCoeffsDirectUseable(){
    //    vector<float*> BFCoeffsToTest=this->GetBFCoeffsDirect();
    //    // 1. Test on the size of the vector
    //    if (BFCoeffsToTest.size()==0) {
    //        return 0;
    //    }
    //    int numbmodal=this->GetNumberModalities();
    //    // 2. Test on the compatibility of the dimensions
    //    if (BFCoeffsToTest.size()!=numbmodal) {
    //        return 0;
    //    }
    //    // 3. Test on the containts of the vector if of the right size.
    //    for (int m=0; m<numbmodal; m++) {
    //        if (BFCoeffsToTest[m]==NULL) {
    //            return 0;
    //        }
    //    }
    if(this->GetBFCoeffsDirect()==NULL){
        return 0;
    }
    return 1;
}

bool TreeEM::IsMRFZero(){
    float * MRFToTest=this->GetMRF();
    int numelmasked=this->GetNumberMaskedElements();
    if(MRFToTest==NULL){
        return 0;
    }
    for(int i=0;i<numelmasked;i++){
        if(MRFToTest[i]>10E-6){
            return 0;
        }
    }
    return 1;
}

/********************** METHODS ON PARAMETERS STRUCTURE ***********************/

int TreeEM::CalculateSizeParameters(){
    int SizeParameters;
    switch (this->GetDistributionType()) {
    case 0:{
        //cout << "We are in the case of a mixture" << endl;
        SizeParameters=0;
    }
        break;
    case 2:{ // case uniform distribution no parameters but can be a leaf
        SizeParameters=0;
    }
        break;
    default:
    {
        //cout<< "We are looking at a Gaussian distribution" << endl;
        int numbmodal=this->GetNumberModalities();
        SizeParameters=(int)(numbmodal*(numbmodal+1)); // Size of the mean vector + size of the variance matrix (numbmodal)
        //cout<< "Size parameters is "<<SizeParameters;
        /* REMARK :
             this might be changed if decided then to consider only number of free parameters, then the Gaussian Parameters calculation function must be changed consequently */
    }

    }
    return SizeParameters;

}

int TreeEM::CalculateSizeParameters(int DistributionTypeInput){
    int SizeParameters;

    switch (DistributionTypeInput) {
    case 0:{
        //cout << "We are in the case of a mixture" << endl;
        SizeParameters=0;
    }
        break;
    case 2:{ // case of uniform distribution without any parameters
        SizeParameters=0;
    }
        break;
    default:
    {
        //cout<< "We are looking at a Gaussian distribution" << endl;
        int numbmodal=this->GetNumberModalities();
        SizeParameters=(int)(numbmodal*(numbmodal+1)); // Size of the mean vector + size of the variance matrix (numbmodal)
        //cout<< "Size parameters is "<<SizeParameters;
        /* REMARK :
             this might be changed if decided then to consider only number of free parameters, then the Gaussian Parameters calculation function must be changed consequently */
    }

    }
    return SizeParameters;
}

void TreeEM::CreateAllocateAndInitializeParameters(int DistributionTypeInput){
    if (this->GetNumberChildren()!=0) { // if there are children, it is necessarily a mixture
        cout<<"We have to modify Distribution type"<<endl;
        DistributionTypeInput=0;
    }
    if (this->GetParameters()==NULL) {
        //cout<<"Parameters are NULL at the momemt"<<endl;
        this->ParametersDistribution=new Parameters();
        this->ParametersDistribution->DistributionType=DistributionTypeInput;
        //cout<<this->GetDistributionType()<< "and " << this->ParametersDistribution->DistributionType;
        this->ParametersDistribution->SizeParameters=this->CalculateSizeParameters();
        if (this->GetSizeParameters()>0) {
            this->GetParameters()->ValueParameters=new float[this->GetSizeParameters()];//{0};
            for(int i=0;i<this->GetSizeParameters();i++){
                this->GetParameters()->ValueParameters[i]=0;
            }
        }
        else{
            this->GetParameters()->ValueParameters=NULL;
        }
        return;
    }
    this->ParametersDistribution->DistributionType=DistributionTypeInput;
    this->ParametersDistribution->SizeParameters=this->CalculateSizeParameters();
    //cout<<"the size of parameters is "<<this->GetSizeParameters()<<endl;
    if (this->GetSizeParameters()>0) {
        if(this->GetParametersValue()!=NULL){
            delete [] this->GetParametersValue();
            this->GetParameters()->ValueParameters=NULL;
        }
        int numbsp=this->GetSizeParameters();
        this->GetParameters()->ValueParameters=new float[this->GetSizeParameters()];//{0};
        for(int i=0;i<numbsp;i++){
            this->GetParameters()->ValueParameters[i]=0;
        }
    }
    else {
        if(this->GetParametersValue()!=NULL){
            delete [] this->GetParametersValue();
        }
        this->GetParameters()->ValueParameters=NULL;
    }
    return;
}

bool TreeEM::CheckForValidityOfParametersStructure(){
    // Incompatibility between type mixture and existence of values for the parameters

    if(this->GetDistributionType()==0 && this->GetParametersValue()!=NULL){
        return 0;
    }
    // Incompatibility between mixture DistributionType and SizeParameters > 0
    if (this->GetDistributionType()==0 && this->GetSizeParameters()>0) {
        return 0;
    }
    // Compatibility between uniform and no parameters
    if(this->GetDistributionType()==2 && this->GetSizeParameters()==0){
        return 1;
    }
    // Incompatibility between simple distribution and SizeParameters=0
    if (this->GetDistributionType()>0 && this->GetDistributionType()!=2 && this->GetSizeParameters()==0) {
        return 0;
    }
    // Incompatibility between simple distribution and non allocation of values for parameters
    if(this->GetDistributionType()>0 && this->GetParametersValue()==NULL){
        return 0;
    }
    // Incompatibility if SizeParameters is 0 and ValueParameters is allocated
    if (this->GetSizeParameters()==0 && this->GetParametersValue()!=NULL) {
        return 0;
    }
    // Incompatibility if SizeParameters is >0 and ValueParameters is not allocated
    if (this->GetSizeParameters()>0 && this->GetParametersValue()==NULL) {
        return 0;
    }
    else {
        return 1;
    }
}

// Check for the compatibility of the data and the size of the parameters to set
bool TreeEM::CheckForSizeParametersValidity(){
    // First check for the validity of the structure of ParametersToSet
    if(!CheckForValidityOfParametersStructure()){
        return 0;
    }
    // Once checked for the validity of the structure check if mixture type or something else (if mixture, ok by definition)
    if(this->GetDistributionType()==0){
        return 1;
    }
    if (this->GetSizeParameters()!=this->CalculateSizeParameters()){
        return 0;
    }
    return 1;
}

bool TreeEM::ArePriorsNormalised(){
    int numel=this->GetNumberElements();
    if (this->IsPriorsVectorFilledWithNULL()) {
        return 1;
    }
//    int numbclasses=this->GetPriorsVector().size();
    int numbchild=this->GetNumberGeneralClasses();
    int numbclasses=this->GetPriorsVector().size();
    // check if there are compatible number of Priors with number of general classes. Otherwise all priors put to NULL;
    if(numbclasses !=numbchild && numbclasses !=0 && this->GetFlagOutliers()==0){ // incompatibility in the number of priors compared to the number of general classes.
        cout<<"Incompatible number of priors files compared to general classes"<<endl;
        return 0;
    }

//    /*Check if some or all of them are NULL pointers in the vector if not all of them, must be normalised (all must be put to NULL)*/
//    int indexFirstNull=numbclasses;
//    for(int c=numbclasses-1;c>=0;c--){
//        if (this->GetPriorNodeVector()[c]==NULL){
//            indexFirstNull=c;
//        }
//    }
//    if(indexFirstNull<numbclasses){// Means one of the pointer is NULL but one at least the first is not NULL
//        cout<<"One of the Pointer to priors is NULL and is not the first one"<<endl;
//        return 0;
//    }
//    if(indexFirstNull==0){// the first pointer is NULL, we must check for all the others in the vector
//        for (int c=0;c<numbclasses;c++){
//            if (this->GetPriorNodeVector()[c]!=NULL){
//                cout<<"One of the Priors is not NULL whereas the first one is NULL"<<endl;
//                return 0;
//            }
//        }
//        return 1; // All the pointers in the vector are NULL
//    }
    /* if get there : all Priors pointers in the vector correspond to valid pointer to nifti images*/
    vector<float *>PriorsVectorData_PTR;
    //cout << "All pointers correspond to a valid nifti_image"<<endl;
    for (int c=0;c<numbclasses;c++){
        PriorsVectorData_PTR.push_back(static_cast<float *>(this->GetPriorsVector()[c]->data));
    }
    PrecisionTYPE sumCheck=0;
    for (int i=0;i<numel;i++){
        sumCheck=0;
        for(int c=0;c<numbclasses;c++){
            sumCheck+=(PrecisionTYPE)*PriorsVectorData_PTR[c];
            PriorsVectorData_PTR[c]++;
        }


        if(fabs(sumCheck-1)>10E-6){
            //cout<<"The priors are not well normalised, the difference to 1 is "<<sumCheck-1<<endl;

            return 0;
        }
    }
    return 1;
}

bool TreeEM::ArePriorsAdaptedNormalised(){
    int numel=this->GetNumberElements();
    if (this->IsPriorsAdaptedVectorFilledWithNULL()) {
        return 1;
    }
    int numbclasses=this->GetPriorsAdaptedVector().size();
    int numbchild=this->GetNumberGeneralClasses(); // Priors can only be adapted when there is anatomical knowledge and appropriate behavior
    if(numbchild != numbclasses && numbclasses !=0 && this->GetFlagOutliers()==0){
        cout<<"Inappropriate number of general classes compared to number of adapted priors arrays"<<endl;
        return 0;
    }
//    /*Check if some or all of them are NULL pointers in the vector if not all of them, must be normalised (all must be put to NULL)*/
//    int indexFirstNull=numbclasses;
//    for(int c=numbclasses-1;c>=0;c--){
//        if (this->GetPriorsAdaptedVector()[c]==NULL){
//            indexFirstNull=c;
//        }
//    }
//    if(indexFirstNull<numbclasses){// Means one of the pointer is NULL but one at least the first is not NULL
//        cout<<"One of the Pointer to priors is NULL and is not the first one"<<endl;
//        return 0;
//    }
//    if(indexFirstNull==0){// the first pointer is NULL, we must check for all the others in the vector
//        for (int c=0;c<numbclasses;c++){
//            if (this->GetPriorsAdaptedVector()[c]!=NULL){
//                cout<<"One of the Priors is not NULL whereas the first one is NULL"<<endl;
//                return 0;
//            }
//        }
//        return 1; // All the pointers in the vector are NULL
//    }
    /* if get there : all Priors pointers in the vector correspond to valid pointer to float array*/
    vector<float *>PriorsVectorData_PTR=this->GetPriorsAdaptedVector();
    PrecisionTYPE sumCheck=0;
    for (int i=0;i<numel;i++){
        sumCheck=0;
        for(int c=0;c<numbclasses;c++){
            sumCheck+=(PrecisionTYPE)*PriorsVectorData_PTR[c];
            PriorsVectorData_PTR[c]++;
        }
        if(fabs(sumCheck-1)>1E-6){
            //cout<<"The priors are not well normalised, the difference to 1 is "<<sumCheck-1<<endl;

            return 0;
        }
    }
    return 1;
}

bool TreeEM::IsPriorsVectorFilledWithNULL(){
//    int numbchild=this->GetNumberChildren();
    vector<nifti_image*> NodePriorsVector=this->GetPriorsVector();
    int numbchild=NodePriorsVector.size();
//    for (int c=0; c<numbchild; c++) {
////        if (this->GetChild(c)->GetPriors()!=NULL) {
////            return 0;
////        }

//    }
    if(numbchild>0){ // If at least one of the priors is not null
        return 0;
    }
    return 1; // By default return 1 if it is a leaf or an initial Node we are looking at
}

bool TreeEM::IsPriorsAdaptedVectorFilledWithNULL(){
    vector<float *> NodePriorsAdaptedVector=this->GetPriorsAdaptedVector();
    int numbchild=NodePriorsAdaptedVector.size();
    if(numbchild>0){
        return 0;
    }
//    int numbchild=this->GetNumberChildren();
//    for (int c=0; c<numbchild; c++) {
//        if (this->GetChild(c)->GetPriorsAdapted()!=NULL) {
//            return 0;
//        }
//    }
    return 1; // By default return 1 if it is a leaf or an initial Node we are looking at
}

bool TreeEM::IsOneOfThePriorsNULL(){
    int numbclasses=this->GetNumberGeneralClasses();
    int numbchild=this->GetPriorsVector().size();
//    int numbchild=this->GetNumberChildren();
    if ((numbclasses !=numbchild && numbclasses !=numbchild-1) && numbchild!=0){ // if the number of general classes is different from the number of priors image to use but at least one of then is not NULL
        return 1;
    }
    return 0;
//    int numbchild=this->GetNumberChildren();
//    for (int c=0; c<numbchild; c++) {
//        if (this->GetChild(c)->GetPriors()==NULL) {
//            return 1;
//        }
//    }
//    return 0;
}

bool TreeEM::IsOneOfThePriorsAdaptedNULL(){
    int numbclasses=this->GetNumberGeneralClasses();
    int numbchild=this->GetPriorsAdaptedVector().size();
    if (numbclasses !=numbchild && numbclasses !=numbchild-1 && numbchild!=0){ // if the number of general classes is different from the number of priors image to use but at least one of then is not NULL
        return 1;
    }
    return 0;
//    int numbchild=this->GetNumberChildren();
//    for (int c=0; c<numbchild; c++) {
//        if (this->GetChild(c)->GetPriorsAdapted()==NULL) {
//            return 1;
//        }
//    }
//    return 0;
}

void TreeEM::NormalisePriorsAdapted(){
//    int numbchild=this->GetNumberChildren();
    int numbclasses=this->GetNumberGeneralClasses();
    int numel=this->GetNumberElements();

    if(this->ArePriorsAdaptedNormalised()){
        cout<<"Priors are already normalised";
        return;
    }
    else{
        vector<TreeEM*> NodePriorsVector=this->GetPriorsNodeVector();
        int numbchild=NodePriorsVector.size();
        // First if one of them is a pointer to NULL, all of the priors are deleted and set to NULL
        if( this->IsOneOfThePriorsAdaptedNULL()){
//            int numbchild=this->GetNumberChildren();
//            for (int c=0; c<numbchild; c++) {
//                if(this->GetChild(c)->PriorsAdapted!=NULL){
//                    delete [] this->GetChild(c)->PriorsAdapted;
//                    this->GetChild(c)->PriorsAdapted=NULL;
//                }
//            }
            for(int c=0;c<numbchild;c++){
                delete [] NodePriorsVector[c]->GetPriorsAdaptedDirect();
                NodePriorsVector[c]->PriorsAdapted=NULL;
            }
        }
        else{ // Priors adapted are not normalised but all pointing to valid float array
            vector<float *> PriorsData_PTR;
            for (int c=0; c<numbchild; c++) { // Pushing pointers to begin of each of the Priors in the vector PriorsData_PTR
//                if (!this->GetChild(c)->IsPriorsAdaptedProbabilityType()) {
//                    cout<<"Priors is not of probability type";
//                    this->GetChild(c)->MakePriorsAdaptedProbabilityType();
//                    cout << "and now is of probability type "<<this->GetChild(c)->IsPriorsAdaptedProbabilityType();
//                }
//                //cout<< "Priors already of probability type"<<endl;
//                PriorsData_PTR.push_back(this->GetChild(c)->GetPriorsAdapted());
//            }
                if(!NodePriorsVector[c]->IsPriorsAdaptedProbabilityType()){
                    NodePriorsVector[c]->MakePriorsAdaptedProbabilityType();
                }
                PriorsData_PTR.push_back(NodePriorsVector[c]->GetPriorsAdapted());
            }
            int * L2S_PTR=this->GetL2S(); // Pointer to beginning of L2S
            PrecisionTYPE tmpSum=0;
            for (int i=0; i<numel; i++,L2S_PTR++) { // for each of the voxels

                tmpSum=0;
                for (int c=0; c<numbchild; c++) { // calculation of the normalising factor
                    tmpSum+=(PrecisionTYPE)(*PriorsData_PTR[c]);
                }
                for (int c=0; c<numbchild; c++) { // Normalisation of the priors data
                    if(tmpSum>0){
                        (*PriorsData_PTR[c])/=tmpSum;
                    }
                    else{ // taking care of breaking case where normalisation factor is 0
                        //cout<<"Normalisation factor is 0"<<endl;
                        (*PriorsData_PTR[c])=(float)1.0/numbchild;
                    }
                }
                for (int c=0; c<numbchild; c++) { // incrementation of the pointer
                    PriorsData_PTR[c]++;
                }
            }
        }
        return;
    }
}

void TreeEM::NormalisePriors(){
//    int numbchild=this->GetNumberChildren();
//    int numbchild=this->GetNumberGeneralClasses();
//    int numbchild=this->GetPriorsNodeVector().size();
    int numel=this->GetNumberElements();

    if(this->ArePriorsNormalised()){
        cout<<"Priors are already normalised";
        return;
    }
    else{
    vector<TreeEM*> NodeVector=this->GetPriorsNodeVector();
    int numbchild=NodeVector.size();
        // First if one of them is a pointer to NULL, all of the priors are deleted and set to NULL
        if( this->IsOneOfThePriorsNULL()){
            for(int c=0;c<numbchild;c++){
                delete NodeVector[c]->Priors;
                NodeVector[c]->Priors=NULL;
            }

//            int numbchild=this->GetNumberChildren();
//            for (int c=0; c<numbchild; c++) {
//                if(this->GetChild(c)->Priors!=NULL){
//                    delete this->GetChild(c)->Priors;
//                    this->GetChild(c)->Priors=NULL;
//                }
//            }
        }
        else{ // Priors are not normalised but all pointing to valid nifti images
            vector<float *> PriorsData_PTR;
            for (int c=0; c<numbchild; c++) { // Pushing pointers to begin of each of the Priors in the vector PriorsData_PTR
//                if (!this->GetChild(c)->IsPriorsProbabilityType()) {
//                    cout<<"Priors is not of probability type";
//                    this->GetChild(c)->MakePriorsProbabilityType();
//                    cout << "and now is of probability type "<<this->GetChild(c)->IsPriorsProbabilityType();
//                }
//                //cout<< "Priors already of probability type"<<endl;
//                PriorsData_PTR.push_back(static_cast<float*>(this->GetChild(c)->GetPriors()->data));
                if(!NodeVector[c]->IsPriorsProbabilityType()){
                    NodeVector[c]->MakePriorsProbabilityType();
                }
                PriorsData_PTR.push_back(static_cast<float*>(NodeVector[c]->GetPriors()->data));
            }
            int * L2S_PTR=this->GetL2S(); // Pointer to beginning of L2S
            PrecisionTYPE tmpSum=0;
            for (int i=0; i<numel; i++,L2S_PTR++) { // for each of the voxels

                tmpSum=0;
                for (int c=0; c<numbchild; c++) { // calculation of the normalising factor
                    tmpSum+=(PrecisionTYPE)(*PriorsData_PTR[c]);
                }
                for (int c=0; c<numbchild; c++) { // Normalisation of the priors data
                    if(tmpSum>0){
                        (*PriorsData_PTR[c])/=tmpSum;
                    }
                    else{ // taking care of breaking case where normalisation factor is 0
                        //cout<<"Normalisation factor is 0"<<endl;
                        (*PriorsData_PTR[c])=(float)1.0/numbchild;
                    }
                }
                for (int c=0; c<numbchild; c++) { // incrementation of the pointer
                    PriorsData_PTR[c]++;
                }
            }
        }
        this->ArePriorsNormalised();
        return;
    }
}


// Check if the priors is a probability type Prior (everything float between 0 and 1)
bool TreeEM::IsPriorsProbabilityType(){
    if (this->GetPriors()==NULL) {
        return 1;
    }
    else{
        int numel=this->GetNumberElements();
        float * Priors_PTR=static_cast<float*>(this->GetPriors()->data);
        float maxPriors=-1E32;
        float minPriors=1E32;
        for (int i=0; i<numel; i++,Priors_PTR++) {
            if ((*Priors_PTR)>maxPriors) {
                maxPriors=(*Priors_PTR);
            }
            if ((*Priors_PTR)<minPriors) {
                minPriors=(*Priors_PTR);
            }
        }
        bool testminPriors=(minPriors>-1E-6);
        bool testmaxPriors=(maxPriors<1+1E-6);
        if (testmaxPriors&&testminPriors) {
            //cout<<"Truly of probability type"<<endl;
            return 1;
        }
    }
    return 0;
}

// Check if the priors adapted is a probability type Prior (everything float between 0 and 1)
bool TreeEM::IsPriorsAdaptedProbabilityType(){
    if (this->GetPriorsAdapted()==NULL) {
        return 1;
    }
    else{
        int numel=this->GetNumberElements();
        float * Priors_PTR=this->GetPriorsAdapted();
        float maxPriors=-1E32;
        float minPriors=1E32;
        for (int i=0; i<numel; i++,Priors_PTR++) {
            if ((*Priors_PTR)>maxPriors) {
                maxPriors=(*Priors_PTR);
            }
            if ((*Priors_PTR)<minPriors) {
                minPriors=(*Priors_PTR);
            }
        }
        bool testminPriors=(minPriors>-1E-6);
        bool testmaxPriors=(maxPriors<1+1E-6);
        if (testmaxPriors&&testminPriors) {
            //cout<<"Truly of probability type"<<endl;
            return 1;
        }
    }
    return 0;
}


bool TreeEM::CheckForValidityOfParametersStructure(Parameters * ParametersToCheck){
    // Incompatibility between type mixture and existence of values for the parameters
    if(ParametersToCheck->DistributionType==0 && ParametersToCheck->ValueParameters!=NULL){
        return 0;
    }
    // Incompatibility between mixture DistributionType and SizeParameters > 0
    if (ParametersToCheck->DistributionType==0 && ParametersToCheck->SizeParameters>0) {
        return 0;
    }
    // Incompatibility between simple distribution and SizeParameters=0
    if (ParametersToCheck->DistributionType>0 && ParametersToCheck->SizeParameters==0) {
        return 0;
    }
    // Incompatibility between simple distribution and non allocation of values for parameters
    if(ParametersToCheck->DistributionType>0 && ParametersToCheck->ValueParameters==NULL){
        return 0;
    }
    // Incompatibility if SizeParameters is 0 and ValueParameters is allocated
    if (ParametersToCheck->SizeParameters==0 && ParametersToCheck->ValueParameters!=NULL) {
        return 0;
    }
    // Incompatibility if SizeParameters is >0 and ValueParameters is not allocated
    if (ParametersToCheck->SizeParameters>0 && ParametersToCheck->ValueParameters==NULL) {
        return 0;
    }
    else {
        return 1;
    }
}

// Check for the compatibility of the data and the size of the parameters to set
bool TreeEM::CheckForSizeParametersValidity( Parameters * ParametersToSet){
    // First check for the validity of the structure of ParametersToSet
    if(!CheckForValidityOfParametersStructure(ParametersToSet)){
        return 0;
    }
    // Once checked for the validity of the structure check if mixture type or something else (if mixture, ok by definition)
    if(ParametersToSet->DistributionType==0){
        return 1;
    }
    if (ParametersToSet->SizeParameters!=this->CalculateSizeParameters(ParametersToSet->DistributionType)){
        return 0;
    }
    return 1;
}


/***************************** CALCULATION OF KLD ********************************/

// Update the histogram of the data weighted with NormResp for each of the modality in a vector of float pointer to each of the monomodal histogram
vector<float *> TreeEM::GetDataHistogram(){
    vector<float *> DataHistogram_VEC;
    // First check if there is Data to make histogram
    if (this->GetDataImage()==NULL) {
        cout<<"No data to make histogram from"<<endl;
        return DataHistogram_VEC;
    }

    // Then if modality wanted in histogram is available
    int numbbins=this->GetNumbbins();
    if (numbbins<=0) {
        cout<<"No proper number of bins asked for"<<endl;
        return DataHistogram_VEC;
    }
    float sizeBin=(float)1.0/numbbins;
    int *L2S_PTR=this->GetL2S();
//    int numel=this->GetNumberElements();
    int numelmasked=this->GetNumberMaskedElements();
    int numbmodal=this->GetNumberModalities();
    //    int numbchild=this->GetNumberChildren();
    float tmpValue=0;
    //    float * Data=this->MakeDataBFCorrected();
    float * Data=this->GetDataBFCorrected();
    //float * Data_PTR =static_cast<float*>(this->GetDataImage()->data);
    //vector<float *> DataHistogram_VEC;
    for (int m=0; m<numbmodal; m++) {
        float * InitDataHistToPushBack=new float[numbbins];
        for(int i=0;i<numbbins;i++){
            InitDataHistToPushBack[i]=0;
        }
        DataHistogram_VEC.push_back(InitDataHistToPushBack);
    }
    float * NormResp_PTR=this->GetNormResp();
    for (int m=0; m<numbmodal; m++) {
        float * DataHistogram_PTR=DataHistogram_VEC[m];
        //float * DataHistogram=new float[numbbins]{0};

        L2S_PTR=this->GetL2S();

        //Data_PTR=&static_cast<float*>(this->GetDataImage()->data)[m*numel];
        float * Data_PTR=&Data[m*numelmasked];
        NormResp_PTR=this->GetNormResp();
        PrecisionTYPE SumNormResp=0;
        tmpValue=0;
        // Case when NormResp is NULL
        if (this->GetNormResp()==NULL) { // then consider that there are 1 everywhere
            for (int i=0; i<numelmasked; i++, Data_PTR++) {
                tmpValue=(float)*Data_PTR/sizeBin-0.5;
                if (tmpValue<0) {
                    DataHistogram_PTR[0]+=1;
                }
                else if (tmpValue>numbbins-1){
                    DataHistogram_PTR[numbbins-1]+=1;
                }
                else{
                    DataHistogram_PTR[(int)floorf(tmpValue)]+=(float)(1-(tmpValue-floorf(tmpValue)));
                    DataHistogram_PTR[(int)floorf(tmpValue)+1]+=(float)(tmpValue-floorf(tmpValue));
                }
                SumNormResp++;
                //NormResp_PTR++;

            }

        }
        else{
            for (int i=0; i<numelmasked; i++, Data_PTR++,NormResp_PTR++) {
                tmpValue=(float)*Data_PTR/sizeBin-0.5;
                if (tmpValue<0) {
                    DataHistogram_PTR[0]+=(float)*NormResp_PTR;
                }
                else if (tmpValue>numbbins-1){
                    DataHistogram_PTR[numbbins-1]+=(float)*NormResp_PTR;
                }
                else{
                    DataHistogram_PTR[(int)floorf(tmpValue)]+=(float)*NormResp_PTR*(1-(tmpValue-floorf(tmpValue)));
                    DataHistogram_PTR[(int)floorf(tmpValue)+1]+=(float)*NormResp_PTR*(tmpValue-floorf(tmpValue));
                }
                SumNormResp+=(float)*NormResp_PTR;
            }

        }
        DataHistogram_PTR=DataHistogram_VEC[m];
        //int numelmasked=this->GetNumberMaskedElements();
//        float SumData=0;
//        for (int i=0; i<numbbins; i++,DataHistogram_PTR++) {
//            *DataHistogram_PTR/=SumNormResp;
//            SumData+=*DataHistogram_PTR;
//        }
        //cout<<"SumData in GetDataHistogram is"<<SumData;
    }
    //    if (Data!=NULL) {
    //        delete [] Data;
    //        Data=NULL;
    //    }
    return DataHistogram_VEC;
}

/* Update the DataHistogram with all modalities weighted by NormResp.
 At the root, it will be the true histogram of the masked image*/
float * TreeEM::GetDataHistogramTotal(){
    // First check if there is Data to make histogram
    if (this->GetDataImage()==NULL) {
        cout<<"No data to make histogram from"<<endl;
        return NULL;
    }
    // then check that the number of bins set for the tree is compatible with
    int numbbins=this->GetNumbbins();
    if (numbbins<=0) {
        return NULL;
    }

    else{
        int numbmodal=this->GetNumberModalities();
        //        int numbchild=this->GetNumberChildren();
        int numelmasked=this->GetNumberMaskedElements();
        int numelhist=(int)pow_int(numbbins,numbmodal);
        //float * DataHistogram=new float[(int)powf(numbbins,numbmodal)]{0};
        float * DataHistogramTotal_PTR=new float[numelhist];//{0};
        for(int i=0;i<numelhist;i++){
            DataHistogramTotal_PTR[i]=0;
        }
        float sizeBin=(float)1.0/numbbins;
        //        int *L2S_PTR=this->GetL2S();

        vector<float *> Data_PTRModalVector;
        //        float * Data=this->MakeDataBFCorrected();
        float * Data=this->GetDataBFCorrected();

        for (int m=0; m<numbmodal; m++) {
            //Data_PTRModalVector.push_back(&static_cast<float*>(this->GetDataImage()->data)[m*numel]);
            Data_PTRModalVector.push_back(&Data[m*numelmasked]);
        }
        float * NormResp_PTR=this->GetNormResp();
        float tmpValue[MaxNumbModal];//{0};

        for(int i=0;i<MaxNumbModal;i++){
            tmpValue[i]=0;
        }
        int twopowmodal=(int)pow_int(2,numbmodal);
        int twopowMaxmodal=(int)pow_int(2,MaxNumbModal);
//        int indexTab[twopowMaxmodal];//{0};
        int *indexTab=new int[twopowMaxmodal];
        for(int i=0;i<twopowMaxmodal;i++){
            indexTab[i]=0;
        }
//        float valuePerc[twopowMaxmodal];//{1};
        float *valuePerc=new float[twopowMaxmodal];
        for(int i=0;i<twopowMaxmodal;i++){
            valuePerc[i]=1;
        }
        int CountNormRespZero=0;
        float SumNormResp=0;
        for (int i=0; i<numelmasked; i++) {
            for (int m=0; m<numbmodal; m++) {
                int twopowm=(int)pow_int(2,m);
                tmpValue[m]=*Data_PTRModalVector[m]/sizeBin-0.5;
                tmpValue[m]=tmpValue[m]>(numbbins-1)?numbbins-1:tmpValue[m];
                tmpValue[m]=tmpValue[m]<0?0:tmpValue[m];
                for (int c=0; c<twopowm; c++) {
                    indexTab[c+twopowm]=indexTab[c];
                    valuePerc[c+twopowm]=valuePerc[c];
                }
                int binspowm=(int)pow_int(numbbins,m);
                for (int c=0; c<twopowm; c++) {
                    indexTab[c]+=binspowm*floorf(tmpValue[m]);
                    indexTab[c+twopowm]+=binspowm*floorf(tmpValue[m])+1;
                    valuePerc[c]*=(1-tmpValue[m]+floorf(tmpValue[m]));
                    valuePerc[c+twopowm]*=(tmpValue[m]-floorf(tmpValue[m]));
                }
            }
            float sumValuePerc=0;
            for (int c=0; c<twopowmodal; c++) {
                sumValuePerc+=valuePerc[c];
            }
            //cout<<"sumValuePerc "<< sumValuePerc<<" for index "<<i<<endl;
            if(NormResp_PTR!=NULL){

                for (int p=0; p<twopowmodal; p++) {
                    DataHistogramTotal_PTR[indexTab[p]]+=(float)valuePerc[p]*(*NormResp_PTR);

                }
                if (*NormResp_PTR<1-1E-6) {
                    CountNormRespZero++;
                }
                SumNormResp+=(PrecisionTYPE)*NormResp_PTR;
                NormResp_PTR++;
            }
            else{
                for(int p=0; p<twopowmodal; p++) {
                    DataHistogramTotal_PTR[indexTab[p]]+=(float)valuePerc[p];

                }
                SumNormResp++;
            }

            for (int p=0; p<twopowmodal; p++) {
                indexTab[p]=0;
                valuePerc[p]=1;
            }
            for (int m=0; m<numbmodal; m++) {
                Data_PTRModalVector[m]++;
            }
        }
        //cout<<"CountNormRespZero is"<<CountNormRespZero<<endl;
                delete [] indexTab;
                indexTab=NULL;
                delete [] valuePerc;
                valuePerc=NULL;
//                delete [] tmpValue;
//                tmpValue=NULL;
        float SumDataHistogram=0;

        //        float NumberDataHistogramZero=0;
        //        int numelmasked=this->GetNumberMaskedElements();
        int binspowmodal=(int)pow_int(numbbins,numbmodal);

        //Normalisation of the DataHistogram in order to properly calculate the KLD afterwards:
        for(int i=0;i<binspowmodal;i++){
            SumDataHistogram+=DataHistogramTotal_PTR[i];
        }
        if(SumDataHistogram>0){
        for(int i=0;i<binspowmodal;i++){
            DataHistogramTotal_PTR[i]/=SumDataHistogram;
        }
        }

//        for (int i=0; i<binspowmodal; i++) {
//            DataHistogramTotal_PTR[i]/=SumNormResp;
//            SumDataHistogram+=DataHistogramTotal_PTR[i];
//        }

        //        cout<<"Sum DataHistogram is "<<SumDataHistogram<<endl;
        //        if (Data!=NULL) {
        //            delete [] Data;
        //            Data=NULL;
        //        }
        return DataHistogramTotal_PTR;
    }
}

// Returns the histogram of the distribution for numbbins
vector<float *> TreeEM::GetDistHistogram(){
    vector<float *> DistHistogramVector;
    // First check if there is Data to make histogram
    if (this->GetDataImage()==NULL) {
        cout<<"No data to make histogram from"<<endl;
        return DistHistogramVector;
    }
    // Then if modality wanted in histogram is available
    //vector<float *> DistHistogramVector;
    int numbchild=this->GetNumberChildren();
    int numbmodal=this->GetNumberModalities();
    int numbbins=this->GetNumbbins();
    if (this->GetNumberChildren()==0) {/* If it is a leaf the distribution is then calculated according to the parameters stored. For the moment, their is only one choice of simple distribution : the Gaussian distribution*/
        switch (this->GetDistributionType()) {
        case 2:{ // case of a uniform distribution
            vector<float*>ChildHistogram_VEC=this->GetUniformDistributionHist();
            for (int m=0; m<numbmodal; m++) {
                DistHistogramVector.push_back(ChildHistogram_VEC[m]);
            }
            break;
        }
        default:
            vector<float*>ChildHistogram_VEC=this->GetGaussianDistributionHist();
            for (int m=0; m<numbmodal; m++) {
                DistHistogramVector.push_back(ChildHistogram_VEC[m]);
            }
            break;
        }
    }
    else{/* It is not a leaf and the distribution is the weighted sum of the distributions below*/
        float * ChildHistogram_PTR;
        float * DistHistogram_PTR;

        float NormWeightChild=1;
        for (int m=0; m<numbmodal; m++) {
            float * InitDistHistToPushBack=new float[numbbins];
            for(int i=0;i<numbbins;i++){
                InitDistHistToPushBack[i]=0;
            }
            DistHistogramVector.push_back(InitDistHistToPushBack);
        }
        for (int c=0; c<numbchild; c++) {
            vector<float*>ChildHistogram_VEC=this->GetChild(c)->GetDistHistogram();
            for (int m=0; m<numbmodal; m++) {
                DistHistogram_PTR=DistHistogramVector[m];
                ChildHistogram_PTR=ChildHistogram_VEC[m];
                NormWeightChild=this->GetChild(c)->GetNormWeight();
                for (int i=0; i<numbbins; i++,DistHistogram_PTR++,ChildHistogram_PTR++) {
                    *DistHistogram_PTR+=NormWeightChild*(*ChildHistogram_PTR);
                }
                // ClearingMemory
                if(ChildHistogram_VEC[m]!=NULL){
                    delete ChildHistogram_VEC[m];
                    ChildHistogram_VEC[m]=NULL;
//                    cout<< " Delete done "<<endl;
                }
                //                if (ChildHistogram_VEC[m]!=NULL && ((*ChildHistogram_VEC[m])!= NULL) ) {
                //                                    delete [] ChildHistogram_VEC[m];
                //                                    cout<< "Delete ChildHistogram_VEC partially done for modality "<<m<<" and child "<<c<<endl;
                //                                   ChildHistogram_VEC[m]=NULL;
                //                }
                //                if (ChildHistogram_VEC[m]==NULL) {
                //                    cout<<"ChildHistogram_VEC pointer already NULL"<<endl;
                //                }
                //                if ((*ChildHistogram_VEC[m])==NULL) {
                //                    cout<<"Pb with memory allocation"<<endl;
                //                }

            }
        }
    }
    return DistHistogramVector;
}

// Returns the vector containing the uniform distribution for each modality in a vector of float pointers
vector<float*> TreeEM::GetUniformDistributionHist(){
    vector<float *> UniformDistributionHist;
    if (this->GetNumberChildren()!=0) {
        cout<< "not leaf so cannot return gaussian distribution"<<endl;
        return UniformDistributionHist;
    }
    if(this->GetDistributionType()!=2){
        cout<<"no uniform distribution type"<<endl;
        return UniformDistributionHist;
    }
    //vector<float *> GaussianDistributionHist;
    int numbmodal=this->GetNumberModalities();
    float sizeBin=1.0/(float)numbbins;
    for(int m=0;m<numbmodal;m++){
        float * UniformDistTmp=new float[numbbins];
        for(int i=0;i<numbbins;i++){
            UniformDistTmp[i]=sizeBin;
        }
        UniformDistributionHist.push_back(UniformDistTmp);
    }
    return UniformDistributionHist;
}

// Returns the Gaussian Distribution for each modality independently in a vector of float pointers
vector<float *> TreeEM::GetGaussianDistributionHist(){
    vector<float *> GaussianDistributionHist;
    if (this->GetNumberChildren()!=0) {
        cout<< "not leaf so cannot return gaussian distribution"<<endl;
        return GaussianDistributionHist;
    }
    if(this->GetDistributionType()!=1){
        cout<<"no gaussian distribution type"<<endl;
        return GaussianDistributionHist;
    }
    //vector<float *> GaussianDistributionHist;
    int numbmodal=this->GetNumberModalities();

    float * VarianceB=this->GetVariance();
    float Variance[MaxNumbModal*MaxNumbModal];
    for(int m1=0;m1<MaxNumbModal;m1++){
        for(int m2=0;m2<MaxNumbModal;m2++){
            if(m1<numbmodal&&m2<numbmodal){
                Variance[m1+m2*MaxNumbModal]=VarianceB[m1+m2*numbmodal];
            }
            else{
                Variance[m1+m2*MaxNumbModal]=0;
            }
        }
    }
    float * MeanB=this->GetMean();
    float Mean[MaxNumbModal];
    for(int m=0;m<MaxNumbModal;m++){
        if(m<numbmodal){
            Mean[m]=MeanB[m];
        }
        else{
            Mean[m]=0;
        }
    }
    float Value=0;
    float sizeBin=1.0/(float)numbbins;
    for (int m=0; m<numbmodal; m++) {
        float * GaussianDist_tmp=new float[numbbins];
        for(int i=0;i<numbbins;i++){
            GaussianDist_tmp[i]=0;
        }
        float * GaussianDist_tmpPTR=GaussianDist_tmp;
        float NormalisationFactor=1.0/powf(2*M_PI*Variance[m+m*MaxNumbModal], 0.5);
        GaussianDist_tmpPTR=GaussianDist_tmp;
        float sumDist=0;
        for (int i=0; i<numbbins; i++,GaussianDist_tmpPTR++) {
            Value=i*sizeBin+0.5*sizeBin;
            *GaussianDist_tmpPTR=NormalisationFactor*expf(-(float)0.5/Variance[m+m*MaxNumbModal]*(Value-Mean[m])*(Value-Mean[m]))*sizeBin;
            sumDist+=*GaussianDist_tmpPTR;
        }
        //        cout<<"sumDist in vectorial Gaussian hist is "<<sumDist<<endl;
        GaussianDistributionHist.push_back(GaussianDist_tmp);
    }
    //    delete [] GaussianDist_tmp;
    //    GaussianDist_tmp=NULL;
    return GaussianDistributionHist;
}

// Returns the complete histogram corresponding to the uniform distribution
float * TreeEM::GetUniformDistributionHistTotal(){
    if (this->GetNumberChildren()!=0) {
        cout<< "not leaf so cannot return uniform distribution"<<endl;
        return NULL;
    }
    if(this->GetDistributionType()!=2){
        cout<<"no uniform distribution type"<<endl;
        return NULL;
    }

    int numbmodal=this->GetNumberModalities();
    int binspowmodal=(int)pow_int(numbbins,numbmodal);
    float * UniformDistributionTotal = new float[binspowmodal];
    float sizeBin=1.0/(float)numbbins;
    float sizeBinpowmodal=pow_int(sizeBin,numbmodal);
    for(int i=0;i<binspowmodal;i++){
        UniformDistributionTotal[i]=sizeBinpowmodal;
    }
    return UniformDistributionTotal;
}

// Returns the complete histogram corresponding to the distribution with the given number of bins. To be used for the derivation of the KLD
float * TreeEM::GetGaussianDistributionHistTotal(){
    if (this->GetNumberChildren()!=0) {
        cout<< "not leaf so cannot return gaussian distribution"<<endl;
        return NULL;
    }
    if(this->GetDistributionType()!=1){
        cout<<"no gaussian distribution type"<<endl;
        return NULL;
    }

    int numbmodal=this->GetNumberModalities();

    // First initialise the Variance Matrix and inverted one
//    matrix<float> VarianceMatrix=matrix<float>(numbmodal);
    float * VarianceToSet=this->GetVariance();

//    for (int m1=0; m1<numbmodal; m1++) {
//        for (int m2=0; m2<numbmodal; m2++) {
//            VarianceMatrix.setvalue(m1, m2, VarianceToSet[m1+m2*numbmodal]);
//        }
//    }

    float * VarianceMatrixToInvert=new float[numbmodal*numbmodal];
    for (int m1=0; m1<numbmodal; m1++) {
        for (int m2=0; m2<numbmodal; m2++) {
            VarianceMatrixToInvert[m1+m2*numbmodal]=VarianceToSet[m1+m2*numbmodal];
        }
    }

    //    matrix<float>VarianceMatrixCopy=matrix<float>(MaxNumbModal);
    //    VarianceMatrixCopy.copymatrix(VarianceMatrix);


    // Calculation of the factor in front of the exponential
//    float DeterminantVariance=VarianceMatrix.determinant();
    float DeterminantVariance=determinant(VarianceMatrixToInvert,numbmodal);
    float * MeanB=this->GetMean();
    float Mean[MaxNumbModal];
    for(int m=0;m<MaxNumbModal;m++){
        if(m<numbmodal){
            Mean[m]=MeanB[m];
        }
        else{
            Mean[m]=0;
        }
    }
    float sizeBin=1.0/(float)numbbins;
    //cout<<"the Variance determinant is "<<DeterminantVariance<<endl;
    float NormalisationFactor=1.0/(float)(powf(2*M_PI , (float)((float)numbmodal/2.0))*powf(DeterminantVariance, 0.5));
    //cout <<"The normalisation factor is "<<NormalisationFactor<<endl;
    //Initialisation of the needed element to calculate the inside of the exponential
    int binspowmodal=(int)pow_int(numbbins,numbmodal);
    float * GaussianDistributionTotal = new float[binspowmodal];//{0};
    for(int i=0;i<binspowmodal;i++){
        GaussianDistributionTotal[i]=0;
    }
    float * GaussianDistribution_PTR=GaussianDistributionTotal;
    float InvertedVariance[MaxNumbModal*MaxNumbModal];
    for(int i=0;i<MaxNumbModal*MaxNumbModal;i++){
        InvertedVariance[i]=0;
    }
    if (numbmodal>1) {
//        VarianceMatrix.invert();
        invertMatrix(VarianceMatrixToInvert,numbmodal);
        //        matrix<float>TestMatrix=matrix<float>(MaxNumbModal);
        //        TestMatrix.settoproduct(VarianceMatrixCopy,VarianceMatrix);
        //        TestMatrix.comparetoidentity();
//        bool success;
        for(int m1=0;m1<numbmodal;m1++){
            for (int m2=0; m2<numbmodal; m2++) {
//                VarianceMatrix.getvalue(m1, m2, InvertedVariance[m1+m2*MaxNumbModal], success);
                InvertedVariance[m1+m2*MaxNumbModal]=VarianceMatrixToInvert[m1+m2*numbmodal];
            }
        }
        delete [] VarianceMatrixToInvert;
        VarianceMatrixToInvert=NULL;
    }
    else{
//        bool success;
//        VarianceMatrix.getvalue(0, 0, InvertedVariance[0], success);
        *InvertedVariance=(float)1.0/DeterminantVariance;
        //*InvertedVariance=1/(*InvertedVariance);
    }
    float temp;
    int  IndexConversion [MaxNumbModal];//{0};
    for(int i=0;i<MaxNumbModal;i++){
        IndexConversion[i]=0;
    }
    int Rem;
    // Calculation of the inside of the exponential
    //cout<<InvertedVariance[0]<<endl;
    int CountNegativeValues=0;
    for (int i=0; i<binspowmodal; i++) {
        Rem=i;

        // Conversion of the index i into the different index for each modality giving then the value to consider
        for (int m=numbmodal-1; m>=0; m--) {
            int binspowm=(int)pow_int(numbbins,m);
            IndexConversion[m]=Rem/binspowm;
            Rem-=IndexConversion[m]*binspowm;
        }
        for (int m1=0; m1<numbmodal; m1++) {
            temp=0;
            for (int m2=0; m2<numbmodal; m2++) {
                temp+=((sizeBin*IndexConversion[m2]+0.5*sizeBin)-Mean[m2])*InvertedVariance[m2+m1*MaxNumbModal];
            }
            *GaussianDistribution_PTR=(*GaussianDistribution_PTR)+temp*((sizeBin*IndexConversion[m1]+0.5*sizeBin)-Mean[m1]);
        }
        if (*GaussianDistribution_PTR<0) {
            CountNegativeValues++;
        }
        GaussianDistribution_PTR++;
    }
    //cout<<"The number of negative values is "<<CountNegativeValues<<endl;
    // Filling the distribution array with the value
    GaussianDistribution_PTR=GaussianDistributionTotal;
    //cout << "Normalisation factor in GetGaussianDistHistTotal is "<<NormalisationFactor<<endl;
    float sumDist=0;
    //    int numelmasked=this->GetNumberMaskedElements();
    float sizeBinpowmodal=pow_int(sizeBin,numbmodal);
    for (int i=0;i<binspowmodal;i++,GaussianDistribution_PTR++){
        *GaussianDistribution_PTR=NormalisationFactor*expf(-0.5*(*GaussianDistribution_PTR))*sizeBinpowmodal;
        sumDist+=*GaussianDistribution_PTR;
    }
    //cout<<"sumDist in GetGaussianTotal is "<<sumDist<<endl;
    // Clearing space needing for inverse of variance


    // Return result
    return GaussianDistributionTotal;
}

float * TreeEM::GetDistHistogramTotal(){

    // Then if modality wanted in histogram is available
    int numbchild=this->GetNumberChildren();
    int numbmodal=this->GetNumberModalities();
    int numbbins=this->GetNumbbins();
    int binspowmodal=(int)pow_int(numbbins,numbmodal);
    float * DistHistogramTotal=new float[binspowmodal];//{0};
    for(int i=0;i<binspowmodal;i++){
        DistHistogramTotal[i]=0;
    }

    if (this->GetNumberChildren()==0) {/* If it is a leaf the distribution is then calculated according to the parameters stored. For the moment, their is only one choice of simple distribution : the Gaussian distribution*/
        switch (this->GetDistributionType()) {
        case 2:{ // case of uniform distribution
            float * LeafHistogramTotal=this->GetUniformDistributionHistTotal();
            float * LeafHistogramTotal_PTR=LeafHistogramTotal;
            for (int i=0; i<binspowmodal; i++,LeafHistogramTotal_PTR++) {
                DistHistogramTotal[i]=*LeafHistogramTotal_PTR;
                //SumLeaf+=*LeafHistogramTotal_PTR;
            }
            //cout << "Sum leaf is "<<SumLeaf<<endl;
            delete [] LeafHistogramTotal;
            LeafHistogramTotal=NULL;
            LeafHistogramTotal_PTR=NULL;
            break;
        }
        default:
            float * LeafHistogramTotal=this->GetGaussianDistributionHistTotal();
            //            vector<float *> LeafVectorHistogram=this->GetGaussianDistributionHist();
            //            float * LeafVectorHistogram_PTR=LeafVectorHistogram[0];
            float * LeafHistogramTotal_PTR=LeafHistogramTotal;
            //            int CountDifference=0;
            //            for (int i=0; i<numbbins; i++,LeafHistogramTotal_PTR++,LeafVectorHistogram_PTR++) {
            //                if (fabsf((*LeafHistogramTotal_PTR)-(*LeafVectorHistogram_PTR))>1E-6) {
            //                    CountDifference++;
            //                }
            //            }
            //            //cout<<"Count difference direct gaussian is "<<CountDifference<<endl;
            //            LeafHistogramTotal_PTR=LeafHistogramTotal;
            //            float SumLeaf=0;
            for (int i=0; i<binspowmodal; i++,LeafHistogramTotal_PTR++) {
                DistHistogramTotal[i]=*LeafHistogramTotal_PTR;
                //SumLeaf+=*LeafHistogramTotal_PTR;
            }
            //cout << "Sum leaf is "<<SumLeaf<<endl;
            delete [] LeafHistogramTotal;
            LeafHistogramTotal=NULL;
            LeafHistogramTotal_PTR=NULL;
            break;
        }
    }
    else{/* It is not a leaf and the distribution is the weighted sum of the distributions below*/

        float * DistHistogram_PTR=DistHistogramTotal;
        float * ChildDistHistogram=NULL;
        float * ChildDistHistogram_PTR=NULL;
        float NormWeightChild;
        for (int c=0; c<numbchild; c++) {

            ChildDistHistogram=(this->GetChild(c))->GetDistHistogramTotal();
            ChildDistHistogram_PTR=ChildDistHistogram;
            DistHistogram_PTR=DistHistogramTotal;
            NormWeightChild=this->GetChild(c)->GetNormWeight();
            for (int i=0; i<binspowmodal; i++) {
                *DistHistogram_PTR+=NormWeightChild*(*ChildDistHistogram_PTR);
                DistHistogram_PTR++;
                ChildDistHistogram_PTR++;
            }
            // Clearing
            ChildDistHistogram_PTR=NULL;
            if (ChildDistHistogram!=NULL) {
                delete [] ChildDistHistogram;
                ChildDistHistogram=NULL;
//                cout<< "Delete done for child "<<c<<endl;
            }
            //        if ((ChildDistHistogram!=NULL) && ((*ChildDistHistogram)!= NULL)) {
            //
            //            delete [] ChildDistHistogram;
            //            ChildDistHistogram=NULL;
            //            cout<<"Delete done for child "<<c<<endl;
            //        }
            //        if (ChildDistHistogram==NULL) {
            //            cout<<"ChildDistHistogram is NULL"<<endl;
            //        }
            //        cout<<ChildDistHistogram<<" and "<< *ChildDistHistogram<<endl;


        }
    }
    return DistHistogramTotal;
}

// Returns the KLD when comparing the total histogram of the distribution modelled and the data
float TreeEM::GetDistCompKLDTotal(){
    float * DistHistTotal=this->GetDistHistogramTotal();
    float * DistHistTotal_PTR=DistHistTotal;
    float * DataHistTotal=this->GetDistHistogramTotal();
    float * DataHistTotal_PTR=DataHistTotal;
    PrecisionTYPE KLD=0;
    int numbmodal=this->GetNumberModalities();
    int binspowmodal=(int)pow_int(numbbins,numbmodal);
    for (int i=0; i<binspowmodal; i++,DataHistTotal_PTR++,DistHistTotal_PTR++) {
        KLD+=(PrecisionTYPE)logf(((*DataHistTotal_PTR)/(*DistHistTotal_PTR)))*(*DataHistTotal_PTR);
    }
    delete [] DistHistTotal;
    DistHistTotal=NULL;
    delete [] DataHistTotal;
    DataHistTotal = NULL;
    return (float)KLD;
}

float * TreeEM::GetDistCompKLD(){
    int numbmodal=this->GetNumberModalities();
    int binspowmodal=(int)pow_int(numbbins,numbmodal);
    float * KLD=new float[numbmodal];//{0};
    for(int i=0;i<numbmodal;i++){
        KLD[i]=0;
    }
    vector<float *> DistHist_VEC=this->GetDistHistogram();
    vector<float *> DataHist_VEC=this->GetDataHistogram();

    for (int m=0; m<numbmodal; m++) {
        PrecisionTYPE KLD_tmp=0;
        float * DistHist_PTR=DistHist_VEC[m];
        float * DataHist_PTR=DataHist_VEC[m];
        for (int i=0; i<binspowmodal; i++, DataHist_PTR++,DistHist_PTR++) {
            KLD_tmp+=(PrecisionTYPE)logf(((*DataHist_PTR)/(*DistHist_PTR)))*(*DataHist_PTR);
        }
        KLD[m]=KLD_tmp;
        delete [] DistHist_VEC[m];
        DistHist_VEC[m]=NULL;
        delete [] DataHist_VEC[m];
        DataHist_VEC[m]=NULL;
    }
    return KLD;
}

float TreeEM::CriterionCalculation(){
    int numbFreeParams=this->GetNumberFreeParameters();
    numbFreeParams+=this->GetNumberAllLeaves()-1;
    int numelmasked=this->GetNumberMaskedElements();
    float IndFactor=this->GetIndFactor();
    float LL=this->GetLogLikelihood();
    float BIC=LL*IndFactor-(float)numbFreeParams/2.0*logf(numelmasked*IndFactor);
    return BIC;
}

float TreeEM::CriterionCalculationDPInd(){

    float LL=this->GetLogLikelihood();
    int numbchild=this->GetNumberChildren();
    int numelmasked=this->GetNumberMaskedElements();
    float DP=1;
    float * DPChildren=this->GetDPChildren();
//    if(this->IsRoot() && DPChildren.size()==numbchild){
    if(this->IsRoot()){
    for(int c=0;c<numbchild;c++ ){
            //float * AlphaResult=AlphaCalculation(DPtemp);
            int numbLeaves=this->GetChild(c)->GetNumberAllLeaves();
            if(numbLeaves >=2){
                DP*=DPChildren[c*MaxSupport+numbLeaves-1];
            }
            else{
                DP*=DPChildren[0+c*MaxSupport];
            }

        }
    }
    int numbFreeParams=this->GetNumberFreeParameters();
    float IndFactor=this->GetIndFactor();
    return (LL+numelmasked*logf(DP))*IndFactor-(float)numbFreeParams/2.0*logf(numelmasked*IndFactor);
}

float TreeEM::CriterionCalculationDPNonInd(int ChildModified){

    float LL=this->GetLogLikelihood();
    int numbchild=this->GetNumberChildren();
    int numelmasked=this->GetNumberMaskedElements();
    float DP=1;
    float logDP=0;
    float * DPChildren=this->GetDPChildren();
//    if(this->IsRoot() && DPChildren.size()==numbchild){
    if(this->IsRoot()){
        int * numbLeaves=new int[numbchild];
        int indexHistogram=0;
        for(int c=0;c<numbchild;c++){
            numbLeaves[c]=this->GetChild(c)->GetNumberAllLeaves();
            indexHistogram+=numbLeaves[c]*pow_int(MaxSupport,c);

        }
        int Shift=pow_int(MaxSupport,ChildModified);
        int indexInitHist=indexHistogram-numbLeaves[ChildModified]*Shift;
        float Denominator=0;
        for (int i=0;i<MaxSupport;i++){
            Denominator+=DPChildren[indexInitHist+i*Shift];
        }
        float Numerator=DPChildren[indexHistogram];
         DP=Numerator/Denominator;
         delete [] numbLeaves;
         DP=Numerator;
         numbLeaves=NULL;
         float * LogDP=LogGaussBlur(DPChildren,pow_int(MaxSupport,numbchild));
//         cout<< LogDP[indexHistogram];
         logDP=LogDP[indexHistogram];

    }
    float IndFactor=this->GetIndFactor();
    cout<<"LL component : "<< LL*IndFactor << "DP Coeff "<<logDP<<endl;
    int numbFreeParams=this->GetNumberFreeParameters();

//    return (LL*IndFactor)+logf(DP)-(float)numbFreeParams/2.0*logf(numelmasked*IndFactor);
    return (LL*IndFactor)+logDP-(float)numbFreeParams/2.0*logf(numelmasked*IndFactor);
}

float TreeEM::CriterionCalculationSplitDP(){

    float LL=this->GetLogLikelihood();
    int numbchild=this->GetNumberChildren();
    int numelmasked=this->GetNumberMaskedElements();
    float DP=1;
    float * DPChildren=this->GetDPChildren();
    if(this->IsRoot() && DPChildren!=NULL){
        for(int c=0;c<numbchild;c++ ){
            float * DPtemp=DPChildren;
            float DPtemp_test[MaxSupport];
            for(int i=0;i<MaxSupport;i++){
                DPtemp_test[i]=0;
            }
            for(int i=0;i<MaxSupport;i++){
                float tmp=0;
                for(int j=MaxSupport-1;j>=i;j--){
                    tmp+=DPtemp[j];
                }
                DPtemp_test[i]=tmp;
            }
            //float * AlphaResult=AlphaCalculation(DPtemp);
            int numbLeaves=this->GetChild(c)->GetNumberAllLeaves();
            if(numbLeaves >=2){
                DP*=DPtemp_test[numbLeaves-1];
            }
            else{
                DP*=1;
            }

        }
    }
    return LL+numelmasked*logf(DP);
}

float TreeEM::CriterionCalculationMergeDP(){

    float LL=this->GetLogLikelihood();
    int numbchild=this->GetNumberChildren();
    int numelmasked=this->GetNumberMaskedElements();
    float DP=1;
    float * DPChildren=this->GetDPChildren();
    if(this->IsRoot() && DPChildren!=NULL){
        for(int c=0;c<numbchild;c++ ){
            float * DPtemp=DPChildren;
            float DPtemp_test[MaxSupport];
            for(int i=0;i<MaxSupport;i++){
                DPtemp_test[i]=0;
            }
            for(int i=0;i<MaxSupport;i++){
                float tmp=0;
                for(int j=0;j<=i;j++){
                    tmp+=DPtemp[j];
                }
                DPtemp_test[i]=tmp;
            }
            //float * AlphaResult=AlphaCalculation(DPtemp);
            int numbLeaves=this->GetChild(c)->GetNumberAllLeaves();
            if(numbLeaves >=2){
                DP*=DPtemp_test[numbLeaves-1];
            }
            else{
                DP*=DPtemp_test[0];
            }

        }
    }
    return LL+numelmasked*logf(DP);
}

/************************** METHODS FOR MERGING OPERATIONS ******************************/

// Returns the comparison between the two modelled distributions if comparable meaning that the two trees compared have the same Parent
MergeKLD * TreeEM::GetKLDMerge(int Child1Input,int Child2Input){
    // Check that the given Children exist
    int numbchild=this->GetNumberChildren();
    if (Child1Input>=numbchild||Child2Input>=numbchild||Child1Input<0||Child2Input<0) {
//        cout<<"Out of bound for the children to merge"<<endl;
        return NULL;
    }
    // Check if the merging try is made on leaves
    if (this->GetChild(Child1Input)->GetNumberChildren()!=0||this->GetChild(Child2Input)->GetNumberChildren()!=0) {
//        cout<< "One of the children is not a leaf, merging cannot be done"<<endl;
        return NULL;
    }
    // Check if the merging is tried on leaves with same priors
    //    cout<< this->GetChild(Child1Input)->GetPriors()<<" and "<<this->GetChild(Child2Input)->GetPriors();
    if (this->GetChild(Child1Input)->GetPriors()!=this->GetChild(Child2Input)->GetPriors()) {
//        cout<<"Impossible to merge leaves with not same priors"<<endl;
        return NULL;
    }
    // Check if merging is done on leaves with same distribution type
    if(this->GetChild(Child1Input)->GetDistributionType()!=this->GetChild(Child2Input)->GetDistributionType()){
        cout<<"No merging for components without same distribution";
        return NULL;
    }
    if(this->IsRoot()){
//        cout<<"Impossible to merge general classes"<<endl;
        return NULL;
    }
    int numbmodal=this->GetNumberModalities();
    int binspowmodal=(int)pow_int(numbbins,numbmodal);
    // Initialisation of ResultMergeKLD
    MergeKLD * ResultMergeKLD=new MergeKLD();
    ResultMergeKLD->Parent=this;
    ResultMergeKLD->Child1=Child1Input;
    ResultMergeKLD->Child2=Child2Input;
    ResultMergeKLD->KLDTot=0;
//    ResultMergeKLD->KLD=new float[MaxNumbModal];//{0};
//    for(int i=0;i<MaxNumbModal;i++){
//        ResultMergeKLD->KLD[i]=0;
//    }

    // Calculation of the needed histograms for the calculation of the KLDcompTot
    float * DistHistTotal1=this->GetChild(Child1Input)->GetDistHistogramTotal();
    float * DistHistTotal1_PTR=DistHistTotal1;
    float * DistHistTotal2=this->GetChild(Child2Input)->GetDistHistogramTotal();
    float * DistHistTotal2_PTR=DistHistTotal2;
    //    int CountNan1=0;
    //    int CountNan2=0;
    //    int CountNanAdd=0;

//    // check sum of distribution in order to calculate KLD
//    PrecisionTYPE sumDist1=0;
//    PrecisionTYPE sumDist2=0;
//    for(int i=0;i<binspowmodal;i++){
//        sumDist1+=(PrecisionTYPE)DistHistTotal1[i];
//        sumDist2+=(PrecisionTYPE)DistHistTotal2[i];
//    }
//    cout<< "sumDist1 is "<<sumDist1<< " and sumDist2 is "<<sumDist2<<endl;


    for (int i=0; i<binspowmodal; i++,DistHistTotal1_PTR++,DistHistTotal2_PTR++) {
        float DistValue1=(*DistHistTotal1_PTR)<=0.00001?0.00001:(*DistHistTotal1_PTR);
        float DistValue2=(*DistHistTotal2_PTR)<=0.00001?0.00001:(*DistHistTotal2_PTR);
        //        if (logf(DistValue1/DistValue2)*(*DistHistTotal1_PTR)!=logf(DistValue1/DistValue2)*(*DistHistTotal1_PTR)) {
        //            CountNan1++;
        //        }
        //        if (logf(DistValue2/DistValue1)*(*DistHistTotal2_PTR)!=logf(DistValue2/DistValue1)*(*DistHistTotal2_PTR)) {
        //            CountNan2++;
        //        }
        //        if (0.5*(logf(DistValue1/DistValue2)*(*DistHistTotal1_PTR)+logf(DistValue2/DistValue1)*(*DistHistTotal2_PTR))!=0.5*(logf(DistValue1/DistValue2)*(*DistHistTotal1_PTR)+logf(DistValue2/DistValue1)*(*DistHistTotal2_PTR))) {
        //            CountNanAdd++;
        //        }
        ResultMergeKLD->KLDTot+=0.5*(logf(((DistValue1)/(DistValue2)))*(*DistHistTotal1_PTR)+logf(((DistValue2)/(DistValue1)))*(*DistHistTotal2_PTR)); // Use of the symmetric KLD between the two distributions
        //cout<<ResultMergeKLD->KLDTot<< "   -   ";
    }
//    cout<<endl;
    delete [] DistHistTotal1;
    DistHistTotal1=NULL;
    delete [] DistHistTotal2;
    DistHistTotal2=NULL;

//    // Calculation of KLDcomp
//    vector<float *> DistHist1_VEC=this->GetChild(Child1Input)->GetDistHistogram();
//    vector<float *> DistHist2_VEC=this->GetChild(Child2Input)->GetDistHistogram();
//    for (int m=0; m<numbmodal; m++) {
//        float * DistHist1_PTR=DistHist1_VEC[m];
//        float * DistHist2_PTR=DistHist2_VEC[m];
//        for (int i=0; i<numbbins; i++,DistHist1_PTR++,DistHist2_PTR++) {
//            float DistValue1=(*DistHist1_PTR)<=0.00001?0.00001:(*DistHist1_PTR);
//            float DistValue2=(*DistHist2_PTR)<=0.00001?0.00001:(*DistHist2_PTR);
//            ResultMergeKLD->KLD[m]+=0.5*(logf(((DistValue1/DistValue2)))*(*DistHist1_PTR)+logf(((DistValue2)/(DistValue1)))*(*DistHist2_PTR)); // Use of the symmetric KLD between the two distributions
//        }
//        delete [] DistHist1_VEC[m];
//        DistHist1_VEC[m]=NULL;
//        delete [] DistHist2_VEC[m];
//        DistHist2_VEC[m]=NULL;
//    }
    return ResultMergeKLD;
}

GenKLD * TreeEM::GetGenKLD(){
    int numbmodal=this->GetNumberModalities();
    int binspowmodal=(int)pow_int(numbbins,numbmodal);
    // Initialisation of ResultMergeKLD
    KLD * ResultKLD=new KLD();

    ResultKLD->KLDTot=0;
//    ResultKLD->KLD=new float[MaxNumbModal];//{0};
//    for(int i=0;i<MaxNumbModal;i++){
//        ResultKLD->KLD[i]=0;
//    }

    // Calculation of the needed histograms for the calculation of the KLDcompTot
    float * DistHistTotal1=this->GetDistHistogramTotal();
    float * DistHistTotal1_PTR=DistHistTotal1;
    float * DistHistTotal2=this->GetDataHistogramTotal();
    float * DistHistTotal2_PTR=DistHistTotal2;
    //    int CountNan1=0;
    //    int CountNan2=0;
    //    int CountNanAdd=0;
    PrecisionTYPE tmpResultKLD=0;
    for (int i=0; i<binspowmodal; i++,DistHistTotal1_PTR++,DistHistTotal2_PTR++) {
        float DistValue1=(*DistHistTotal1_PTR)<=0.00001?0.00001:(*DistHistTotal1_PTR);
        float DistValue2=(*DistHistTotal2_PTR)<=0.00001?0.00001:(*DistHistTotal2_PTR);
        //        if (logf(DistValue1/DistValue2)*(*DistHistTotal1_PTR)!=logf(DistValue1/DistValue2)*(*DistHistTotal1_PTR)) {
        //            CountNan1++;
        //        }
        //        if (logf(DistValue2/DistValue1)*(*DistHistTotal2_PTR)!=logf(DistValue2/DistValue1)*(*DistHistTotal2_PTR)) {
        //            CountNan2++;
        //        }
        //        if (0.5*(logf(DistValue1/DistValue2)*(*DistHistTotal1_PTR)+logf(DistValue2/DistValue1)*(*DistHistTotal2_PTR))!=0.5*(logf(DistValue1/DistValue2)*(*DistHistTotal1_PTR)+logf(DistValue2/DistValue1)*(*DistHistTotal2_PTR))) {
        //            CountNanAdd++;
        //        }
        tmpResultKLD+=(PrecisionTYPE)0.5*(log(((DistValue1)/(DistValue2)))*(*DistHistTotal1_PTR)+log(((DistValue2)/(DistValue1)))*(*DistHistTotal2_PTR)); // Use of the symmetric KLD between the two distributions
        //cout<<ResultMergeKLD->KLDTot<< "   -   ";
    }
    ResultKLD->KLDTot=(float)tmpResultKLD;
//    cout<<endl;
    delete [] DistHistTotal1;
    DistHistTotal1=NULL;
    delete [] DistHistTotal2;
    DistHistTotal2=NULL;

//    // Calculation of KLDcomp
//    vector<float *> DistHist1_VEC=this->GetDistHistogram();
//    vector<float *> DistHist2_VEC=this->GetDataHistogram();
//    for (int m=0; m<numbmodal; m++) {
//        float * DistHist1_PTR=DistHist1_VEC[m];
//        float * DistHist2_PTR=DistHist2_VEC[m];
//        for (int i=0; i<numbbins; i++,DistHist1_PTR++,DistHist2_PTR++) {
//            float DistValue1=(*DistHist1_PTR)<=0.00001?0.00001:(*DistHist1_PTR);
//            float DistValue2=(*DistHist2_PTR)<=0.00001?0.00001:(*DistHist2_PTR);
//            ResultKLD->KLD[m]+=0.5*(logf(((DistValue1/DistValue2)))*(*DistHist1_PTR)+logf(((DistValue2)/(DistValue1)))*(*DistHist2_PTR)); // Use of the symmetric KLD between the two distributions
//        }
//        delete [] DistHist1_VEC[m];
//        DistHist1_VEC[m]=NULL;
//        delete [] DistHist2_VEC[m];
//        DistHist2_VEC[m]=NULL;
//    }
    return ResultKLD;
}

// Returns for all children of a TreeEM the pairwise MergeKLD in a vector;
vector<MergeKLD *> TreeEM::GetVectorKLDMergeChildren(){
    int numbchild=this->GetNumberChildren();
    vector<MergeKLD *> ResultCompChildrenKLD;
    for (int c1=0; c1<numbchild; c1++) {
        for (int c2=c1+1; c2<numbchild; c2++) {
            MergeKLD * MergeKLDToAdd=this->GetKLDMerge(c1, c2);
            if (MergeKLDToAdd!=NULL) {
                ResultCompChildrenKLD.push_back(MergeKLDToAdd);
            }

        }
    }
    return ResultCompChildrenKLD;
}

vector<MergeKLD*> TreeEM::GetVectorKLDMergeLeaves(){
    vector<MergeKLD *> VectorKLDMergeLeaves;
    vector<MergeKLD *> VectorKLDMergeChildren=this->GetVectorKLDMergeChildren();
    for (int i=0; i<VectorKLDMergeChildren.size(); i++) {
        VectorKLDMergeLeaves.push_back(VectorKLDMergeChildren[i]);
    }
    //VectorKLDMergeChildren.clear();
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        VectorKLDMergeChildren=(this->GetChild(c))->GetVectorKLDMergeLeaves();
        for (int i=0; i<VectorKLDMergeChildren.size(); i++) {
            VectorKLDMergeLeaves.push_back(VectorKLDMergeChildren[i]);
        }
        //VectorKLDMergeChildren.clear();
    }
    return VectorKLDMergeLeaves;
}

MergeKLD * TreeEM::GetToMerge(){
    vector<MergeKLD *> VectorKLDMergeLeaves=this->GetVectorKLDMergeLeaves();
    if (VectorKLDMergeLeaves.size()==0) {
        return NULL;
    }
    else{
        int numbmodal=this->GetNumberModalities();
        MergeKLD * MergeMin=new MergeKLD();
        MergeMin->Parent=NULL;
        MergeMin->Child1=-1;
        MergeMin->Child2=-1;
        MergeMin->KLDTot=1E32;
//        MergeMin->KLD=new float[numbmodal];//{0};
//        for(int i=0;i<numbmodal;i++){
//            MergeMin->KLD[i]=0;
//        }

        for (int i=0; i<VectorKLDMergeLeaves.size(); i++) {
            int Child1=VectorKLDMergeLeaves[i]->Child1;
            int Child2=VectorKLDMergeLeaves[i]->Child2;
            int numbDirectLeaves=VectorKLDMergeLeaves[i]->Parent->GetNumberDirectLeaves();
            float KLDTotTry=VectorKLDMergeLeaves[i]->KLDTot;
//            float * KLDTry=VectorKLDMergeLeaves[i]->KLD;
            // Update if not already checked and is lowest KLDTot
            bool MergeCheckTest=VectorKLDMergeLeaves[i]->Parent->GetMergeCheck()[Child1+Child2*numbDirectLeaves];
            bool MergingPossibility=VectorKLDMergeLeaves[i]->Parent->GetChild(Child1)->GetPriors()==VectorKLDMergeLeaves[i]->Parent->GetChild(Child2)->GetPriors();
            if (MergingPossibility) {
                //                cout<<"Possible merging"<<endl;
                if (!MergeCheckTest&& KLDTotTry<MergeMin->KLDTot) {
                    MergeMin->Parent=VectorKLDMergeLeaves[i]->Parent;
                    MergeMin->Child1=Child1;
                    MergeMin->Child2=Child2;
                    MergeMin->KLDTot=KLDTotTry;
//                    for (int m=0; m<numbmodal; m++) {
//                        MergeMin->KLD[m]=KLDTry[m];
//                    }
                }
            }

            // clear space needed by ith element of VectorKLDMergeLeaves once it has been looked at
            delete VectorKLDMergeLeaves[i];
            VectorKLDMergeLeaves[i]=NULL;
        }
        // if MergeMin not modified <=> no possible merging return NULL
        if (MergeMin->Child1==-1) {
            delete MergeMin;
            return NULL;
        }
        else{
            return MergeMin;
        }
    }
}



// Returns the vector of indices of the children ordered according to increasing KLDTot and availability of leaves to merge
vector<int> TreeEM::OrderingMergingChildren(){
    int numbchild=this->GetNumberChildren();
    vector<GenKLD*> MergeOrder;
    vector<int> OrderMerge;
    for(int c=0;c<numbchild;c++){
        if(this->GetChild(c)->GetNumberAllLeaves()>0){
            MergeOrder.push_back(this->GetChild(c)->GetGenKLD());
            OrderMerge.push_back(c);
        }
    }
    bool flag_swap=1;
    int Csize=MergeOrder.size();
    while(flag_swap==1){
        flag_swap=0;
        for(int c=0;c<Csize-1;c++){
            if(MergeOrder[c]->KLDTot>MergeOrder[c+1]->KLDTot){
                GenKLD * tmp1=MergeOrder[c]->CopyGenKLD();
                GenKLD * tmp2=MergeOrder[c+1]->CopyGenKLD();
                int tmpInd=OrderMerge[c];
                delete MergeOrder[c];
                delete MergeOrder[c+1];
                MergeOrder[c]=tmp2;
                OrderMerge[c]=OrderMerge[c+1];
                MergeOrder[c+1]=tmp1;
                OrderMerge[c+1]=tmpInd;
                flag_swap=1;
            }
        }
    }
    for(int i=0;i<Csize;i++){
        delete MergeOrder[i];
        MergeOrder[i]=NULL;
    }
return OrderMerge;
}

vector<MergeKLD *> TreeEM::OrderingMergingLeaves(){
    int numbleaves=this->GetNumberAllLeaves();
    vector<MergeKLD*> MergeToOrder;
    if(numbleaves>0){
    MergeToOrder=this->GetVectorKLDMergeLeaves();
    int MergeSize=MergeToOrder.size();
    bool flag_swap=1;
    while(flag_swap==1){
        flag_swap=0;
        for(int c=0;c<MergeSize-1;c++){
            if(MergeToOrder[c]->KLDTot>MergeToOrder[c+1]->KLDTot){
                MergeKLD * tmp1=MergeToOrder[c]->CopyMergeKLD();
                MergeKLD * tmp2=MergeToOrder[c+1]->CopyMergeKLD();
                delete MergeToOrder[c];
                MergeToOrder[c]=NULL;
                delete MergeToOrder[c+1];
                MergeToOrder[c+1]=NULL;
                MergeToOrder[c]=tmp2;
                MergeToOrder[c+1]=tmp1;
                flag_swap=1;
            }
        }
    }
    }
    return MergeToOrder;
}

vector<MergeKLD*> TreeEM::GetMergeOrderVertical(){
    vector<int> MergeOrderedChildren=this->OrderingMergingChildren();
    int numbchild=this->GetNumberChildren();
    vector<MergeKLD*> VerticalMergeVector;
    int numbchildmerging=MergeOrderedChildren.size();
    for(int c=0;c<numbchildmerging;c++){
        // if the considered element is a leaf, put directly into the vector

            vector<MergeKLD*> ChildMergeVector=this->GetChild(MergeOrderedChildren[c])->OrderingMergingLeaves();
            VerticalMergeVector.insert(VerticalMergeVector.end(),ChildMergeVector.begin(),ChildMergeVector.end());
            delete ChildMergeVector[c];
            ChildMergeVector[c]=NULL;
    }
    return VerticalMergeVector;
}

//vector<MergeKLD *> TreeEM::GetMergeDiagonal(){
//    vector<MergeKLD *> VectorMergeKLD=this->GetVectorKLDMergeLeaves();
//    bool flag_swap=1;
//    int Csize=VectorMergeKLD.size();
//    while(flag_swap==1){
//        flag_swap=0;
//        for(int c=0;c<Csize-1;c++){
//            if(VectorMergeKLD[c]->KLDTot>VectorMergeKLD[c+1]->KLDTot){
//                MergeKLD * tmp1=VectorMergeKLD[c]->CopyMergeKLD();
//                MergeKLD * tmp2=VectorMergeKLD[c+1]->CopyMergeKLD();
//                delete VectorMergeKLD[c];
//                VectorMergeKLD[c]=NULL;
//                delete VectorMergeKLD[c+1];
//                VectorMergeKLD[c+1]=NULL;
//                VectorMergeKLD[c]=tmp2;
//                VectorMergeKLD[c+1]=tmp1;
//                flag_swap=1;
//            }
//        }
//    }
//    return VectorMergeKLD;
//}

vector<MergeKLD *> TreeEM::GetMergeMoreVertical(int numbclasses){
    int numbchild=this->GetNumberChildren();
    int MaxNumbClasses=6;
    vector<MergeKLD*> ResultMergeVectorList;
    vector<int > ChildrenMergeFirstLevelOrdered=this->OrderingMergingChildren(); // Ordering the general class to split
    int numbposchild=ChildrenMergeFirstLevelOrdered.size();
    vector<MergeKLD *> * ChildrenMergeKLDVec=new vector<MergeKLD*>[numbposchild];
    if(numbclasses>numbposchild || numbclasses<=0){
        return ResultMergeVectorList;
    }
    for(int c=0;c<numbposchild;c++){
        int Child=ChildrenMergeFirstLevelOrdered[c];
        ChildrenMergeKLDVec[c]=this->GetChild(Child)->OrderingMergingLeaves();
    }
    int* ChildrenCombinations=Combination(numbposchild,numbclasses);
    int numbChildrenCombinations=NumbComb(numbclasses,numbposchild);
    for(int i=0;i<numbChildrenCombinations;i++){
        int * TabSize=new int[numbclasses];
        for(int c=0;c<numbclasses;c++){
            TabSize[c]=ChildrenMergeKLDVec[ChildrenCombinations[i*numbclasses+c]].size();
        }
        int * CombinationLeaves=this->CombinationBis(TabSize,numbclasses);
        int prodTabSize=1;
        for(int c=0;c<numbclasses;c++){
            prodTabSize*=TabSize[c];
        }
        for(int j=0;j<prodTabSize;j++){
            for(int c=0;c<numbclasses;c++){
                ResultMergeVectorList.push_back(ChildrenMergeKLDVec[ChildrenCombinations[i*numbclasses+c]][CombinationLeaves[j*numbclasses+c]]);
//                delete [] ChildrenMergeKLDVec[ChildrenCombinations[i*numbclasses+c]][CombinationLeaves[j*numbclasses+c]];
//                ChildrenMergeKLDVec[ChildrenCombinations[i*numbclasses+c]][CombinationLeaves[j*numbclasses+c]]=NULL;
            }
        }
        delete[]TabSize;
        TabSize=NULL;
        delete[]CombinationLeaves;
        CombinationLeaves=NULL;
    }
//    for(int c=0;c<numbposchild;c++){
//        ChildrenMergeKLDVec[c].clear();
//    }
//    delete [] ChildrenMergeKLDVec;
//    ChildrenMergeKLDVec=NULL;
    delete [] ChildrenCombinations;
    ChildrenCombinations=NULL;
return ResultMergeVectorList;
}

void TreeEM::MergeOperation(int Child1, int Child2,SEG_PARAMETERS* segment_param){
    int numbchild=this->GetNumberChildren();
    if (Child1>=numbchild||Child1<0||Child2<0||Child2>=numbchild) {
        cout<<"Index children out of bound, no merging possible"<<endl;
        return;
    }
    if ((this->GetChild(Child1))->GetNumberChildren()!=0|| (this->GetChild(Child2))->GetNumberChildren()!=0) {
        cout << "Nothing other than leaf can be merged"<<endl;
        return;
    }
    if (Child1==Child2) {
        cout<<"Impossible to merge the same child"<<endl;
        return;
    }
    if (Child1>Child2) {// Reorder Child index in order to properly erase them afterwards
        int tmp=Child1;
        Child1=Child2;
        Child2=tmp;
    }
    this->CreateAndAddChildPriors(segment_param,this->GetChild(Child1)->GetPriorsDirect(),1);
    int newnumbchild=this->GetNumberChildren();
    this->GetChild(newnumbchild-1)->NormWeight=(this->GetChild(Child1))->GetNormWeight()+(this->GetChild(Child2))->GetNormWeight();
    Parameters * ParametersMerge=this->ParametersForMerging(Child1,Child2);
    this->GetChild(newnumbchild-1)->SetParameters(ParametersMerge);
    //    delete[] ParametersMerge->ValueParameters;
    //    ParametersMerge->ValueParameters=NULL;
    delete ParametersMerge;
    ParametersMerge=NULL;
    delete this->GetChild(Child1);
    //this->GetChild(Child1)=NULL;
    delete this->GetChild(Child2);
    //this->GetChild(Child2)=NULL;
    //cout<<this->GetChildren().begin()+Child1<<" is the place to look at in the vector"<<endl;
    (this->Children).erase(this->Children.begin()+Child1);
    this->Children.erase(this->Children.begin()+Child2-1);
    newnumbchild=this->GetNumberChildren();
    if (newnumbchild==1) {
        this->CollapseOnlyChild();
    }

}

Parameters * TreeEM::ParametersForMerging(int Child1, int Child2){
    int numbchild=this->GetNumberChildren();
    int numbmodal=this->GetNumberModalities();
    if (Child1>=numbchild||Child1<0||Child2<0||Child2>=numbchild) {
        cout<<"Index children out of bound, no merging possible"<<endl;
        return NULL;
    }
    if ((this->GetChild(Child1))->GetNumberChildren()!=0|| (this->GetChild(Child2))->GetNumberChildren()!=0) {
        cout << "Nothing other than leaf can be merged"<<endl;
        return NULL;
    }
    if (Child1==Child2) {
        cout<<"Impossible to merge the same child"<<endl;
        return NULL;
    }
    if(this->GetChild(Child1)->GetDistributionType()!=this->GetChild(Child2)->GetDistributionType()){
        cout<< "Impossible to merge leaves of different distribution type"<<endl;
        return NULL;
    }
    Parameters * ParametersMerge=new Parameters;
    ParametersMerge->DistributionType=(this->GetChild(Child1))->GetDistributionType();
    ParametersMerge->SizeParameters=(this->GetChild(Child1))->GetSizeParameters();
    ParametersMerge->ValueParameters=new float[ParametersMerge->SizeParameters];//{0};
    for(int i=0;i<ParametersMerge->SizeParameters;i++){
        ParametersMerge->ValueParameters[i]=0;
    }
    // As checked before it is a leaf, necessary of simple distribution type so there is no problem with SizeParameters that is necessarily strictly positive
    switch (ParametersMerge->DistributionType) {


    default: // Gaussian case as default
        float * Mean1_PTR=(this->GetChild(Child1))->GetMean();
        float * Mean2_PTR=(this->GetChild(Child2))->GetMean();
        float * Variance1_PTR=(this->GetChild(Child1))->GetVariance();
        float * Variance2_PTR=(this->GetChild(Child2))->GetVariance();
        float Weight1=this->GetChild(Child1)->GetNormWeight();
        float Weight2=this->GetChild(Child2)->GetNormWeight();
        // Mean for ParametersMerge
        for (int m=0; m<numbmodal; m++) {
            ParametersMerge->ValueParameters[m]=(Weight1*Mean1_PTR[m]+Weight2*Mean2_PTR[m])/(Weight1+Weight2);
        }
        // Variance for ParametersMerge
        for (int m1=0; m1<numbmodal; m1++) {
            for (int m2=0; m2<numbmodal; m2++) {
                ParametersMerge->ValueParameters[numbmodal+m1+m2*numbmodal]=(Weight1*(Variance1_PTR[m1+m2*numbmodal]+(Mean1_PTR[m1]-ParametersMerge->ValueParameters[m1])*(Mean1_PTR[m2]-ParametersMerge->ValueParameters[m2]))+Weight2*(Variance2_PTR[m1+m2*numbmodal]+(Mean2_PTR[m1]-ParametersMerge->ValueParameters[m1])*(Mean2_PTR[m2]-ParametersMerge->ValueParameters[m2])))/(Weight1+Weight2);
            }
        }
        break;
    }
    return ParametersMerge;
}

// When testing a S operation returns the best of the two models obtained.
TreeEM* TreeEM::RunMergeOperation(MergeKLD * MergeTry, bool & AcceptanceDecision,SEG_PARAMETERS * segment_param){
    //    MergeKLD * MergeTry=this->GetToMerge();
    if (MergeTry==NULL) {
        cout<<"No merging to try"<<endl;
        return this;
    }
    // Updating of the checking
    cout<< "Trying merging"<<endl;
    int numbDirectLeaves=(MergeTry->Parent)->GetNumberDirectLeaves();
    (MergeTry->Parent)->GetMergeCheck()[MergeTry->Child1+numbDirectLeaves*MergeTry->Child2]=1;
    (MergeTry->Parent)->GetMergeCheck()[MergeTry->Child2+numbDirectLeaves*MergeTry->Child1]=1;
    // Copy of the current tree
    TreeEM * CopiedTree=this->CopyTree(NULL);
    (MergeTry->Parent)->MergeOperation(MergeTry->Child1,MergeTry->Child2,segment_param);
    this->CollapseOnlyChildTot();
    this->ClearSMChecks();

    // Clear memory for MergeTry
    //    delete [] MergeTry->KLD;
    //    MergeTry->KLD=NULL;


    // Partial EM first
    //    float PartCLL=0;
    //    float PartOldCLL=0;
    //    int PartIteration=0;
    //(MergeTry->Parent)->RunFullEM(PartCLL, PartOldCLL, PartIteration);
    // Run EM
    //    this->UpdateDistribution();
    //    this->UpdateNonNormResp();
    this->UpdateNonNormResp(segment_param);
    this->UpdateNormResp();
    this->UpdateNormRespRoot();
    if(segment_param->flag_MRF){
        this->ClearMRFModel();
        float * GMatrixToSet;
        if(segment_param->flag_GMatrixIn){
            GMatrixToSet=this->PrepareGMatrixFromFile(segment_param->filename_GMatrix);
            if(GMatrixToSet==NULL){
                cout<< "Pb in the preparation of GMatrix"<<endl;
                segment_param->flag_GMatrix=0;
            }
        }
        else{
            GMatrixToSet=this->MRFOptSolveLS();
        }
        this->SetGMatrix(GMatrixToSet);
        if(GMatrixToSet!=NULL){
            delete [] GMatrixToSet;
            GMatrixToSet=NULL;
        }
        this->UpdateMRF();
    }
//    this->SaveAllClasses("/Users/Carole/Documents/PhD/TestMergeResult");
//    this->UpdateNonNormWeights();
//    this->UpdateNormWeights();
    //    CopiedTree->UpdateNormResp();
    //    CopiedTree->UpdateNonNormWeights();
    //    CopiedTree->UpdateNormWeights();
    float CompleteLogLikelihood=0;
    float OldCompleteLogLikelihood=0;
    int Iteration=0;
    if(segment_param->flag_CEM){
        this->RunFullCEM(CompleteLogLikelihood, OldCompleteLogLikelihood,Iteration,segment_param);
    }
    else{
    this->RunFullEM(CompleteLogLikelihood, OldCompleteLogLikelihood,Iteration,segment_param);
    }
//    this->SaveAllClasses("/Users/Carole/Documents/PhD/TestMergeResult");
    // Test of the operation
    float OldCriterion=0;
    float NewCriterion=0;
    if(this->GetDPChildrenDirect()!=NULL){
        if(segment_param->flag_Countmod){
            OldCriterion=CopiedTree->CriterionCalculationMergeDP();
            NewCriterion=this->CriterionCalculationMergeDP();
        }
        else{
            if(!segment_param->flag_DistClassInd){
                OldCriterion=CopiedTree->CriterionCalculationDPInd();
                NewCriterion=this->CriterionCalculationDPInd();
            }
            else{
                OldCriterion=CopiedTree->CriterionCalculationDPNonInd(MergeTry->Parent->FindGeneralClass());
                NewCriterion=this->CriterionCalculationDPNonInd(MergeTry->Parent->FindGeneralClass());
            }

        }
    }
    else{
     OldCriterion=CopiedTree->CriterionCalculation();
     NewCriterion=this->CriterionCalculation();
    }
    delete MergeTry;
    MergeTry=NULL;
    cout<< "The new criterion is "<<NewCriterion<<" and the old criterion is "<<OldCriterion<<endl;
    if (OldCriterion>=NewCriterion) { // Criterion of change not met, delete the modified tree and return the copied one
        delete this;
        AcceptanceDecision=0;
        cout<<"No merging accepted"<<endl;
        return CopiedTree;
    }
    else{ // Criterion of change met delete the not modified tree stored in CopiedTree and return the new model.
        delete CopiedTree;
        AcceptanceDecision=1;
        cout<<"Merging accepted"<<endl;
        this->CollapseOnlyChildTot();
        this->ClearSMChecks();
        return this;
    }
}


// When testing a M operation returns the best of the two models obtained.
TreeEM* TreeEM::RunMergeOperation(vector<MergeKLD *> MergeTry, bool & AcceptanceDecision,SEG_PARAMETERS * segment_param){
    //    MergeKLD * MergeTry=this->GetToMerge();
    int Msize=MergeTry.size();
    if (Msize==0) {
        cout<<"No merging to try"<<endl;
        return this;
    }
    // Updating of the checking
    cout<< "Trying merging"<<endl;

//    int numbDirectLeaves=(MergeTry->Parent)->GetNumberDirectLeaves();
//    (MergeTry->Parent)->GetMergeCheck()[MergeTry->Child1+numbDirectLeaves*MergeTry->Child2]=1;
//    (MergeTry->Parent)->GetMergeCheck()[MergeTry->Child2+numbDirectLeaves*MergeTry->Child1]=1;
    // Copy of the current tree
    TreeEM * CopiedTree=this->CopyTree(NULL);
    for(int c=0;c<Msize;c++){
        (MergeTry[c]->Parent)->MergeOperation(MergeTry[c]->Child1,MergeTry[c]->Child2,segment_param);
        this->CollapseOnlyChildTot();
    }
    this->ClearSMChecks();

    // Clear memory for MergeTry
    //    delete [] MergeTry->KLD;
    //    MergeTry->KLD=NULL;


    // Partial EM first
    //    float PartCLL=0;
    //    float PartOldCLL=0;
    //    int PartIteration=0;
    //(MergeTry->Parent)->RunFullEM(PartCLL, PartOldCLL, PartIteration);
    // Run EM
    //    this->UpdateDistribution();
    //    this->UpdateNonNormResp();
    this->UpdateNonNormResp(segment_param);
    this->UpdateNormResp();
    this->UpdateNormRespRoot();
    if(segment_param->flag_MRF){
        this->ClearMRFModel();
        float * GMatrixToSet;
        if(segment_param->flag_GMatrixIn){
            GMatrixToSet=this->PrepareGMatrixFromFile(segment_param->filename_GMatrix);
            if(GMatrixToSet==NULL){
                cout<< "Pb in the preparation of GMatrix"<<endl;
                segment_param->flag_GMatrix=0;
            }
        }
        else{
            GMatrixToSet=this->MRFOptSolveLS();
        }
        this->SetGMatrix(GMatrixToSet);
        if(GMatrixToSet!=NULL){
            delete [] GMatrixToSet;
            GMatrixToSet=NULL;
        }
        this->UpdateMRF();
    }
//    this->SaveAllClasses("/Users/Carole/Documents/PhD/TestMergeResult");
//    this->UpdateNonNormWeights();
//    this->UpdateNormWeights();
    //    CopiedTree->UpdateNormResp();
    //    CopiedTree->UpdateNonNormWeights();
    //    CopiedTree->UpdateNormWeights();
    float CompleteLogLikelihood=0;
    float OldCompleteLogLikelihood=0;
    int Iteration=0;
    if(segment_param->flag_CEM){
        this->RunFullCEM(CompleteLogLikelihood, OldCompleteLogLikelihood,Iteration,segment_param);
    }
    else{
    this->RunFullEM(CompleteLogLikelihood, OldCompleteLogLikelihood,Iteration,segment_param);
    }
//    this->SaveAllClasses("/Users/Carole/Documents/PhD/TestMergeResult");
    // Test of the operation
    float OldCriterion=0;
    float NewCriterion=0;
    if(this->GetDPChildrenDirect()!=NULL){
        if(segment_param->flag_Countmod){
            OldCriterion=CopiedTree->CriterionCalculationMergeDP();
            NewCriterion=this->CriterionCalculationMergeDP();
        }
        else{
            if(segment_param->flag_DistClassInd){
                OldCriterion=CopiedTree->CriterionCalculationDPInd();
                NewCriterion=this->CriterionCalculationDPInd();
            }
//            else{
//                OldCriterion=CopiedTree->CriterionCalculationDPNonInd(MergeTry->Parent->FindGeneralClass());
//                NewCriterion=this->CriterionCalculationDPNonInd(MergeTry->Parent->FindGeneralClass());
//            }

        }
    }
    else{
     OldCriterion=CopiedTree->CriterionCalculation();
     NewCriterion=this->CriterionCalculation();
    }
//    delete MergeTry;
//    MergeTry=NULL;
    cout<< "The new criterion is "<<NewCriterion<<" and the old criterion is "<<OldCriterion<<endl;
    if (OldCriterion>=NewCriterion) { // Criterion of change not met, delete the modified tree and return the copied one
        delete this;
        AcceptanceDecision=0;
        cout<<"No merging accepted"<<endl;
        return CopiedTree;
    }
    else{ // Criterion of change met delete the not modified tree stored in CopiedTree and return the new model.
        delete CopiedTree;
        AcceptanceDecision=1;
        cout<<"Merging accepted"<<endl;
        this->CollapseOnlyChildTot();
        this->ClearSMChecks();
        return this;
    }
}

/******** METHODS FOR SPLITTING OPERATIONS *************************/
SplitKLD * TreeEM::GetKLDSplit(int ChildInput){
    // First check if input integer is correct
    int numbchild=this->GetNumberChildren();
    if (ChildInput<0||ChildInput>=numbchild) {
        cout<<"Out of bound for child to split"<<endl;
        return NULL;
    }
    // then check if can be split <=> check if split a leaf
    if (this->GetChild(ChildInput)->GetNumberChildren()!=0) {
        cout<<"It is not a leaf so cannot be split"<<endl;
        return NULL;
    }
    SplitKLD * SplitKLDResult=new SplitKLD();
    SplitKLDResult->Parent=this;
    SplitKLDResult->ChildToSplit=ChildInput;
    SplitKLDResult->KLDTot=0;
    int numbmodal=this->GetNumberModalities();
    int binspowmodal=(int)pow_int(numbbins,numbmodal);
    SplitKLDResult->KLD=new float[MaxNumbModal];//{0};
    for(int i=0;i<MaxNumbModal;i++){
        SplitKLDResult->KLD[i]=0;
    }

    // Calculating KLDTot and deleting the corresponding float arrays needed after use
    float * DistHistTotal=(this->GetChild(ChildInput))->GetDistHistogramTotal();
    float * DistHistTotal_PTR=DistHistTotal;
    float * DataHistTotal=(this->GetChild(ChildInput))->GetDataHistogramTotal();
    float * DataHistTotal_PTR=DataHistTotal;
    float SumDist=0;
    float SumData=0;

    PrecisionTYPE tmpSplitKLDTot=0;
    int CountNanSplit=0;

    for (int i=0; i<binspowmodal; i++,DistHistTotal_PTR++,DataHistTotal_PTR++) {
        SumDist+=*DistHistTotal_PTR;
        SumData+=*DataHistTotal_PTR;
        float DistValue=(*DistHistTotal_PTR)<=0.00001?0.00001:(*DistHistTotal_PTR);
        float DataValue=(*DataHistTotal_PTR)<=0.00001?0.00001:(*DataHistTotal_PTR);

        tmpSplitKLDTot+=(PrecisionTYPE)0.5*(log(DataValue/DistValue)*(*DataHistTotal_PTR)+log(DistValue/DataValue)*(*DistHistTotal_PTR)); // Use of the symmetric KLD between the data and the distribution
        if (tmpSplitKLDTot!=tmpSplitKLDTot) {
//            cout<<"Pb at index "<<i<<" : "<<DataValue<<" and "<<DistValue<<endl;
            CountNanSplit++;
        }
    }
//    cout<<"Number of Nan for split KLD is "<<CountNanSplit<<endl;
    SplitKLDResult->KLDTot=(float)tmpSplitKLDTot;
//    cout<<SplitKLDResult->KLDTot<< " and SumData "<<SumData<< "and SumDist are "<<SumDist<<endl;



//    // Calculation of KLDcomp
//    vector<float *> DistHist_VEC=this->GetChild(ChildInput)->GetDistHistogram();
//    vector<float *> DataHist_VEC=this->GetChild(ChildInput)->GetDataHistogram();
//    DistHistTotal_PTR=DistHistTotal;
//    DataHistTotal_PTR=DataHistTotal;

//    //    float * DiffDistHist=new float[numbbins];//{0};
//    //    for(int i=0;i<numbbins;i++){
//    //        DiffDistHist[i]=0;
//    //    }
//    //    float * DiffDataHist=new float[numbbins];//{0};
//    //    for(int i=0;i<numbbins;i++){
//    //        DiffDataHist[i]=0;
//    //    }

//    //    float * DiffDistHist_PTR=DiffDistHist;
//    //    float * DiffDataHist_PTR=DiffDataHist;
//    //    float * DistHist_PTR=DistHist_VEC[0];
//    //    float * DataHist_PTR=DataHist_VEC[0];
//    //    int CountNonZeroDiffDist=0;
//    //    int CountNonZeroDiffData=0;
//    //      SumData=0;
//    //    for (int i=0; i<numbbins; i++,DistHist_PTR++,DistHistTotal_PTR++,DataHistTotal_PTR++,DataHist_PTR++,DiffDataHist_PTR++,DiffDistHist_PTR++) {

//    //        *DiffDistHist_PTR=(*DistHist_PTR)-(*DistHistTotal_PTR);
//    //        *DiffDataHist_PTR=(*DataHist_PTR)-(*DataHistTotal_PTR);
//    //        SumData+=*DistHistTotal_PTR;
//    //        if (fabsf(*DiffDistHist_PTR)>1E-6) {
//    //            CountNonZeroDiffDist++;
//    //            //cout<<*DiffDistHist_PTR<<endl;
//    //        }
//    //        if (fabsf(*DiffDataHist_PTR)>1E-6) {
//    //            CountNonZeroDiffData++;
//    //        }
//    //    }
//    //    delete [] DiffDataHist;
//    //    DiffDataHist=NULL;
//    //    delete [] DiffDistHist;
//    //    DiffDistHist=NULL;

//    for (int m=0; m<numbmodal; m++) {
//        float * DistHist_PTR=DistHist_VEC[m];
//        float * DataHist_PTR=DataHist_VEC[m];
//        for (int i=0; i<numbbins; i++,DistHist_PTR++,DataHist_PTR++) {
//            float DistValue=(*DistHist_PTR)<=0.00001?0.00001:(*DistHist_PTR);

//            float DataValue=(*DataHist_PTR)<=0.00001?0.00001:(*DataHist_PTR);

//            SplitKLDResult->KLD[m]+=0.5*(logf(DataValue/DistValue)*(*DataHist_PTR)+logf(DistValue/DataValue)*(*DistHist_PTR)); // Use of the symmetric KLD between the two distributions
//        }
//        delete [] DistHist_VEC[m];
//        DistHist_VEC[m]=NULL;
//        delete [] DataHist_VEC[m];
//        DataHist_VEC[m]=NULL;
//    }
    delete [] DistHistTotal;
    DistHistTotal=NULL;
    delete [] DataHistTotal;
    DataHistTotal=NULL;

    return SplitKLDResult;
}

// Returns the vector containing the pointers to all the SplitKLD of the children of the considered tree. Used to obtain ordered list of Splitting to try.
vector<SplitKLD *> TreeEM::GetVectorKLDSplitChildren(SEG_PARAMETERS * segment_param){
    int numbchild=this->GetNumberChildren();
    vector<SplitKLD *> KLDSplitChildren;
    for (int c=0; c<numbchild; c++) {
        SplitKLD * SplitKLDToPush=this->GetKLDSplit(c);
        if(this->GetChild(c)->GetDistributionType()==2){
            SplitKLDToPush->TypeChange=segment_param->uniformTypeChange;
        }

        KLDSplitChildren.push_back(SplitKLDToPush);
    }
    return KLDSplitChildren;
}

// Returns the vector of pointers to all SplitKLD structures for the leaves of a given tree.
vector<SplitKLD *> TreeEM::GetVectorKLDSplitLeaves(SEG_PARAMETERS * segment_param){
    vector<SplitKLD *> KLDSplitLeaves;
    //cout<<KLDSplitLeaves.size()<<endl;
    int numbchild=this->GetNumberChildren();
    for (int c=0; c<numbchild; c++) {
        if (this->GetChild(c)->GetNumberChildren()==0) {
            //cout<<"Leaf for KLD calculation"<<endl;
            SplitKLD * SplitKLDToPush=this->GetKLDSplit(c);
            if(this->GetChild(c)->GetDistributionType()==2){ // case where we try to change a uniform distribution split into 1 gaussian (type change 0) or 2 gaussians (type change 2)
//                SplitKLD * SplitKLDToPush_2=SplitKLDToPush->CopySplitKLD();
//                SplitKLDToPush_2->TypeChange=2;
//                SplitKLDToPush->TypeChange=0;
//                KLDSplitLeaves.push_back(SplitKLDToPush_2);
                SplitKLDToPush->TypeChange=segment_param->uniformTypeChange;
            }
            if(this->GetChild(c)->GetNormWeight()>10E-6){ // Possible to split only if weight is high enough
                KLDSplitLeaves.push_back(SplitKLDToPush);
            }

        }
        else{
            vector<SplitKLD*> KLDSplitLeavesChild=this->GetChild(c)->GetVectorKLDSplitLeaves(segment_param);
            for (int i=0; i<KLDSplitLeavesChild.size(); i++) {
                KLDSplitLeaves.push_back(KLDSplitLeavesChild[i]);
            }
        }
    }
    //cout<<"Size of KLDVector before return "<<KLDSplitLeaves.size()<<endl;
    return KLDSplitLeaves;
}

vector<SplitKLD*> TreeEM::GetToSplitMore(vector<int> ClassesToSplit, SEG_PARAMETERS * segment_param){
    int numbchild=this->GetNumberChildren();
    int sizeVec=ClassesToSplit.size();
    vector<SplitKLD*> SplitKLDVec;
    if(sizeVec>numbchild){
        return SplitKLDVec;
    }
    // first check if the classes we are trying to split really do exist
    for(int i=0;i<sizeVec;i++){
        if (ClassesToSplit[i]>=numbchild){
            return SplitKLDVec;
        }
    }
    // Now for each of the classes, get the max SplitKLD or the split KLD of this class itself if it is a leaf
    for(int i=0;i<sizeVec;i++){
        if(this->GetChild(ClassesToSplit[i])->IsLeaf()){
            SplitKLDVec.push_back(this->GetKLDSplit(ClassesToSplit[i]));
        }
        else{
            SplitKLDVec.push_back(this->GetChild(ClassesToSplit[i])->GetToSplit(segment_param));
        }
    }
    return SplitKLDVec;

}

// Returns the vector of SplitKLD pointers ordered according to decreasing KLDTot
vector<int> TreeEM::OrderingSplittingChildren(){
    int numbchild=this->GetNumberChildren();
    vector<GenKLD*> SplitTryToOrder;
    vector<int> OrderedChildren;
    for(int c=0;c<numbchild;c++){
        if(this->GetChild(c)->GetNumberAllLeaves()<MaxSupport && this->GetChild(c)->GetNormWeight()>WeightThreshold){
            SplitTryToOrder.push_back(this->GetChild(c)->GetGenKLD());
            OrderedChildren.push_back(c);
        }
    }
    bool flag_swap=1;
    int Csize=SplitTryToOrder.size();
    while(flag_swap==1){
        flag_swap=0;
        for(int c=0;c<Csize-1;c++){
            if(SplitTryToOrder[c]->KLDTot<SplitTryToOrder[c+1]->KLDTot){
                GenKLD * tmp1=SplitTryToOrder[c]->CopyGenKLD();
                GenKLD * tmp2=SplitTryToOrder[c+1]->CopyGenKLD();
                int tmpInd=OrderedChildren[c];
                delete SplitTryToOrder[c];
                delete SplitTryToOrder[c+1];
                SplitTryToOrder[c]=tmp2;
                OrderedChildren[c]=OrderedChildren[c+1];
                SplitTryToOrder[c+1]=tmp1;
                OrderedChildren[c+1]=tmpInd;
                flag_swap=1;
            }
        }
    }
    for(int c=0;c<Csize;c++){
        delete SplitTryToOrder[c];
        SplitTryToOrder[c]=NULL;
    }
return OrderedChildren;
}

vector<SplitKLD *> TreeEM::OrderingSplittingLeaves(SEG_PARAMETERS * segment_param){
    int numbleaves=this->GetNumberAllLeaves();
    vector<SplitKLD*> SplitTryToOrder;
    if(numbleaves>0){
    SplitTryToOrder=this->GetVectorKLDSplitLeaves(segment_param);
    int numbSplit=SplitTryToOrder.size();
    bool flag_swap=1;
    while(flag_swap==1){
        flag_swap=0;
        for(int c=0;c<numbSplit-1;c++){
            if(SplitTryToOrder[c]->KLDTot<SplitTryToOrder[c+1]->KLDTot){
                SplitKLD * tmp1=SplitTryToOrder[c]->CopySplitKLD();
                SplitKLD * tmp2=SplitTryToOrder[c+1]->CopySplitKLD();
                delete SplitTryToOrder[c];
                SplitTryToOrder[c]=NULL;
                delete SplitTryToOrder[c+1];
                SplitTryToOrder[c+1]=NULL;
                SplitTryToOrder[c]=tmp2;
                SplitTryToOrder[c+1]=tmp1;
                flag_swap=1;
            }
        }
    }
    }
    else{
        int Child=this->GetParent()->FindIndex(this);
        SplitTryToOrder.push_back(this->GetParent()->GetKLDSplit(Child));
    }
    return SplitTryToOrder;
}

vector<SplitKLD*> TreeEM::GetSplitOrderVertical(SEG_PARAMETERS * segment_param){
    vector<int> SplitOrderedChildren=this->OrderingSplittingChildren();
    int numbchild=this->GetNumberChildren();
    vector<SplitKLD*> VerticalSplitVector;
    int numbposchild=SplitOrderedChildren.size();
    for(int c=0;c<numbposchild;c++){
        // if the considered element is a leaf, put directly into the vector
        if(this->GetChild(SplitOrderedChildren[c])->IsLeaf()){
            VerticalSplitVector.push_back(this->GetKLDSplit(SplitOrderedChildren[c]));
//            delete SplitOrderedChildren[c];
//            SplitOrderedChildren[c]=NULL;
        }
        else {
            vector<SplitKLD*> ChildSplitVector=this->GetChild(SplitOrderedChildren[c])->OrderingSplittingLeaves(segment_param);
            VerticalSplitVector.insert(VerticalSplitVector.end(),ChildSplitVector.begin(),ChildSplitVector.end());
//            delete SplitOrderedChildren[c];
//            SplitOrderedChildren[c]=NULL;
        }
    }
    return VerticalSplitVector;
}


vector<SplitKLD *> TreeEM::GetSplitMoreVertical(int numbclasses,SEG_PARAMETERS * segment_param){
    int numbchild=this->GetNumberChildren();
    int MaxNumbClasses=6;
    vector<SplitKLD*> ResultVectorSplitList;
    vector<int > ChildrenSplitFirstLevelOrdered=this->OrderingSplittingChildren(); // Ordering the general class to split
    int numbposchild=ChildrenSplitFirstLevelOrdered.size();
    vector<SplitKLD *> * ChildrenSplitKLDVec=new vector<SplitKLD*>[numbposchild];
    if(numbclasses>numbposchild || numbclasses<=0){
        return ResultVectorSplitList;
    }
    for(int c=0;c<numbchild;c++){
        int Child=ChildrenSplitFirstLevelOrdered[c];
        ChildrenSplitKLDVec[c]=this->GetChild(Child)->OrderingSplittingLeaves(segment_param);
    }
    int* ChildrenCombinations=Combination(numbposchild,numbclasses);
    int numbChildrenCombinations=NumbComb(numbclasses,numbposchild);
    for(int i=0;i<numbChildrenCombinations;i++){
        int * TabSize=new int[numbclasses];
        for(int c=0;c<numbclasses;c++){
            TabSize[c]=ChildrenSplitKLDVec[ChildrenCombinations[i*numbclasses+c]].size();
        }
        int * CombinationLeaves=this->CombinationBis(TabSize,numbclasses);
        int prodTabSize=1;
        for(int c=0;c<numbclasses;c++){
            prodTabSize*=TabSize[c];
        }
        for(int j=0;j<prodTabSize;j++){
            for(int c=0;c<numbclasses;c++){
//                cout<< ChildrenCombinations[i*numbclasses+c]<< " "<<CombinationLeaves[j*numbclasses+c]<<endl;
                ResultVectorSplitList.push_back(ChildrenSplitKLDVec[ChildrenCombinations[i*numbclasses+c]][CombinationLeaves[j*numbclasses+c]]);
//                delete [] ChildrenSplitKLDVec[ChildrenCombinations[i*numbclasses+c]][CombinationLeaves[j*numbclasses+c]];
//                ChildrenSplitKLDVec[ChildrenCombinations[i*numbclasses+c]][CombinationLeaves[j*numbclasses+c]]=NULL;
            }
        }

        delete[]TabSize;
        TabSize=NULL;
        delete[]CombinationLeaves;
        CombinationLeaves=NULL;
    }
//    for(int c=0;c<numbchild;c++){
//        if (ChildrenSplitKLDVec[c]!=NULL){
//            delete [] ChildrenSplitKLDVec[c];
//            ChildrenSplitKLDVec[c]=NULL;
//        }
//    }
//    for(int c=0;c<numbposchild;c++){
//        ChildrenSplitKLDVec[c].clear();
//    }
//    delete [] ChildrenSplitKLDVec;
//    ChildrenSplitKLDVec=NULL;
    if(ChildrenCombinations!=NULL){
        delete[] ChildrenCombinations;
        ChildrenCombinations=NULL;
    }
return ResultVectorSplitList;
}

int * TreeEM::Combination(int n,int k){
    if(n<k){
        cout<<"Impossible to take more than what is available"<<endl;
        return NULL;
    }
    int * Combinations=new int[k*NumbComb(k,n)];
    int SizeComb=k*NumbComb(k,n);
//    int Comb_tmp[k];
    int * Comb_tmp=new int[k];
    for(int i=0;i<k;i++){
        Comb_tmp[i]=i;
    }
    int Ind=0;
    for(int i=0;i<k;i++){
        Combinations[Ind+i]=Comb_tmp[i];
    }
    Ind=Ind+k;
    while(Comb_tmp[0]<=n-k && Ind<SizeComb){
        int i=k-1;
        bool flag_change=0;
        while(i>=0 && !flag_change){
            if(Comb_tmp[i]<n-k+i){
                Comb_tmp[i]++;
                flag_change=1;
            }
            i--;
        }
        for(int j=0;j<k;j++){
            Combinations[Ind+j]=Comb_tmp[j];
        }
        Ind+=k;
    }
    if(Ind<SizeComb){
    Comb_tmp[0]++;
    for(int j=0;j<k;j++){
        Combinations[Ind+j]=Comb_tmp[j];
    }
    }
    if(Comb_tmp!=NULL){
        delete[] Comb_tmp;
        Comb_tmp=NULL;
    }
    return Combinations;
}

int * TreeEM::CombinationBis(int * TabSize,int numbclasses){
    // Calculation of the number of possible combinations and initialisation of the array
    int prodTabSize=1;
//    int Shift[numbclasses];
    int * Shift=new int[numbclasses];
    for(int i=0;i<numbclasses;i++){
        prodTabSize*=TabSize[i];
        Shift[i]=1;
    }
    for(int c=numbclasses-1;c>=0;c--){
        if(c<numbclasses-1){
            Shift[c]=Shift[c+1]*TabSize[c+1];
        }
    }
    int * Combination=new int[numbclasses*prodTabSize];
    int Ind=0;
    for(int i=0;i<prodTabSize;i++){
        int tmpind=i;
        for(int c=0;c<numbclasses;c++){
            Combination[Ind+c]=tmpind/Shift[c];
            tmpind-=Combination[Ind+c]*Shift[c];
        }
        Ind+=numbclasses;
    }
    if(Shift!=NULL){
        delete [] Shift;
        Shift=NULL;
    }
cout<<Ind;
return Combination;
}

inline int TreeEM::Factorial(int x) {
  return (x <= 1 ? 1 : x * Factorial(x - 1));
}

int TreeEM::NumbComb(int k, int n){
    if(k>n){
        return 0;
    }
    cout<< n << " "<<k<<" "<<Factorial(n-k);
//    printf("the factorial components are %d %d %d \n",Factorial(n),Factorial(n-k),Factorial(k));
    return Factorial(n)/(Factorial(n-k)*Factorial(k));
}

SplitKLD * TreeEM::GetToSplit(SEG_PARAMETERS * segment_param){
    vector<SplitKLD *> VectorSplit( this->GetVectorKLDSplitLeaves(segment_param));
    int numbLeaves=VectorSplit.size();
    //cout<<VectorSplit.size()<<"is the size of Vector Split"<<endl;
    if (numbLeaves<=0) {
        return NULL;
    }
    // Initialisation of the result
    SplitKLD * SplitMax=new SplitKLD;
    SplitMax->KLDTot=-1E32;
    int numbmodal=this->GetNumberModalities();
    SplitMax->KLD=new float[numbmodal];//{-1E32};
    for(int i=0;i<numbmodal;i++){
        SplitMax->KLD[i]=-1E32;
    }
    for (int i=0; i<numbLeaves; i++) {
        bool SplitCheckTest=(VectorSplit[i]->Parent)->GetSplitCheck()[VectorSplit[i]->ChildToSplit];
        if (!SplitCheckTest&& SplitMax->KLDTot<VectorSplit[i]->KLDTot) { // if the leave has not been checked yet for the split and the KLD is more than the previous one, we update SplitMax by copying the structure of VectorSplit[i]
            SplitMax->Parent=VectorSplit[i]->Parent;
            SplitMax->ChildToSplit=VectorSplit[i]->ChildToSplit;
            SplitMax->KLDTot=VectorSplit[i]->KLDTot;
            for (int m=0; m<numbmodal; m++) {
                SplitMax->KLD[m]=VectorSplit[i]->KLD[m];
            }
        }
        // Delete and free the considered element of VectorSplit.
        //delete [] VectorSplit[i]->KLD;
        delete VectorSplit[i];
        VectorSplit[i]=NULL;
    }
    if (SplitMax->KLDTot<0) {
        delete SplitMax;
        SplitMax=NULL;
    }
    return SplitMax;
}

// Performs the Split operation and modify accordingly the tree
void TreeEM::SplitOperation(SplitKLD * SplitTry, int choiceInit,SEG_PARAMETERS * segment_param){
    // Check if it is effectively a child that can be split
    int numbchild=SplitTry->Parent->GetNumberChildren();
    int numbmodal=this->GetNumberModalities();
    int ChildToSplit=SplitTry->ChildToSplit;
    TreeEM* ParentToSplit=SplitTry->Parent;
    if (ChildToSplit>=numbchild||ChildToSplit<0) {
        cout<<"Child to split out of bounds"<<endl;
        return;
    }
    if (ParentToSplit->GetChild(ChildToSplit)->GetNumberChildren()!=0) {
        cout<<"Impossible to split something else than leaf"<<endl;
        return;
    }
    switch(SplitTry->TypeChange){
    case 0:{ // case where the split transforms the uniform distribution into a unique Gaussian
        float * Mean=ParentToSplit->GetChild(ChildToSplit)->GetMeanDirect();
        float * Variance=ParentToSplit->GetChild(ChildToSplit)->GetDiagVarianceDirect();
        // Building of the parameters of newly formed Gaussian
        Parameters * ParametersSplit=new Parameters();
        ParametersSplit->DistributionType=1;
        ParametersSplit->SizeParameters=CalculateSizeParameters(ParametersSplit->DistributionType);
        ParametersSplit->ValueParameters=new float[ParametersSplit->SizeParameters];//{0};
        for(int i=0;i<ParametersSplit->SizeParameters;i++){
            ParametersSplit->ValueParameters[i]=0;
        }
        for(int m=0;m<numbmodal;m++){
            ParametersSplit->ValueParameters[m]=Mean[m];
        }
        for(int m1=0;m1<numbmodal;m1++){
            for(int m2=0;m2<numbmodal;m2++){
                if(m1==m2){
                    ParametersSplit->ValueParameters[numbmodal+m1*numbmodal+m2]=Variance[m1];
                }
            }
        }
        ParentToSplit->GetChild(ChildToSplit)->SetParameters(ParametersSplit);
        // Clearing memory used for parameters construction
        delete [] Mean;
        Mean=NULL;
        delete [] Variance;
        Variance = NULL;
        delete ParametersSplit;
        ParametersSplit=NULL;
        break;
    }
    case 2:{//or uniform giving two gaussian at first step
                Parameters ** ParametersSplit=(ParentToSplit->GetChild(ChildToSplit))->ParametersDoubleForSplittingUniform( choiceInit);
                (ParentToSplit->GetChild(ChildToSplit))->CreateAndAddChildWeight(segment_param,0.5,1);
                (ParentToSplit->GetChild(ChildToSplit))->CreateAndAddChildWeight(segment_param,0.5,1);
                (ParentToSplit->GetChild(ChildToSplit)->GetChild(0))->SetParameters(ParametersSplit[0]);
                (ParentToSplit->GetChild(ChildToSplit)->GetChild(1))->SetParameters(ParametersSplit[1]);

                // Clearing memory dedicated to Parameters split
                delete ParametersSplit[0];
                ParametersSplit[0]=NULL;
                delete ParametersSplit[1];
                ParametersSplit[1]=NULL;
                delete ParametersSplit;
                ParametersSplit=NULL;
        //        this->GetChild(ChildToSplit)->GetChild(0)->SetNormWeight(0.5);
        //        this->GetChild(ChildToSplit)->GetChild(1)->SetNormWeight(0.5);
                //    this->PutAllLeavesToChildrenLevel();
                break;
                return;
    }
    case 3:{ // uniform giving one gaussian + uniform
        Parameters * ParametersSplit=(ParentToSplit->GetChild(ChildToSplit))->ParametersForSplittingUniform3();
        nifti_image * PriorsToAdd=ParentToSplit->GetChild(ChildToSplit)->CreatePriorsFromAdaptedPriors();
        nifti_image * PriorsConstant=ParentToSplit->BuildConstantPriors(segment_param->OutliersWeight);
        ParentToSplit->GetChild(ChildToSplit)->SetPriors(PriorsToAdd);
        ParentToSplit->GetChild(ChildToSplit)->SetParameters(ParametersSplit);
        ParentToSplit->CreateAndAddChildPriors(segment_param,PriorsConstant,2);
        ParentToSplit->NormalisePriors();
        ParentToSplit->NormalisePriorsAdapted();






        // Clearing memory dedicated to Parameters split
        delete ParametersSplit;
        ParametersSplit=NULL;
        break;
        return;

    }
    case 4:{ // uniform giving one gaussian and one uniform under same priors
       Parameters * ParametersSplit=(ParentToSplit->GetChild(ChildToSplit))->ParametersForSplittingUniform3();
       (ParentToSplit->GetChild(ChildToSplit))->CreateAndAddChildWeight(segment_param,0.9999,1);
       (ParentToSplit->GetChild(ChildToSplit)->GetChild(0))->SetParameters(ParametersSplit);
               (ParentToSplit->GetChild(ChildToSplit))->CreateAndAddChildWeight(segment_param,0.0001,2);

       // Clearing memory dedicated to Parameters split
       delete ParametersSplit;
       ParametersSplit=NULL;
       break;
       return;
    }

    default :{ // classical case where the split happens on a gaussian distributed class
        Parameters ** ParametersSplit=(ParentToSplit->GetChild(ChildToSplit))->ParametersForSplitting( choiceInit);
        (ParentToSplit->GetChild(ChildToSplit))->CreateAndAddChildWeight(segment_param,0.5,1);
        (ParentToSplit->GetChild(ChildToSplit))->CreateAndAddChildWeight(segment_param,0.5,1);
        (ParentToSplit->GetChild(ChildToSplit)->GetChild(0))->SetParameters(ParametersSplit[0]);
        (ParentToSplit->GetChild(ChildToSplit)->GetChild(1))->SetParameters(ParametersSplit[1]);

        // Clearing memory dedicated to Parameters split
        delete ParametersSplit[0];
        ParametersSplit[0]=NULL;
        delete ParametersSplit[1];
        ParametersSplit[1]=NULL;
        delete ParametersSplit;
        ParametersSplit=NULL;
//        this->GetChild(ChildToSplit)->GetChild(0)->SetNormWeight(0.5);
//        this->GetChild(ChildToSplit)->GetChild(1)->SetNormWeight(0.5);
        //    this->PutAllLeavesToChildrenLevel();
        return;
    }
    }

    return;
}

// Performs the Split operation and modify accordingly the tree
void TreeEM::SplitOperation(int ChildToSplit, int choiceInit,SEG_PARAMETERS * segment_param){
    // Check if it is effectively a child that can be split
    int numbchild=this->GetNumberChildren();
    if (ChildToSplit>=numbchild||ChildToSplit<0) {
        cout<<"Child to split out of bounds"<<endl;
        return;
    }
    if (this->GetChild(ChildToSplit)->GetNumberChildren()!=0) {
        cout<<"Impossible to split something else than leaf"<<endl;
        return;
    }
    Parameters ** ParametersSplit=(this->GetChild(ChildToSplit))->ParametersForSplitting( choiceInit);
    (this->GetChild(ChildToSplit))->CreateAndAddChildPriors(segment_param,NULL,1);
    (this->GetChild(ChildToSplit))->CreateAndAddChildPriors(segment_param,NULL,1);
    (this->GetChild(ChildToSplit)->GetChild(0))->SetParameters(ParametersSplit[0]);
    (this->GetChild(ChildToSplit)->GetChild(1))->SetParameters(ParametersSplit[1]);

    // Clearing memory dedicated to Parameters split
    delete ParametersSplit[0];
    ParametersSplit[0]=NULL;
    delete ParametersSplit[1];
    ParametersSplit[1]=NULL;
    delete ParametersSplit;
    ParametersSplit=NULL;
    this->GetChild(ChildToSplit)->GetChild(0)->SetNormWeight(0.5);
    this->GetChild(ChildToSplit)->GetChild(1)->SetNormWeight(0.5);
    //    this->PutAllLeavesToChildrenLevel();
    return;
}

Parameters * TreeEM::ParametersForSplittingUniform3(){
    if(this->GetDistributionType()!=2){ // check if we are trying to split a uniform distribution
        return NULL;
    }
    int numbmodal=this->GetNumberModalities();
    float * DataHistogram=this->GetDataHistogramTotal();
    vector<int> dimHist;
    for(int m=0;m<numbmodal;m++){
        dimHist.push_back(numbbins);
    }
    float sizeBins=1.0/numbbins;
    float * DataHistogramBlurred=GaussianBlurring(DataHistogram,1,dimHist);
    float maxHist=-1;
    int numberHistElements=pow_int(numbbins,numbmodal);
    int maxInd=numbbins/2;
    for(int i=0;i<numberHistElements;i++){
        if(DataHistogramBlurred[i]>maxHist){
            maxInd=i;
            maxHist=DataHistogramBlurred[i];
        }
    }
    // Temporary writing of data histogram into image
//    nifti_image * TmpHistImage=nifti_copy_nim_info(this->GetDataImage());
//    TmpHistImage->dim[0]=3;
//    TmpHistImage->dim[4]=1;
//    TmpHistImage->dim[1]=numbbins;
//    TmpHistImage->dim[2]=numbbins;
//    TmpHistImage->dim[3]=1;
//    nifti_update_dims_from_array(TmpHistImage);
//    TmpHistImage->data=(void *)calloc(TmpHistImage->nvox,sizeof(float));
//    float * TmpHistImageData_PTR=static_cast<float *>(TmpHistImage->data);
//    for(int i=0;i<numbbins*numbbins;i++,TmpHistImageData_PTR++){
//        *TmpHistImageData_PTR=DataHistogramBlurred[i];
//    }
//    nifti_set_filenames(TmpHistImage,"/Users/Carole/Documents/PhD/BlurredHist.nii.gz",0,0);
//    nifti_image_write(TmpHistImage);
//    nifti_image_free(TmpHistImage);

    // Transform maxInd found into Index for each modality
    int * MaxIndPerModal=new int[numbmodal];
    for(int m=0;m<numbmodal;m++){
        MaxIndPerModal[m]=0;
    }
    int tmp=maxInd;
    for(int m=numbmodal-1;m>=0;m--){
        MaxIndPerModal[m]=tmp/(int)pow_int(numbbins,m);
        tmp=tmp-MaxIndPerModal[m]*pow_int(numbbins,m);
    }
    // create the parameters corresponding to the gaussian to create
    Parameters * ParametersSplit=new Parameters();
    ParametersSplit->DistributionType=1;
    ParametersSplit->SizeParameters=CalculateSizeParameters(ParametersSplit->DistributionType);
    ParametersSplit->ValueParameters=new float[ParametersSplit->SizeParameters];
    //Initialise parameters at value 0:
    for(int m=0;m<numbmodal+numbmodal*numbmodal;m++){
        ParametersSplit->ValueParameters[m]=0;
    }
//    MaxIndPerModal[0]=15;
//    MaxIndPerModal[1]=55;

    // Transform the index max to form the mean of the parameters
    float * Variance=this->GetVarianceDirect();
    float * Mean = new float[numbmodal];
    for(int m=0;m<numbmodal;m++){
        Mean[m]=MaxIndPerModal[m]*sizeBins;
    }
    // copy values into value parameters
    for(int m1=0;m1<numbmodal;m1++){
        ParametersSplit->ValueParameters[m1]=Mean[m1];
        for(int m2=0;m2<numbmodal;m2++){
            ParametersSplit->ValueParameters[numbmodal+m1+m2*numbmodal]=Variance[m1+numbmodal*m2]/10;
        }
    }
    // clearing memory
    delete [] Variance;
    Variance = NULL;
    delete [] Mean;
    Mean = NULL;
    delete[] MaxIndPerModal;
    MaxIndPerModal=NULL;
    delete [] DataHistogram;
    DataHistogram=NULL;
    delete [] DataHistogramBlurred;
    DataHistogramBlurred = NULL;
    dimHist.clear();

    return ParametersSplit;

}

Parameters ** TreeEM::ParametersDoubleForSplittingUniform(int choiceInit){
    // Check if it is really a leaf to split and if it is really a uniform distribution
    if (this->GetNumberChildren()!=0 || this->GetDistributionType()!=2) {
        return NULL;
    }
    Parameters ** ParametersSplit=new Parameters* [2];
    ParametersSplit[0]=new Parameters();
    ParametersSplit[1]=new Parameters();
    ParametersSplit[0]->DistributionType=1;
    ParametersSplit[1]->DistributionType=1;
    ParametersSplit[0]->SizeParameters=CalculateSizeParameters(ParametersSplit[0]->DistributionType);
    ParametersSplit[1]->SizeParameters=CalculateSizeParameters(ParametersSplit[1]->DistributionType);
    ParametersSplit[0]->ValueParameters=new float[ParametersSplit[0]->SizeParameters];//{0};
    for(int i=0;i<ParametersSplit[0]->SizeParameters;i++){
        ParametersSplit[0]->ValueParameters[i]=0;
    }
    ParametersSplit[1]->ValueParameters=new float[ParametersSplit[1]->SizeParameters];//{0};
    for(int i=0;i<ParametersSplit[1]->SizeParameters;i++){
        ParametersSplit[1]->ValueParameters[i]=0;
    }
    //    float * ParametersMean=this->GetMean();
    //    float * ParametersVariance=this->GetVariance();

    // Temporary values for splitting initialisation :
    int numbmodal=this->GetNumberModalities();
    //    for (int m=0; m<numbmodal; m++) {
    //
    //        //Mean Initialisation
    //        ParametersSplit[0]->ValueParameters[m]=ParametersMean[m]-0.5*powf(ParametersVariance[m+numbmodal*m], 0.5);
    //        ParametersSplit[1]->ValueParameters[m]=ParametersMean[m]+0.5*powf(ParametersVariance[m+numbmodal*m], 0.5);
    //
    //        //Variance Initialisation
    //        ParametersSplit[0]->ValueParameters[numbmodal+m+numbmodal*m]=ParametersVariance[m+numbmodal*m]/2.0;
    //        ParametersSplit[1]->ValueParameters[numbmodal+m+numbmodal*m]=ParametersVariance[m+numbmodal*m]/2.0;
    //
    //    }

    // Splitting initialisation with Richardson initialisation
//    float * VarianceDiag=this->GetDiagVarianceDirect();

//    // Possibly to change if we have GetVarianceDirect instead of GetDiagVarianceDirect
//    float * VarianceToStudy=new float[numbmodal*numbmodal];
//    for(int m1=0;m1<numbmodal;m1++){
//        for(int m2=0;m2<numbmodal;m2++){
//            if(m1==m2){
//                VarianceToStudy[m1+m2*numbmodal]=VarianceDiag[m1];
//            }
//            else{
//                VarianceToStudy[m1+m2*numbmodal]=0;
//            }

//        }
//    }

    float * VarianceToStudy=this->GetVarianceDirect();
    SVD SVDForSplit=SVD(VarianceToStudy, numbmodal,numbmodal);
    float * SVDEigen=SVDForSplit.getSingularValues();
    float * SVDVec=SVDForSplit.getV();
    int IndMaxEigen=0;
    float MaxEigen=1E-6;
    float * MaxVec=SVDVec;
    int * Index=new int[numbmodal];
    for(int m=0;m<numbmodal;m++){
        Index[m]=m;
    }
    this->SortingSVDEigen(SVDEigen,Index);
    if(choiceInit<=numbmodal){
        MaxEigen=SVDEigen[choiceInit-1];
        MaxVec=&SVDVec[Index[choiceInit-1]*numbmodal];

    }

    else {

    for (int i=0; i<numbmodal; i++) {
        if (SVDEigen[i]>MaxEigen) {
            IndMaxEigen=i;
            MaxEigen=SVDEigen[i];
            MaxVec=&SVDVec[i*numbmodal];
        }
    }
    }
    delete [] Index;
    Index=NULL;
    delete [] VarianceToStudy;
    VarianceToStudy=NULL;

    float * MeanInitSplit=this->MeanForSplitInitialisationUniform(MaxEigen,MaxVec);
    float * VarianceInitSplit=this->VarianceForSplitInitialisationUniform(MaxEigen,MaxVec);

    // Filling in the parameters
    for (int m1=0; m1<numbmodal; m1++) {
        ParametersSplit[0]->ValueParameters[m1]=MeanInitSplit[m1];
        ParametersSplit[1]->ValueParameters[m1]=MeanInitSplit[m1+numbmodal];
        for (int m2=0; m2<numbmodal; m2++) {
            ParametersSplit[0]->ValueParameters[numbmodal+m1+m2*numbmodal]=VarianceInitSplit[m1+m2*numbmodal];
            ParametersSplit[1]->ValueParameters[numbmodal+m1+m2*numbmodal]=VarianceInitSplit[numbmodal*numbmodal+m1+m2*numbmodal];
        }
    }
    // Memory clearing
    delete [] MeanInitSplit;
    MeanInitSplit=NULL;
    delete [] VarianceInitSplit;
    VarianceInitSplit=NULL;

    return ParametersSplit;
}

// Returns the pointer to the array of Parameters pointer initialising the classes formed by the split
Parameters ** TreeEM::ParametersForSplitting(int choiceInit){

    // Check if it is really a leaf to split
    if (this->GetNumberChildren()!=0) {
        return NULL;
    }
    Parameters ** ParametersSplit=new Parameters* [2];
    ParametersSplit[0]=new Parameters();
    ParametersSplit[1]=new Parameters();
    ParametersSplit[0]->DistributionType=this->GetDistributionType();
    ParametersSplit[1]->DistributionType=this->GetDistributionType();
    ParametersSplit[0]->SizeParameters=this->GetSizeParameters();
    ParametersSplit[1]->SizeParameters=this->GetSizeParameters();
    ParametersSplit[0]->ValueParameters=new float[ParametersSplit[0]->SizeParameters];//{0};
    for(int i=0;i<ParametersSplit[0]->SizeParameters;i++){
        ParametersSplit[0]->ValueParameters[i]=0;
    }
    ParametersSplit[1]->ValueParameters=new float[ParametersSplit[1]->SizeParameters];//{0};
    for(int i=0;i<ParametersSplit[1]->SizeParameters;i++){
        ParametersSplit[1]->ValueParameters[i]=0;
    }
    //    float * ParametersMean=this->GetMean();
    //    float * ParametersVariance=this->GetVariance();

    // Temporary values for splitting initialisation :
    int numbmodal=this->GetNumberModalities();
    //    for (int m=0; m<numbmodal; m++) {
    //
    //        //Mean Initialisation
    //        ParametersSplit[0]->ValueParameters[m]=ParametersMean[m]-0.5*powf(ParametersVariance[m+numbmodal*m], 0.5);
    //        ParametersSplit[1]->ValueParameters[m]=ParametersMean[m]+0.5*powf(ParametersVariance[m+numbmodal*m], 0.5);
    //
    //        //Variance Initialisation
    //        ParametersSplit[0]->ValueParameters[numbmodal+m+numbmodal*m]=ParametersVariance[m+numbmodal*m]/2.0;
    //        ParametersSplit[1]->ValueParameters[numbmodal+m+numbmodal*m]=ParametersVariance[m+numbmodal*m]/2.0;
    //
    //    }

    // Splitting initialisation with Richardson initialisation
    SVD SVDForSplit=SVD(this->GetVariance(), numbmodal,numbmodal);
    float * SVDEigen=SVDForSplit.getSingularValues();
    float * SVDVec=SVDForSplit.getV();
    int IndMaxEigen=0;
    float MaxEigen=1E-6;
    float * MaxVec=SVDVec;
    int * Index=new int[numbmodal];
    for(int m=0;m<numbmodal;m++){
        Index[m]=m;
    }
    this->SortingSVDEigen(SVDEigen,Index);
    if(choiceInit<=numbmodal){
        MaxEigen=SVDEigen[choiceInit-1];
        MaxVec=&SVDVec[Index[choiceInit-1]*numbmodal];

    }

    else {

    for (int i=0; i<numbmodal; i++) {
        if (SVDEigen[i]>MaxEigen) {
            IndMaxEigen=i;
            MaxEigen=SVDEigen[i];
            MaxVec=&SVDVec[i*numbmodal];
        }
    }
    }
    delete [] Index;
    Index=NULL;

    float * MeanInitSplit=this->MeanForSplitInitialisation(MaxEigen,MaxVec);
    float * VarianceInitSplit=this->VarianceForSplitInitialisation(MaxEigen,MaxVec);

    // Filling in the parameters
    for (int m1=0; m1<numbmodal; m1++) {
        ParametersSplit[0]->ValueParameters[m1]=MeanInitSplit[m1];
        ParametersSplit[1]->ValueParameters[m1]=MeanInitSplit[m1+numbmodal];
        for (int m2=0; m2<numbmodal; m2++) {
            ParametersSplit[0]->ValueParameters[numbmodal+m1+m2*numbmodal]=VarianceInitSplit[m1+m2*numbmodal];
            ParametersSplit[1]->ValueParameters[numbmodal+m1+m2*numbmodal]=VarianceInitSplit[numbmodal*numbmodal+m1+m2*numbmodal];
        }
    }
    // Memory clearing
    delete [] MeanInitSplit;
    MeanInitSplit=NULL;
    delete [] VarianceInitSplit;
    VarianceInitSplit=NULL;

    return ParametersSplit;
}

// Applied on the leaf to split returns a pointer to the float array containing the initialisation to use for the mean
float * TreeEM::MeanForSplitInitialisation(float MaxEigen, float * MaxVec){
    if (this->GetNumberChildren()!=0) {
        cout<<"Nothing to split since not leaf"<<endl;
        return NULL;
    }
    int numbmodal=this->GetNumberModalities();
    float * MeanSplitInit=new float[2*numbmodal];
    float * MeanToSplit=this->GetMean();

    for (int m=0; m<numbmodal; m++) {
        MeanSplitInit[m]=MeanToSplit[m]-0.5*MaxVec[m]*sqrtf(MaxEigen);
        MeanSplitInit[m+numbmodal]= MeanToSplit[m]+0.5*MaxVec[m]*sqrtf(MaxEigen);
    }
    //    delete [] SVDForSplit.getU();
    //    delete [] SVDForSplit.getV();
    //    delete [] SVDForSplit.getSingularValues();
    return MeanSplitInit;
}

float * TreeEM::MeanForSplitInitialisationUniform(float MaxEigen,float*MaxVec){
    if (this->GetNumberChildren()!=0 || this->GetDistributionType()!=2) {
        cout<<"Nothing to split since not leaf"<<endl;
        return NULL;
    }
    int numbmodal=this->GetNumberModalities();
    float * MeanSplitInit=new float[2*numbmodal];
    float * MeanToSplit=this->GetMeanDirect();

    for (int m=0; m<numbmodal; m++) {
        MeanSplitInit[m]=MeanToSplit[m]-0.5*MaxVec[m]*sqrtf(MaxEigen);
        MeanSplitInit[m+numbmodal]= MeanToSplit[m]+0.5*MaxVec[m]*sqrtf(MaxEigen);
    }
    //    delete [] SVDForSplit.getU();
    //    delete [] SVDForSplit.getV();
    //    delete [] SVDForSplit.getSingularValues();
    delete[]MeanToSplit;
    MeanToSplit=NULL;
    return MeanSplitInit;
}

void TreeEM::SortingSVDEigen(float * SVDEigen, int * Index){
    int numbmodal=this->GetNumberModalities();
    bool flag_swap=1;
    while(flag_swap){
        flag_swap=0;
        for(int i=0;i<numbmodal-1;i++){
            if(SVDEigen[i+1]>SVDEigen[i]){
                float tmp=SVDEigen[i];
                SVDEigen[i]=SVDEigen[i+1];
                SVDEigen[i+1]=tmp;

                int tmp_ind=Index[i];
                Index[i]=Index[i+1];
                Index[i+1]=tmp_ind;

                flag_swap=1;
            }
        }
    }
}

float * TreeEM::VarianceForSplitInitialisation(float MaxEigen,float * MaxVec){
    if (this->GetNumberChildren()!=0) {
        cout<<"Nothing to split since not leaf"<<endl;
        return NULL;
    }
    int numbmodal=this->GetNumberModalities();
    float * VarianceSplitInit=new float[2*numbmodal*numbmodal];
    float * VarianceToSplit=this->GetVariance();

    for (int m1=0; m1<numbmodal; m1++) {
        for (int m2=0; m2<numbmodal; m2++) {
            VarianceSplitInit[m1+m2*numbmodal]=VarianceToSplit[m1+m2*numbmodal]+2*(0.5-1-0.5*0.5*0.5)*MaxVec[m1]*MaxVec[m2]*MaxEigen+MaxVec[m1]*MaxVec[m2]*MaxEigen;
            VarianceSplitInit[m1+m2*numbmodal+numbmodal*numbmodal]=VarianceToSplit[m1+m2*numbmodal]+2*(0.5*0.5*0.5-0.5-0.5*0.5)*MaxVec[m1]*MaxVec[m2]*MaxEigen+MaxVec[m1]*MaxVec[m2]*MaxEigen;
        }
    }
    //delete &SVDForSplit;
    //    delete [] SVDForSplit.getU();
    //    delete [] SVDForSplit.getV();
    //    delete [] SVDForSplit.getSingularValues();
    return VarianceSplitInit;
}

float * TreeEM::VarianceForSplitInitialisationUniform(float MaxEigen,float * MaxVec){
    if (this->GetNumberChildren()!=0) {
        cout<<"Nothing to split since not leaf"<<endl;
        return NULL;
    }
    int numbmodal=this->GetNumberModalities();
    float * VarianceSplitInit=new float[2*numbmodal*numbmodal];
    float * VarianceToSplit=this->GetVarianceDirect();

    for (int m1=0; m1<numbmodal; m1++) {
        for (int m2=0; m2<numbmodal; m2++) {
            VarianceSplitInit[m1+m2*numbmodal]=VarianceToSplit[m1+m2*numbmodal]+2*(0.5-1-0.5*0.5*0.5)*MaxVec[m1]*MaxVec[m2]*MaxEigen+MaxVec[m1]*MaxVec[m2]*MaxEigen;
            VarianceSplitInit[m1+m2*numbmodal+numbmodal*numbmodal]=VarianceToSplit[m1+m2*numbmodal]+2*(0.5*0.5*0.5-0.5-0.5*0.5)*MaxVec[m1]*MaxVec[m2]*MaxEigen+MaxVec[m1]*MaxVec[m2]*MaxEigen;
        }
    }
    //delete &SVDForSplit;
    //    delete [] SVDForSplit.getU();
    //    delete [] SVDForSplit.getV();
    //    delete [] SVDForSplit.getSingularValues();
    delete [] VarianceToSplit;
    VarianceToSplit=NULL;
    return VarianceSplitInit;
}




// When testing a S operation returns the best of the two models obtained.
TreeEM* TreeEM::RunSplitOperation(vector<SplitKLD *> SplitTry, bool & AcceptanceDecision, SEG_PARAMETERS * segment_param){
    //    SplitKLD * SplitTry=this->GetToSplit();
    int numbclassesSplit=SplitTry.size();
    if (numbclassesSplit==0) {
        return this;
    }
    // Updating of the checking

    TreeEM * GeneralClassToSplit=SplitTry[0]->Parent->GetChild(SplitTry[0]->ChildToSplit)->FindMainNode();
    int IndGenClassSplit=-1;
    if(GeneralClassToSplit->GetParent()!=NULL){
     IndGenClassSplit=GeneralClassToSplit->GetParent()->FindIndex(GeneralClassToSplit);
    }
    cout<< "Trying split operation on node "<<IndGenClassSplit<<endl;
//    SplitTry->Parent->GetSplitCheck()[SplitTry->ChildToSplit]=1;
    // Copy of the current tree
    TreeEM * CopiedTree=this->CopyTree(NULL);
    for(int c=0;c<numbclassesSplit;c++){
//        (SplitTry[c]->Parent)->SplitOperation(SplitTry[c]->ChildToSplit,segment_param->choiceInitSplit,segment_param);
        this->SplitOperation(SplitTry[c],segment_param->choiceInitSplit,segment_param);
    }

    this->ClearSMChecks();

    // Partial EM first
    float PartCLL=0;
    float PartOldCLL=0;
    int PartIteration=0;
    //    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->UpdateDistribution();
    //    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->UpdateNonNormResp();
    //    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->UpdateNormResp();
    //    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->UpdateNormRespRoot();
    //    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->UpdateNonNormWeights();
    //    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->UpdateNormWeights();
    //    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->RunFullEM(PartCLL, PartOldCLL, PartIteration);
    // Run EM
    float CompleteLogLikelihood=0;
    float OldCompleteLogLikelihood=0;
    int Iteration=0;
    //    this->UpdateDistribution();
    this->UpdateNonNormResp(segment_param);
    this->UpdateNormResp();
//    this->SaveAllClasses("/Users/Carole/Documents/PhD/TestSplitResult.nii.gz",segment_param);
    this->UpdateNormRespRoot();
    if(segment_param->flag_MRF){
        this->ClearMRFModel();
        float * GMatrixToSet=NULL;
        if(segment_param->flag_GMatrixIn){
            GMatrixToSet=this->PrepareGMatrixFromFile(segment_param->filename_GMatrix);
            if(GMatrixToSet==NULL){
                cout<< "Pb in the preparation of GMatrix"<<endl;
                segment_param->flag_GMatrix=0;
            }
        }
        else{
            GMatrixToSet=this->MRFOptSolveLS();
        }
        this->SetGMatrix(GMatrixToSet);
        if(GMatrixToSet!=NULL){
            delete [] GMatrixToSet;
            GMatrixToSet=NULL;
        }
        this->UpdateMRF();
    }
//    this->UpdateNonNormWeights();
//    this->UpdateNormWeights();

    //this->RunFullEM(PartCLL, PartOldCLL, PartIteration);
    if(segment_param->flag_CEM){
        this->RunFullCEM(CompleteLogLikelihood,OldCompleteLogLikelihood,Iteration,segment_param);
    }
    else{
    this->RunFullEM(CompleteLogLikelihood, OldCompleteLogLikelihood,Iteration,segment_param);
    }
//    this->SaveAllClasses("/Users/Carole/Documents/PhD/TestSplitResult");
    // Test of the operation
    float OldCriterion=0;
    float NewCriterion=0;
    if(this->GetDPChildrenDirect()!=NULL){
        if(segment_param->flag_Countmod){
            OldCriterion=CopiedTree->CriterionCalculationSplitDP();
            NewCriterion=this->CriterionCalculationSplitDP();
        }
        else{
            if(segment_param->flag_DistClassInd){
         OldCriterion=CopiedTree->CriterionCalculationDPInd();
         NewCriterion=this->CriterionCalculationDPInd();
            }
            else{
//                OldCriterion=CopiedTree->CriterionCalculationDPNonInd(SplitTry[c]->Parent->GetChild(SplitTry->ChildToSplit)->FindGeneralClass());
//                NewCriterion=this->CriterionCalculationDPNonInd(SplitTry[c]->Parent->GetChild(SplitTry->ChildToSplit)->FindGeneralClass());
            }
        }
    }
    else{
     OldCriterion=CopiedTree->CriterionCalculation();
     NewCriterion=this->CriterionCalculation();
    }
    cout<< "The new criterion is "<<NewCriterion<<" and the old criterion is "<<OldCriterion<<endl;

//    // Clear SplitTry
//    if (SplitTry!=NULL) {
//        delete SplitTry;
//        SplitTry=NULL;
//    }

    bool TestNormResp=this->AreNormRespValid();
    bool TestNormWeight=this->AreWeightsValid();
    if(!TestNormResp){
        cout<<"Pb of norm resp in Split operation"<<endl;
    }
    if(!TestNormWeight){
        cout<<"Pb of norm weights in Split operation"<<endl;
    }

    if (OldCriterion>=NewCriterion) { // Criterion of change not met, delete the modified tree and return the copied one
        delete this;
        AcceptanceDecision=0;
        cout<<"Split not accepted"<<endl;
        return CopiedTree;
    }
    else{ // Criterion of change met delete the not modified tree stored in CopiedTree and return the new model.
        delete CopiedTree;
        AcceptanceDecision=1;
        cout<<"Split accepted"<<endl;
//        int numbchild=this->GetNumberChildren();
//        for(int c=0;c<numbchild;c++){
//            this->GetChild(c)->PutAllLeavesToChildrenLevel();
//        }
        this->PutAllLeavesToMainNodesChildrenLevel();
        TestNormWeight=this->AreWeightsValid();
        if(!TestNormWeight){
            cout<<"Pb when putting all leaves to main nodes children level at split operation"<<endl;
        }
        this->ClearSMChecks();
        return this;
    }
}


// When testing a S operation returns the best of the two models obtained.
TreeEM* TreeEM::RunSplitOperation(SplitKLD * SplitTry, bool & AcceptanceDecision, SEG_PARAMETERS * segment_param){
    //    SplitKLD * SplitTry=this->GetToSplit();
    if (SplitTry==NULL) {
        return this;
    }
    // Updating of the checking
    cout<< "Trying split operation"<<endl;
    SplitTry->Parent->GetSplitCheck()[SplitTry->ChildToSplit]=1;
    // Copy of the current tree
    TreeEM * CopiedTree=this->CopyTree(NULL);
    (SplitTry->Parent)->SplitOperation(SplitTry->ChildToSplit,segment_param->choiceInitSplit,segment_param);
    this->ClearSMChecks();

    // Partial EM first
    float PartCLL=0;
    float PartOldCLL=0;
    int PartIteration=0;
    //    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->UpdateDistribution();
    //    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->UpdateNonNormResp();
    //    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->UpdateNormResp();
    //    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->UpdateNormRespRoot();
    //    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->UpdateNonNormWeights();
    //    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->UpdateNormWeights();
    //    (SplitTry->Parent->GetChild(SplitTry->ChildToSplit))->RunFullEM(PartCLL, PartOldCLL, PartIteration);
    // Run EM
    float CompleteLogLikelihood=0;
    float OldCompleteLogLikelihood=0;
    int Iteration=0;
    //    this->UpdateDistribution();
    this->UpdateNonNormResp(segment_param);
    this->UpdateNormResp();
//    this->SaveAllClasses("/Users/Carole/Documents/PhD/TestSplitResult");
    this->UpdateNormRespRoot();
    if(segment_param->flag_MRF){
        this->ClearMRFModel();
        float * GMatrixToSet=NULL;
        if(segment_param->flag_GMatrixIn){

          GMatrixToSet=this->PrepareGMatrixFromFile(segment_param->filename_GMatrix);
        if(GMatrixToSet==NULL){
            cout<< "Pb in the preparation of GMatrix"<<endl;
            segment_param->flag_GMatrix=0;
        }
        }
        else{
            GMatrixToSet=this->MRFOptSolveLS();
        }
        this->SetGMatrix(GMatrixToSet);
        if(GMatrixToSet!=NULL){
            delete [] GMatrixToSet;
            GMatrixToSet=NULL;
        }
        this->UpdateMRF();
    }
//    this->UpdateNonNormWeights();
//    this->UpdateNormWeights();

    //this->RunFullEM(PartCLL, PartOldCLL, PartIteration);
    if(segment_param->flag_CEM){
        this->RunFullCEM(CompleteLogLikelihood, OldCompleteLogLikelihood,Iteration,segment_param);
    }
    else{
    this->RunFullEM(CompleteLogLikelihood, OldCompleteLogLikelihood,Iteration,segment_param);
    }
//    this->SaveAllClasses("/Users/Carole/Documents/PhD/TestSplitResult");
    // Test of the operation
    float OldCriterion=0;
    float NewCriterion=0;
    if(this->GetDPChildrenDirect()!=NULL){
        if(segment_param->flag_Countmod){
            OldCriterion=CopiedTree->CriterionCalculationSplitDP();
            NewCriterion=this->CriterionCalculationSplitDP();
        }
        else{
            if(segment_param->flag_DistClassInd){
         OldCriterion=CopiedTree->CriterionCalculationDPInd();
         NewCriterion=this->CriterionCalculationDPInd();
            }
            else{
                OldCriterion=CopiedTree->CriterionCalculationDPNonInd(SplitTry->Parent->GetChild(SplitTry->ChildToSplit)->FindGeneralClass());
                NewCriterion=this->CriterionCalculationDPNonInd(SplitTry->Parent->GetChild(SplitTry->ChildToSplit)->FindGeneralClass());
            }
        }
    }
    else{
     OldCriterion=CopiedTree->CriterionCalculation();
     NewCriterion=this->CriterionCalculation();
    }
    cout<< "The new criterion is "<<NewCriterion<<" and the old criterion is "<<OldCriterion<<endl;

    // Clear SplitTry
    if (SplitTry!=NULL) {
        delete SplitTry;
        SplitTry=NULL;
    }

    if (OldCriterion>=NewCriterion) { // Criterion of change not met, delete the modified tree and return the copied one
        delete this;
        AcceptanceDecision=0;
        cout<<"Split not accepted"<<endl;
        return CopiedTree;
    }
    else{ // Criterion of change met delete the not modified tree stored in CopiedTree and return the new model.
        delete CopiedTree;
        AcceptanceDecision=1;
        cout<<"Split accepted"<<endl;
//        int numbchild=this->GetNumberChildren();
//        for(int c=0;c<numbchild;c++){
//            this->GetChild(c)->PutAllLeavesToChildrenLevel();
//        }
        this->PutAllLeavesToMainNodesChildrenLevel();
        this->ClearSMChecks();
        return this;
    }
}


/*************** METHODS FOR ATLAS ADAPTATION *******************/
float* TreeEM::AdaptPriors(SEG_PARAMETERS * segment_param){
    // Check if there are some priors to adapt
    if(this->GetPriors()==NULL){
        cout<<"No need for adaptation of statistical atlas since no atlas there"<<endl;
        return this->GetPriorsAdapted();
    }
    int numel=this->GetNumberElements();

    // Obtention of the needed dimensions to then blurr the responsabilities

    vector<int> DimImage;
    DimImage.push_back(this->GetPriors()->nx);
    DimImage.push_back(this->GetPriors()->ny);
    DimImage.push_back(this->GetPriors()->nz);

    // Temperation by the atlas with weight segment_param->AtlasWeight
    float * Priors_PTR=static_cast<float *>(this->GetPriors()->data);
    float * SmoothedPriors=new float[numel];
    float Weight=segment_param->AtlasWeight;
    Weight=Weight<=0?0:Weight;
    Weight=Weight>=1?1:Weight;
    float NormalisedWeight=this->GetPartNormWeightAbovePriors();
    for(int i=0;i<numel;i++,Priors_PTR++){
        SmoothedPriors[i]=NormalisedWeight*(1-Weight)*(*Priors_PTR);
    }

    // Smoothing of the responsabilities (function developed in DirichletPriors) if there is some to do
    if(segment_param->AtlasWeight<=0){
        return SmoothedPriors;
    }
    else {
        // Initialise float * such that contains values of NormResp on active voxels
        float * NormResp_PTR=this->GetNormResp();
        float * ImageResp=new float[numel];
        int * L2S_PTR=this->GetL2S();
        for(int i=0;i<numel;i++,L2S_PTR++){
            if(*L2S_PTR>=0){
                ImageResp[i]=*NormResp_PTR;
                NormResp_PTR++;
            }
            else{
                ImageResp[i]=0;
            }
        }
    float * SmoothedPriors_tmp=GaussianBlurring(ImageResp, segment_param->AtlasSmoothing, DimImage);
//    SaveTmpResult(SmoothedPriors_tmp,"/Users/Carole/Documents/PhD/SmoothResult.nii.gz");
    delete [] ImageResp;
    ImageResp=NULL;
    DimImage.clear();

    for(int i=0;i<numel;i++,Priors_PTR++){
        SmoothedPriors[i]+=Weight*SmoothedPriors_tmp[i];
    }
    Priors_PTR=static_cast<float *>(this->GetPriors()->data);
    delete [] SmoothedPriors_tmp;
    SmoothedPriors_tmp=NULL;

//    nifti_image * Result=nifti_copy_nim_info(this->GetDataImage());
//    Result->data = (void *) calloc(Result->nvox, sizeof(float));
//    float * Result_PTRtmp=static_cast<float *>(Result->data);
//    float * SmoothedPriors_PTR=SmoothedPriors;
//    for (int i=0; i<numel; i++, Result_PTRtmp++,SmoothedPriors_PTR++) {
//        *Result_PTRtmp=*SmoothedPriors_PTR;
//    }
//    nifti_set_filenames(Result, "/Users/Carole/Documents/PhD/TestAtlasSmooth.nii.gz", 0, 0);
//    nifti_image_write(Result);
//    nifti_image_free(Result);
    }

    return SmoothedPriors;
}

nifti_image * TreeEM::BuildConstantPriors(float WeightInput){
    nifti_image * ConstantPriors=nifti_copy_nim_info(this->GetDataImage());
    ConstantPriors->dim[0]=3;
    ConstantPriors->dim[4]=1;
    int numel=this->GetNumberElements();
    nifti_update_dims_from_array(ConstantPriors);
    ConstantPriors->data = (void *) calloc(ConstantPriors->nvox, sizeof(float));
    float * ConstantPriors_PTRtmp=static_cast<float *>(ConstantPriors->data);
    for(int i=0;i<numel;i++,ConstantPriors_PTRtmp++){
        *ConstantPriors_PTRtmp=WeightInput;
    }
    return ConstantPriors;
}

nifti_image * TreeEM::CreatePriorsFromAdaptedPriors(){
    float * PriorsAdaptedToUse=this->GetPriorsAdapted();
    if(PriorsAdaptedToUse==NULL){
        return NULL;
    }
    nifti_image * PriorsCreated=nifti_copy_nim_info(this->GetPriors());
    PriorsCreated->data=(void *) calloc(PriorsCreated->nvox,sizeof(float));
    float * PriorsCreated_PTR=static_cast<float *>(PriorsCreated->data);
    int numel=PriorsCreated->nvox;
    for(int i=0;i<numel;i++,PriorsCreated_PTR++,PriorsAdaptedToUse++){
        *PriorsCreated_PTR=*PriorsAdaptedToUse;
    }
    return PriorsCreated;
}

nifti_image * TreeEM::TransformNormRespIntoPriors(SEG_PARAMETERS * segment_param){
    //Create image in which Priors will be stored.
    nifti_image * PriorsResult = nifti_copy_nim_info(this->GetDataImage());
    //Make in only one time point 3D image
    PriorsResult->dim[0]=3;
    PriorsResult->dim[4]=1;
    PriorsResult->dim[5]=1;
    nifti_update_dims_from_array(PriorsResult);
    PriorsResult->data=(void *)calloc(PriorsResult->nvox,sizeof(float));
    float * PriorsResult_PTR=static_cast<float*>(PriorsResult->data);
    // Transform NormResp into float pointer of size of image
    int numel = this->GetNumberElements();
    float * NormRespImage=new float[numel];
    // Initialisation by 0
    for(int i=0;i<numel;i++){
        NormRespImage[i]=0;
    }
    float * NormResp_PTR=this->GetNormResp();
    int * L2S_PTR=this->GetL2S();
    for(int i=0;i<numel;i++,L2S_PTR++){
        if(*L2S_PTR>=0){
            NormRespImage[i]=*NormResp_PTR;
            *NormResp_PTR++;
        }
    }
    vector<int> Dim;
    for(int d=1;d<=3;d++){
        Dim.push_back(PriorsResult->dim[d]);
    }
    float * NormRespImageSmoothed=GaussianBlurring(NormRespImage,segment_param->AtlasSmoothing,Dim);
    delete [] NormRespImage;
    NormRespImage=NULL;
    // Copy result into PriorsResultData
    for(int i=0;i<numel;i++,PriorsResult_PTR++){
        *PriorsResult_PTR=NormRespImageSmoothed[i];
    }
    delete [] NormRespImageSmoothed;
    NormRespImageSmoothed=NULL;
    return PriorsResult;
}


/***************** METHODS FOR BIAS FIELD CORRECTION *************/

// Returns a float array with all the basis functions to be used in the bias field correction
/** REMARK : To be thought again to see if we want to get it at each iteration and delete it afterwards or store it at the root as L2S or S2L **/
float * TreeEM::MakeBasisFunctions(){
    int numbFunctions=(int)((BForder+1)*(BForder+2)*(BForder+3))/6;
    int numelmasked=this->GetNumberMaskedElements();
    float * BasisFunctions=new float[numelmasked*numbFunctions];//{0};
    for(int i=0;i<numelmasked*numbFunctions;i++){
        BasisFunctions[i]=0;
    }
    // Find the factors for the affine transformation of the indices so that the power functions are applied between 0 and 1
    float FactorX=(float)(2.0/(this->GetDataImage()->nx-1));
    float FactorY=(float)(2.0/(this->GetDataImage()->ny-1));
    float FactorZ=(float)(2.0/(this->GetDataImage()->nz-1));

    int SizePlane=this->GetDataImage()->nx*this->GetDataImage()->ny;
    int SizeColumn=this->GetDataImage()->nx;
    int zValue=0;
    int yValue=0;
    int xValue=0;

    int l=0;//Represent the index of the considered basis function
    for (int order=0; order<=BForder; order++) {
        for (int xorder=0; xorder<=order; xorder++) {
            for (int yorder=0; yorder<=order-xorder; yorder++) {
                int zorder=order-xorder-yorder; // the sum of the orders for the different directions must be at most BForder.
                // Reinitialisation of the pointers
                int * S2L_PTR=this->GetS2L();
                float * BasisFunctions_PTR=&BasisFunctions[numelmasked*l];
                for (int i=0; i<numelmasked; i++,S2L_PTR++,BasisFunctions_PTR++) {
                    // conversion of S2L value into x,y and z component
                    zValue=*S2L_PTR/SizePlane;
                    yValue=(*S2L_PTR-zValue*SizePlane)/SizeColumn;
                    xValue=(*S2L_PTR)-zValue*SizePlane-yValue*SizeColumn;

                    // Filling the BasisFunctions matrix with the proper values modified to be symmetric and between -1 and 1.
                    *BasisFunctions_PTR=pow_int(xValue*FactorX-1, xorder)*pow_int(yValue*FactorY-1, yorder)*pow_int(zValue*FactorZ-1, zorder);
                }
                l++;
            }
        }
    }
    return BasisFunctions;
}

// Returns the vector of numbmodal pointers to float array containing the values for the W matrix used in VL BF correction.
float * TreeEM::MakeWMatrixChildren(){
    int numbmodal=this->GetNumberModalities();
    int numelmasked=this->GetNumberMaskedElements();
    int numbchild=this->GetNumberChildren();

    float * WMatrix=new float[numelmasked*numbmodal];//{0};
    if(WMatrix==NULL){
        cout<<"Did not allocate Wmatrix properly"<<endl;
    }
    for(int i=0;i<numelmasked*numbmodal;i++){
        WMatrix[i]=0;
    }
    for(int m=0;m<numbmodal;m++){
        for (int c=0; c<numbchild; c++) {
            //Reinitialisation of the pointers
            float * WMatrix_PTR=&WMatrix[m*numelmasked];
            float * NormResp_PTR=this->GetChild(c)->GetNormResp();
            float * DiagVarianceValue=this->GetChild(c)->GetDiagVarianceDirect();
            float invVariance=DiagVarianceValue[m]>1E-6?1.0/DiagVarianceValue[m]:1E6;
            for (int i=0; i<numelmasked; i++,WMatrix_PTR++,NormResp_PTR++) {
                // Filling of the matrix as sum of normresp for the classes divided by the value of the diagonal in the covariance matrix.
                *WMatrix_PTR+=(float)*NormResp_PTR*invVariance;
            }
            delete [] DiagVarianceValue;
        }
    }

    for(int m=0;m<numbmodal;m++){
       int  CountWMatrixZero=0;
        for(int i=0;i<numelmasked;i++){
            if(WMatrix[i+m*numelmasked]==0){
                CountWMatrixZero++;
            }
        }
        if(CountWMatrixZero>0){
            cout<<"Pb with sum NormResp in calculation WMatrix for modality "<<m<<" is "<<CountWMatrixZero<<endl;
        }
    }
    return WMatrix;
}

// Returns in a vector of pointers  the RMatrix needed for VL BF correction according to the different modalities.
float * TreeEM::MakeRMatrixChildren(){


    int numbmodal=this->GetNumberModalities();
    int numelmasked=this->GetNumberMaskedElements();
    float * RMatrix=new float[numelmasked*numbmodal];
    int numbchild=this->GetNumberChildren();
    int numel=this->GetNumberElements();
    float * WMatrix=this->MakeWMatrixChildren();
    float * Data=static_cast<float*>(this->GetDataImage()->data);
    float * Data_PTR=Data;
    int * L2S_PTR=this->GetL2S();
    for(int i=0;i<numelmasked*numbmodal;i++){
        RMatrix[i]=0;
    }
    for(int m=0;m<numbmodal;m++){
        L2S_PTR=this->GetL2S();
        Data_PTR=&Data[m*numel];
        float * RMatrix_PTR=&RMatrix[m*numelmasked];
        for (int c=0; c<numbchild; c++) {
            //Reinitialisation of the pointers
            RMatrix_PTR=&RMatrix[m*numelmasked];
            float * NormResp_PTR=this->GetChild(c)->GetNormResp();
            float * WMatrix_PTR=&WMatrix[m*numelmasked];
            float *  Mean=this->GetChild(c)->GetMeanDirect();
            float * DiagVarianceValue=this->GetChild(c)->GetDiagVarianceDirect();
            float invVariance=DiagVarianceValue[m]>1E-6?1.0/DiagVarianceValue[m]:1E6;
            for (int i=0; i<numelmasked; i++,WMatrix_PTR++,NormResp_PTR++,RMatrix_PTR++) {
                // Filling of the matrix as normalised weighted sum of the means
                *RMatrix_PTR+=(float)(*NormResp_PTR*Mean[m]*invVariance)/(*WMatrix_PTR);
            }
            delete [] Mean;
            delete [] DiagVarianceValue;
        }
        RMatrix_PTR=&RMatrix[m*numelmasked];
        for (int i=0; i<numel; i++, Data_PTR++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
                *RMatrix_PTR=*Data_PTR-*RMatrix_PTR;
                RMatrix_PTR++;
            }
        }
    }
    if(WMatrix!=NULL){
        delete [] WMatrix;
        WMatrix=NULL;
    }
    return RMatrix;
}

// Returns a vector of pointers to the AtWA matrices used for the BF correction. The obtained matrix will need to be inverted afterwards.
float * TreeEM::MakeAtWAMatrixChildren(){
    float * WMatrix=this->MakeWMatrixChildren();
    //    float * BFFunctions=this->MakeBasisFunctions();
    int numbmodal=this->GetNumberModalities();
    int numelmasked=this->GetNumberMaskedElements();
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;

    float FactorX=(float)(2.0/(this->GetDataImage()->nx-1));
    float FactorY=(float)(2.0/(this->GetDataImage()->ny-1));
    float FactorZ=(float)(2.0/(this->GetDataImage()->nz-1));

    int SizePlane=this->GetDataImage()->nx*this->GetDataImage()->ny;
    int SizeColumn=this->GetDataImage()->nx;

    int XorderBF[MaxNumbBF];
    int YorderBF[MaxNumbBF];
    int ZorderBF[MaxNumbBF];
    int tmpnumbBF=0;
    for(int order=0;order<=BForder;order++){
        for(int zorder=0;zorder<=order;zorder++){
            for(int yorder=0;yorder<=order-zorder;yorder++){
                int xorder=order-zorder-yorder;
                XorderBF[tmpnumbBF]=xorder;
                YorderBF[tmpnumbBF]=yorder;
                ZorderBF[tmpnumbBF]=zorder;
                tmpnumbBF++;
            }
        }
    }
    float* AtWAMatrixResult=new float[numbBF*numbBF*numbmodal];
    for(int i=0;i<numbBF*numbBF*numbmodal;i++){
        AtWAMatrixResult[i]=0;
    }
    float AtWAMatrix[MaxNumbBF*MaxNumbBF];
    for(int i=0;i<MaxNumbBF*MaxNumbBF;i++){
        AtWAMatrix[i]=0;
    }
    // Not optimal since all calculated whereas symmetric matrix obtained
    for(int m=0;m<numbmodal;m++){

        int * S2L_PTR=this->GetS2L();
        for(int i=0;i<numelmasked;i++,S2L_PTR++){
            int zValue=*S2L_PTR/SizePlane;
            int yValue=(*S2L_PTR-zValue*SizePlane)/SizeColumn;
            int xValue=(*S2L_PTR)-zValue*SizePlane-yValue*SizeColumn;

            float xUsed=xValue*FactorX-1;
            float yUsed=yValue*FactorY-1;
            float zUsed=zValue*FactorZ-1;


            for (int j1=0; j1<numbBF; j1++) {
                float Value1=pow_int(xUsed,XorderBF[j1])*pow_int(yUsed,YorderBF[j1])*pow_int(zUsed,ZorderBF[j1]);
                if(Value1>1||Value1<-1){
                    cout<<"Pb with Value1 "<<endl;
                }
                for (int j2=0; j2<numbBF; j2++) {
                    float Value2=pow_int(xUsed,XorderBF[j2])*pow_int(yUsed,YorderBF[j2])*pow_int(zUsed,ZorderBF[j2]);
                    if(Value2>1||Value2<-1){
                        cout<<"Pb with Value2 "<<endl;
                    }
                    float WMatrix_Val=WMatrix[m*numelmasked+i];
                    AtWAMatrix[j1+numbBF*j2]+=(float)Value1*Value2*(WMatrix_Val);
                }
            }
        }
        for(int l=0;l<numbBF*numbBF;l++){
            AtWAMatrixResult[l+m*numbBF*numbBF]=AtWAMatrix[l];
        }
        for(int l=0;l<MaxNumbBF*MaxNumbBF;l++){
            AtWAMatrix[l]=0;
        }
    }
    //clearing WMatrix
    if(WMatrix!=NULL){
        delete [] WMatrix;
        WMatrix=NULL;
    }
    cout<<"AtWAMatrix done"<<endl;
    return AtWAMatrixResult;
}

// Returns a vector of float pointers to the AtWR vectors needed for the VL BF correction
float * TreeEM::MakeAtWRVectorChildren(){
    // Initialisation and getting all needed presteps with other matrices (A, W and R)

    int numelmasked=this->GetNumberMaskedElements();
    int numbmodal=this->GetNumberModalities();
    //    float * BFFunctions=this->MakeBasisFunctions();
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    float * WMatrix=this->MakeWMatrixChildren();
    float * RMatrix=this->MakeRMatrixChildren();


    float * AtWRVector=new float[numbBF*numbmodal];//{0};
    for(int i=0;i<numbBF*numbmodal;i++){
        AtWRVector[i]=0;
    }

    float FactorX=(float)(2.0/(this->GetDataImage()->nx-1));
    float FactorY=(float)(2.0/(this->GetDataImage()->ny-1));
    float FactorZ=(float)(2.0/(this->GetDataImage()->nz-1));

    int SizePlane=this->GetDataImage()->nx*this->GetDataImage()->ny;
    int SizeColumn=this->GetDataImage()->nx;

    int XorderBF[MaxNumbBF];
    int YorderBF[MaxNumbBF];
    int ZorderBF[MaxNumbBF];
    int tmpnumbBF=0;
    for(int order=0;order<=BForder;order++){
        for(int zorder=0;zorder<=order;zorder++){
            for(int yorder=0;yorder<=order-zorder;yorder++){
                int xorder=order-zorder-yorder;
                XorderBF[tmpnumbBF]=xorder;
                YorderBF[tmpnumbBF]=yorder;
                ZorderBF[tmpnumbBF]=zorder;
                tmpnumbBF++;
            }
        }
    }
    int * S2L_PTR=this->GetS2L();
    for(int i=0;i<numelmasked;i++,S2L_PTR++){


        int zValue=*S2L_PTR/SizePlane;
        int yValue=(*S2L_PTR-zValue*SizePlane)/SizeColumn;
        int xValue=(*S2L_PTR)-zValue*SizePlane-yValue*SizeColumn;

        float xUsed=xValue*FactorX-1;
        float yUsed=yValue*FactorY-1;
        float zUsed=zValue*FactorZ-1;

        for (int j=0; j<numbBF; j++) {
            for (int m=0; m<numbmodal; m++) {
                //            float * BFFunctions_PTR=&BFFunctions[j*numelmasked];
                //            float * WMatrix_PTR=&WMatrix[m*numelmasked];
                //            float * RMatrix_PTR=&RMatrix[m*numelmasked];
                float WMatrix_Val=WMatrix[m*numelmasked+i];
                float RMatrix_Val=RMatrix[m*numelmasked+i];
                float BFValue=pow_int(xUsed,XorderBF[j])*pow_int(yUsed,YorderBF[j])*pow_int(zUsed,ZorderBF[j]);
                AtWRVector[j+m*numbBF]+=(float)BFValue*WMatrix_Val*RMatrix_Val;
            }
        }
    }
    // Clearing memory allocation
    delete [] WMatrix;
    WMatrix=NULL;
    delete [] RMatrix;
    RMatrix=NULL;
    //        delete [] BFFunctions;
    //        BFFunctions=NULL;

    return AtWRVector;
}

// Returns a vector of float pointers to the AtWR vectors needed for the VL BF correction
float * TreeEM::MakeBFCoeffsDirect(){
    // Initialisation and getting all needed presteps with other matrices (A, W and R)

    int numelmasked=this->GetNumberMaskedElements();
    int numbmodal=this->GetNumberModalities();
    //    float * BFFunctions=this->MakeBasisFunctions();
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    float * WMatrix=this->MakeWMatrixChildren();
    float * RMatrix=this->MakeRMatrixChildren();
    float  AtWRResult[MaxNumbModal*MaxNumbBF];
    float * AtWAResult=new float[numbBF*numbBF*numbmodal];
    for(int i=0;i<MaxNumbModal*MaxNumbBF;i++){
        AtWRResult[i]=0;
    }
    for(int i=0;i<numbBF*numbBF*numbmodal;i++){
        AtWAResult[i]=0;
    }

    float FactorX=(float)(2.0/(this->GetDataImage()->nx-1));
    float FactorY=(float)(2.0/(this->GetDataImage()->ny-1));
    float FactorZ=(float)(2.0/(this->GetDataImage()->nz-1));

    int SizePlane=this->GetDataImage()->nx*this->GetDataImage()->ny;
    int SizeColumn=this->GetDataImage()->nx;

    int XorderBF[MaxNumbBF];
    int YorderBF[MaxNumbBF];
    int ZorderBF[MaxNumbBF];
    int tmpnumbBF=0;
    for(int order=0;order<=BForder;order++){
        for(int zorder=0;zorder<=order;zorder++){
            for(int yorder=0;yorder<=order-zorder;yorder++){
                int xorder=order-zorder-yorder;
                XorderBF[tmpnumbBF]=xorder;
                YorderBF[tmpnumbBF]=yorder;
                ZorderBF[tmpnumbBF]=zorder;
                tmpnumbBF++;
            }
        }
    }
    int * S2L_PTR=this->GetS2L();

    int CountNanRMatrix=0;
    int CountNanWMatrix=0;
    for(int i=0;i<numelmasked;i++,S2L_PTR++){


        int zValue=*S2L_PTR/SizePlane;
        int yValue=(*S2L_PTR-zValue*SizePlane)/SizeColumn;
        int xValue=(*S2L_PTR)-zValue*SizePlane-yValue*SizeColumn;

        float xUsed=xValue*FactorX-1;
        float yUsed=yValue*FactorY-1;
        float zUsed=zValue*FactorZ-1;


        for (int j1=0; j1<numbBF; j1++) {
            float BFValue1=pow_int(xUsed,XorderBF[j1])*pow_int(yUsed,YorderBF[j1])*pow_int(zUsed,ZorderBF[j1]);
            for (int m=0; m<numbmodal; m++) {
                float WMatrix_Val=WMatrix[m*numelmasked+i];
                float RMatrix_Val=RMatrix[m*numelmasked+i];
                // Checking non Nan of WMatrix and RMatrix
                if(WMatrix_Val!=WMatrix_Val){
                    CountNanWMatrix++;
                }
                if(RMatrix_Val!=RMatrix_Val){
                    CountNanRMatrix++;
                }
                AtWRResult[j1+m*numbBF]+=(float)BFValue1*WMatrix_Val*RMatrix_Val;
                for (int j2=0; j2<numbBF; j2++) {
                    float BFValue2=pow_int(xUsed,XorderBF[j2])*pow_int(yUsed,YorderBF[j2])*pow_int(zUsed,ZorderBF[j2]);

                    AtWAResult[j1+numbBF*j2+m*numbBF*numbBF]+=(float)BFValue1*BFValue2*(WMatrix_Val);
                }
            }
        }

    }
    if(CountNanWMatrix >0 || CountNanRMatrix >0){
    cout<<"Nan Wmatrix is "<<CountNanWMatrix<<" and Nan RMatrix is "<<CountNanRMatrix<<endl;
    }

//    for(int m=0;m<numbmodal;m++){
//        for(int j1=0;j1<numbBF;j1++){
//            for(int j2=0;j2<numbBF;j2++){
//                cout<< AtWAResult[j1+numbBF*j2+m*numbBF*numbBF]<<"  ";
//            }
//            cout<<endl;
//        }
//        cout<<endl;
//    }



    // Clearing memory allocation
    delete [] WMatrix;
    WMatrix=NULL;
    delete [] RMatrix;
    RMatrix=NULL;

    // Calculate Inverse of AtWAResult
    float * invAtWAMatrix=new float[numbBF*numbBF*numbmodal];//{0};
    for(int i=0;i<numbBF*numbBF*numbmodal;i++){
        invAtWAMatrix[i]=0;
    }
    // Allocation of the matrix
    for (int m=0; m<numbmodal; m++) {
//        matrix<float> AtWAMatrix_mat=matrix<float>(numbBF);
        float * AtWAMatrix_mat=new float[numbBF*numbBF];
        float * AtWAMatrix_tmp=&AtWAResult[m*numbBF*numbBF];
//        bool success=0;
        for (int j1=0; j1<numbBF; j1++) {
            for (int j2=0; j2<numbBF; j2++) {
//                AtWAMatrix_mat.setvalue(j1, j2, AtWAMatrix_tmp[j1+j2*numbBF]);
                AtWAMatrix_mat[j1+j2*numbBF]=AtWAMatrix_tmp[j1+j2*numbBF];
            }
        }
//        AtWAMatrix_mat.invert();
        invertMatrix(AtWAMatrix_mat,numbBF);
        int CountNanBF=0;
        for (int j1=0; j1<numbBF; j1++) {
            for (int j2=0; j2<numbBF; j2++) {
//                AtWAMatrix_mat.getvalue(j1, j2, invAtWAMatrix[j1+j2*numbBF+m*numbBF*numbBF], success);
                invAtWAMatrix[j1+j2*numbBF+m*numbBF*numbBF]=AtWAMatrix_mat[j1+j2*numbBF];
                if(invAtWAMatrix[j1+j2*numbBF+m*numbBF*numbBF]!=invAtWAMatrix[j1+j2*numbBF+m*numbBF*numbBF]){
                    CountNanBF++;
                }
            }
        }
        if(CountNanBF!=0){
            cout<<"number of nans is "<<CountNanBF<<endl;
        }
        delete [] AtWAMatrix_mat;
        AtWAMatrix_mat=NULL;
    }
    // Clearing memory
    delete [] AtWAResult;
    AtWAResult=NULL;

    // Now InvAtWA is obtained and we also have AtWR so we can directly calculate the coeffs
    float * FinalBFCoeffsResult=new float[numbBF*numbmodal];
    for(int l=0;l<numbBF*numbmodal;l++){
        FinalBFCoeffsResult[l]=0;
    }
    for (int m=0; m<numbmodal; m++) {
        float * invAtWA=&invAtWAMatrix[m*numbBF*numbBF];
        float * AtWRV_PTR=&AtWRResult[m*numbBF];
        float * invAtWA_PTR=invAtWA;
        for (int j1=0; j1<numbBF; j1++) {
            invAtWA_PTR=&invAtWA[j1*numbBF];
            AtWRV_PTR=&AtWRResult[m*numbBF];
            float FinalBFCoeffsResult_tmp=0;
            for (int j2=0; j2<numbBF; j2++,AtWRV_PTR++,invAtWA_PTR++) {
                FinalBFCoeffsResult_tmp+=(PrecisionTYPE)(*AtWRV_PTR)*(*invAtWA_PTR);
            }
            FinalBFCoeffsResult[j1+m*numbBF]=(float)FinalBFCoeffsResult_tmp;
        }

    }
    // Clearing memory
    delete [] invAtWAMatrix;
    invAtWAMatrix=NULL;
    //    delete [] AtWRVector;
    //    AtWRVector=NULL;
    return FinalBFCoeffsResult;

}

// Returns in a vector of pointers the inverted matrix AtWA for each modality.
float* TreeEM::MakeInvAtWAMatrixChildren(){
    float * AtWAMatrix=this->MakeAtWAMatrixChildren();
    int numbBF=(int)((BForder+1)*(BForder+2)*(BForder+3))/6;
    int numbmodal=this->GetNumberModalities();
    // float * invAtWAMatrixResult;
    float * invAtWAMatrix=new float[numbBF*numbBF*numbmodal];//{0};
    for(int i=0;i<numbBF*numbBF*numbmodal;i++){
        invAtWAMatrix[i]=0;
    }
    // Allocation of the matrix
    for (int m=0; m<numbmodal; m++) {
        matrix<float> AtWAMatrix_mat=matrix<float>(numbBF);

        float * AtWAMatrix_tmp=&AtWAMatrix[m*numbBF*numbBF];
        bool success=0;
        for (int j1=0; j1<numbBF; j1++) {
            for (int j2=0; j2<numbBF; j2++) {
                AtWAMatrix_mat.setvalue(j1, j2, AtWAMatrix_tmp[j1+j2*numbBF]);
            }
        }
        AtWAMatrix_mat.invert();
        for (int j1=0; j1<numbBF; j1++) {
            for (int j2=0; j2<numbBF; j2++) {
                AtWAMatrix_mat.getvalue(j1, j2, invAtWAMatrix[j1+j2*numbBF+m*numbBF*numbBF], success);
            }
        }
    }
    // Clearing memory
    delete [] AtWAMatrix;
    AtWAMatrix=NULL;
    return invAtWAMatrix;
}

// Returns in a vector of pointers the coefficients for the bias field
float * TreeEM::MakeFinalBFCoeffsChildren(){
    // Initialisation and obtention of the needed values

    int numbmodal=this->GetNumberModalities();
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    float * FinalBFCoeffsResult = new float[numbBF*numbmodal];
    for(int i=0;i<numbBF*numbmodal;i++){
        FinalBFCoeffsResult[i]=0;
    }
    float * invAtWAMatrixVector=this->MakeInvAtWAMatrixChildren();
    float * AtWRVector=this->MakeAtWRVectorChildren();
    for (int m=0; m<numbmodal; m++) {
        float * invAtWA=&invAtWAMatrixVector[m*numbBF*numbBF];
        float * AtWRV_PTR=&AtWRVector[m*numbBF];
        float * invAtWA_PTR=invAtWA;
        for (int j1=0; j1<numbBF; j1++) {
            invAtWA_PTR=&invAtWA[j1*numbBF];
            AtWRV_PTR=&AtWRVector[m*numbBF];
            PrecisionTYPE FinalBFCoeffsResult_tmp=0;
            for (int j2=0; j2<numbBF; j2++,AtWRV_PTR++,invAtWA_PTR++) {
                FinalBFCoeffsResult_tmp+=(PrecisionTYPE)(*AtWRV_PTR)*(*invAtWA_PTR);
            }
            FinalBFCoeffsResult[j1+m*numbBF]=FinalBFCoeffsResult_tmp;
        }
    }
    // Clearing memory
    delete [] invAtWAMatrixVector;
    invAtWAMatrixVector=NULL;
    delete [] AtWRVector;
    AtWRVector=NULL;
//    float * FinalTestBFCoeffs=this->MakeBFCoeffsDirect();
    return FinalBFCoeffsResult;
}

float * TreeEM::MakeDataBFCorrected(){
    int numelmasked=this->GetNumberMaskedElements();
    int numel=this->GetNumberElements();
    int numbmodal=this->GetNumberModalities();
    float * DataBFCorrected=new float[numelmasked*numbmodal];//{0};
    //    for(int i=0;i<numelmasked*numbmodal;i++){
    //        DataBFCorrected[i]=0;
    //    }
    float * DataBFCorrected_PTR=DataBFCorrected;
    float * DataInit=static_cast<float *>(this->GetDataImage()->data);
    float * DataInit_PTR=DataInit;
    int * L2S_PTR= this->GetL2S();
    //    int CountActive=0;
    //    for(int i=0;i<numel;i++,L2S_PTR++){
    //        if(*L2S_PTR>=0){
    //            CountActive++;
    //        }
    //    }
    float * BFCoeffsUsed=this->GetBFCoeffs();
    for (int m=0; m<numbmodal; m++) {
        L2S_PTR=this->GetL2S();
        //        CountActive=0;
        DataBFCorrected_PTR=&DataBFCorrected[m*numelmasked];
        DataInit_PTR=&DataInit[m*numel];
        for (int i=0; i<numel; i++,DataInit_PTR++,L2S_PTR++) {
            if(*L2S_PTR>=0){
                //                CountActive++;
                *DataBFCorrected_PTR=*DataInit_PTR;
                DataBFCorrected_PTR++;
            }
        }
    }
    if(BFCoeffsUsed!=NULL){
        // If there is a correction of the bias field to be done

        float * BFCorrection_new=this->MakeBFCorrection();
        for (int m=0; m<numbmodal; m++) {
            float * BFCorrection_PTR=&BFCorrection_new[m*numelmasked];
            //                    int * L2S_PTR=this->GetL2S();
            DataBFCorrected_PTR=&DataBFCorrected[m*numelmasked];
            for (int i=0; i<numelmasked; i++,DataBFCorrected_PTR++,BFCorrection_PTR++) {
                *DataBFCorrected_PTR-=(*BFCorrection_PTR);
            }
        }
        if(BFCorrection_new!=NULL){
            delete [] BFCorrection_new;
            BFCorrection_new=NULL;
        }
    }
    //        int CountDataCorrectedZero=0;
    //        DataBFCorrected_PTR=DataBFCorrected;
    //        for(int i=0;i<numelmasked*numbmodal;i++,DataBFCorrected_PTR++){
    //            if(*DataBFCorrected_PTR==0){
    //                CountDataCorrectedZero++;
    //            }
    //        }
    return DataBFCorrected;
}


//float * TreeEM::MakeBFCorrection(){

//    float * BFCoeffs_newB=this->GetBFCoeffs();
//    if(BFCoeffs_newB==NULL){
//        return NULL;
//    }
//    float * BFBasisFunctions=this->MakeBasisFunctions();
//    int numelmasked=this->GetNumberMaskedElements();
//    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
//    float BFCoeffs_new[numbBF*MaxNumbModal];
//    for(int i=0;i<numbBF*MaxNumbModal;i++){
//        BFCoeffs_new[i]=BFCoeffs_newB[i];
//    }
//    int numbmodal=this->GetNumberModalities();
//    float * BFCorrectionResult=new float[numbmodal*numelmasked];
//    for(int i=0;i<numbmodal*numelmasked;i++){
//        BFCorrectionResult[i]=0;
//    }
//    if(BFCoeffs_new!=NULL){
//    for(int m=0;m<numbmodal;m++){
//        for(int j=0;j<numbBF;j++){
//            float * BasisFunctions_PTR=&BFBasisFunctions[j*numelmasked];
////            float * BFCoeffs_PTR=&BFCoeffs_new[m*numbBF];
//            float * BFCorrectionResult_PTR=&BFCorrectionResult[m*numelmasked];
//            for(int i=0;i<numelmasked;i++,BFCorrectionResult_PTR++,BasisFunctions_PTR++){
//                *BFCorrectionResult_PTR+=BFCoeffs_new[j+m*numbBF]*(*BasisFunctions_PTR);
//            }
//        }
//    }
//    }
//    delete[] BFBasisFunctions;
//    BFBasisFunctions=NULL;
//    return BFCorrectionResult;
//}



float * TreeEM::MakeBFCorrection(){

    float * BFCoeffs_newB=this->GetBFCoeffs();
    if(BFCoeffs_newB==NULL){
        return NULL;
    }
    //    float * BFBasisFunctions=this->MakeBasisFunctions();
    int numelmasked=this->GetNumberMaskedElements();
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    int numbmodal=this->GetNumberModalities();
    float BFCoeffs_new[MaxNumbBF*MaxNumbModal];
    for(int i=0;i<MaxNumbBF*MaxNumbModal;i++){
        if(i<numbBF*numbmodal){
            BFCoeffs_new[i]=BFCoeffs_newB[i];
        }
        else{
            BFCoeffs_new[i]=0;
        }
    }

    float * BFCorrectionResult=new float[numbmodal*numelmasked];
    for(int i=0;i<numbmodal*numelmasked;i++){
        BFCorrectionResult[i]=0;
    }
    if(BFCoeffs_newB!=NULL){
        float FactorX=(float)(2.0/(this->GetDataImage()->nx-1));
        float FactorY=(float)(2.0/(this->GetDataImage()->ny-1));
        float FactorZ=(float)(2.0/(this->GetDataImage()->nz-1));

        int SizePlane=this->GetDataImage()->nx*this->GetDataImage()->ny;
        int SizeColumn=this->GetDataImage()->nx;

        int XorderBF[MaxNumbBF];
        int YorderBF[MaxNumbBF];
        int ZorderBF[MaxNumbBF];
        int tmpnumbBF=0;
        for(int order=0;order<=BForder;order++){
            for(int zorder=0;zorder<=order;zorder++){
                for(int yorder=0;yorder<=order-zorder;yorder++){
                    int xorder=order-zorder-yorder;
                    XorderBF[tmpnumbBF]=xorder;
                    YorderBF[tmpnumbBF]=yorder;
                    ZorderBF[tmpnumbBF]=zorder;
                    tmpnumbBF++;
                }
            }
        }
        int * S2L_PTR=this->GetS2L();
        for(int i=0;i<numelmasked;i++,S2L_PTR++){


            int zValue=*S2L_PTR/SizePlane;
            int yValue=(*S2L_PTR-zValue*SizePlane)/SizeColumn;
            int xValue=(*S2L_PTR)-zValue*SizePlane-yValue*SizeColumn;

            float xUsed=xValue*FactorX-1;
            float yUsed=yValue*FactorY-1;
            float zUsed=zValue*FactorZ-1;
            for(int m=0;m<numbmodal;m++){
                for(int j=0;j<numbBF;j++){
                    float BFValue=pow_int(xUsed,XorderBF[j])*pow_int(yUsed,YorderBF[j])*pow_int(zUsed,ZorderBF[j]);
                    BFCorrectionResult[i+m*numelmasked]+=BFCoeffs_new[j+m*numbBF]*BFValue;

                }
            }
        }
    }

    return BFCorrectionResult;
}





/******************* SAVING RESULTS ******************************/
void TreeEM::SaveAllClasses(string filenameOut,SEG_PARAMETERS * segment_param){
    int numbTotalClasses=0;
    //this->GetNumberLeaves(numbTotalClasses);
    numbTotalClasses=this->GetNumberAllLeaves();
    vector<TreeEM *> LeavesVector;
    //this->GetLeaves(LeavesVector);
    LeavesVector=this->GetAllLeaves();
    nifti_image * Result=nifti_copy_nim_info(this->GetDataImage());
    Result->dim[0]=4;
    Result->dim[4]=numbTotalClasses;
    Result->dim[5]=1;
    nifti_update_dims_from_array(Result);
    Result->data = (void *) calloc(Result->nvox, sizeof(float));
    int numel=this->GetNumberElements();
    float * Result_PTR=static_cast<float *>(Result->data);
    float * Result_PTRtmp=Result_PTR;
    //float * Input_PTR=LeavesVector[0]->GetPartialResult(NORMRESP);
    for (int c=0; c<numbTotalClasses; c++) {
        Result_PTRtmp=&Result_PTR[c*numel];
        float * Input=LeavesVector[c]->GetPartialResult(NORMRESP,segment_param);
        float * Input_PTR=Input;
        for (int i=0; i<numel; i++, Result_PTRtmp++,Input_PTR++) {
            *Result_PTRtmp=*Input_PTR;
        }
        if(Input!=NULL){
            delete[]Input;
            Input=NULL;
        }
    }
    nifti_set_filenames(Result, filenameOut.c_str(), 0, 0);
    nifti_image_write(Result);
    nifti_image_free(Result);
//    LeavesVector.clear();
}

void TreeEM::SaveGeneralClasses(string filenameOut,SEG_PARAMETERS * segment_param){
    int numbGeneralClasses=this->GetNumberChildren();
    vector<TreeEM *> ChildrenVector=this->GetChildren();
    nifti_image * Result=nifti_copy_nim_info(this->GetDataImage());
    Result->dim[0]=4;
    Result->dim[4]=numbGeneralClasses;
    Result->dim[5]=1;
    nifti_update_dims_from_array(Result);
    Result->data = (void *) calloc(Result->nvox, sizeof(float));
    int numel=this->GetNumberElements();
    float * Result_PTR=static_cast<float *>(Result->data);
    float * Result_PTRtmp=Result_PTR;
    //float * Input_PTR=ChildrenVector[0]->GetPartialResult(NORMRESP);
    for (int c=0; c<numbGeneralClasses; c++) {
        Result_PTRtmp=&Result_PTR[c*numel];
        float *Input=ChildrenVector[c]->GetPartialResult(NORMRESP,segment_param);
        float * Input_PTR=Input;
        for (int i=0; i<numel; i++, Result_PTRtmp++,Input_PTR++) {
            *Result_PTRtmp=*Input_PTR;
        }
        if(Input!=NULL){
            delete [] Input;
            Input=NULL;
        }
    }
    nifti_set_filenames(Result, filenameOut.c_str(), 0, 0);
    nifti_image_write(Result);
    nifti_image_free(Result);
    ChildrenVector.clear();
}

void TreeEM::SaveBFBasisFunctions(string filenameBF){
    int numelmasked=this->GetNumberMaskedElements();
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    float * BasisFunctionsToSave=this->MakeBasisFunctions();
    nifti_image * Result=nifti_copy_nim_info(this->GetDataImage());
    Result->dim[0]=4;
    Result->dim[4]=numbBF;
    Result->dim[5]=1;
    nifti_update_dims_from_array(Result);
    Result->data = (void *) calloc(Result->nvox, sizeof(float));
    int numel=this->GetNumberElements();
    float * Result_PTR=static_cast<float *>(Result->data);
    float * Result_PTRtmp=Result_PTR;
    float * Input_PTR=BasisFunctionsToSave;
    int * L2S_PTR=this->GetL2S();
    for (int j=0; j<numbBF; j++) {
        Result_PTRtmp=&Result_PTR[j*numel];
        Input_PTR=&BasisFunctionsToSave[j*numelmasked];
        L2S_PTR=this->GetL2S();
        for (int i=0; i<numel; i++, Result_PTRtmp++,L2S_PTR++) {
            if (*L2S_PTR>=0) {
                *Result_PTRtmp=*Input_PTR;
                Input_PTR++;
            }
            else{
                *Result_PTRtmp=0;
            }

        }
    }
    if(BasisFunctionsToSave!=NULL){
        delete[] BasisFunctionsToSave;
        BasisFunctionsToSave=NULL;
    }
    nifti_set_filenames(Result, filenameBF.c_str(), 0, 0);
    nifti_image_write(Result);
    nifti_image_free(Result);
}

void TreeEM::SaveBFBasisFunctions(int BFnumb,string filenameBF){
    int numelmasked=this->GetNumberMaskedElements();
    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    if (BFnumb>=numbBF) {
        return;
    }
    float * BasisFunctionsToSave=this->MakeBasisFunctions();
    nifti_image * Result=nifti_copy_nim_info(this->GetDataImage());
    Result->dim[0]=3;
    Result->dim[4]=1;
    Result->dim[5]=1;
    nifti_update_dims_from_array(Result);
    Result->data = (void *) calloc(Result->nvox, sizeof(float));
    int numel=this->GetNumberElements();
    float * Result_PTR=static_cast<float *>(Result->data);
    float * Result_PTRtmp=Result_PTR;
    float * Input_PTR=&BasisFunctionsToSave[BFnumb*numelmasked];
    int * L2S_PTR=this->GetL2S();


    for (int i=0; i<numel; i++, Result_PTRtmp++,L2S_PTR++) {
        if (*L2S_PTR>=0) {
            *Result_PTRtmp=*Input_PTR;
            Input_PTR++;
        }
        else{
            *Result_PTRtmp=0;
        }

    }
    if(BasisFunctionsToSave!=NULL){
        delete[]BasisFunctionsToSave;
        BasisFunctionsToSave=NULL;
    }
    nifti_set_filenames(Result, filenameBF.c_str(), 0, 0);
    nifti_image_write(Result);
    nifti_image_free(Result);
}

void TreeEM::SaveBFCorrection(string filenameBF){
    int numelmasked=this->GetNumberMaskedElements();
    //    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    int numbmodal=this->GetNumberModalities();

        if (this->AreBFCoeffsDirectUseable()) {
            //                float * BFCoeffs_tmp=BFCoeffsToSave;
            float * BFCorrection=this->MakeBFCorrection();
            nifti_image * Result=nifti_copy_nim_info(this->GetDataImage());
            Result->dim[0]=4;
            Result->dim[4]=numbmodal;
            Result->dim[5]=1;
            nifti_update_dims_from_array(Result);
            Result->data = (void *) calloc(Result->nvox, sizeof(float));
            int numel=this->GetNumberElements();
            int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
            float * Result_PTR=static_cast<float *>(Result->data);
            float * Result_PTRtmp=Result_PTR;
            for (int i=0; i<numel; i++, Result_PTRtmp++) {
                *Result_PTRtmp=0;
            }
            //                float * Input_PTR=BasisFunctionsToSave;
            int * L2S_PTR=this->GetL2S();

            for (int m=0; m<numbmodal; m++) {
                for (int j=0; j<numbBF; j++) {
                    float * BFCorrection_PTR=&BFCorrection[m*numelmasked];
                    Result_PTRtmp=&Result_PTR[m*numel];
                    //                        Input_PTR=&BasisFunctionsToSave[j*numelmasked];
                    L2S_PTR=this->GetL2S();
                    for (int i=0; i<numel; i++, Result_PTRtmp++,L2S_PTR++) {
                        if (*L2S_PTR>=0) {
                            *Result_PTRtmp=*BFCorrection_PTR;
                            BFCorrection_PTR++;
                        }
                        else{
                            *Result_PTRtmp=0;
                        }

                    }

                }
            }
            nifti_set_filenames(Result, filenameBF.c_str(), 0, 0);
            nifti_image_write(Result);
            nifti_image_free(Result);
            if(BFCorrection!=NULL){
                delete[] BFCorrection;
                BFCorrection=NULL;
            }
        }



}

void TreeEM::SaveBFCorrectedData(string filenameBF){
    //    int numelmasked=this->GetNumberMaskedElements();
    //    int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
    int numbmodal=this->GetNumberModalities();
    int numelmasked=this->GetNumberMaskedElements();
    if (this->AreBFCoeffsDirectUseable() || !this->AreBFCoeffsDirectUseable()) {
        nifti_image * Result=nifti_copy_nim_info(this->GetDataImage());
        Result->dim[0]=4;
        Result->dim[4]=numbmodal;
        Result->dim[5]=1;
        nifti_update_dims_from_array(Result);
        Result->data = (void *) calloc(Result->nvox, sizeof(float));
        int numel=this->GetNumberElements();
        float * Result_PTR=static_cast<float *>(Result->data);
        float * Result_PTRtmp=Result_PTR;
        for (int i=0; i<numel; i++, Result_PTRtmp++) {
            *Result_PTRtmp=0;
        }
        float * Input=this->GetDataBFCorrected();
        float * Input_PTR=Input;
        int * L2S_PTR=this->GetL2S();

        for (int m=0; m<numbmodal; m++) {
            L2S_PTR=this->GetL2S();
            Result_PTRtmp=&Result_PTR[m*numel];
            Input_PTR=&Input[m*numelmasked];

            for (int i=0; i<numel; i++, Result_PTRtmp++,L2S_PTR++) {
                if(*L2S_PTR>=0){
                    *Result_PTRtmp=*Input_PTR;
                    Input_PTR++;
                }


            }

        }

        nifti_set_filenames(Result, filenameBF.c_str(), 0, 0);
        nifti_image_write(Result);
        nifti_image_free(Result);
    }

}

//void TreeEM::SaveTreeInTextFile(string filenameOut){
//    std::ofstream fs(filenameOut);

//    if(!fs)
//    {
//        std::cout<<"Cannot open the output file."<<std::endl;
//        return;
//    }
//    int numbchild=this->GetNumberChildren();
//    int numbmodal=this->GetNumberModalities();
//    float indFactor=this->GetIndFactor();
//    fs<<"There are "<<numbchild<<" general classes \n";
//    fs<<"The studied data contains "<<numbmodal<<" modalities\n";
//    fs<<"The factor for independence is "<<indFactor<<" \n";
//    // Writing BF coefficients if exist
//    if (this->AreBFCoeffsDirectUseable()) {
//        int numbBF=((BForder+1)*(BForder+2)*(BForder+3))/6;
//        float * BFCoeffs_new=this->GetBFCoeffs();
//        fs<<" The BF coefficients are : \n";
//        for (int m=0; m<numbmodal; m++) {
//            float * BFCoeffs_PTR=&BFCoeffs_new[m*numbBF];
//            fs << "For modality "<<m<<" \n";
//            for (int l=0; l<numbBF; l++,BFCoeffs_PTR++) {
//                fs<<*BFCoeffs_PTR<<"   -    ";
//            }
//            fs<<"\n";
//        }
//    }
//    // Writing parameters and weights of the different classes
//    for (int c=0; c<numbchild; c++) {
//        int numbLeaves=this->GetChild(c)->GetNumberAllLeaves();
//        vector<TreeEM *> LeavesVector=this->GetChild(c)->GetAllLeaves();
//        fs<<"General class "<<c<<" contains "<<numbLeaves<<" subclasses and its weight is "<<this->GetChild(c)->GetNormWeight()<<" \n";
//        if (numbLeaves==0) {
//            fs<<"The general class "<<c<<" is not a mixture and its parameters are : \n Mean :\n";
//            float * Mean=this->GetChild(c)->GetMean();
//            float * Variance=this->GetChild(c)->GetVariance();
//            for (int m=0; m<numbmodal; m++) {
//                fs<<Mean[m]<<"    ";
//            }
//            fs<<"\n Variance \n";
//            for (int m1=0; m1<numbmodal; m1++) {
//                for (int m2=0; m2<numbmodal; m2++) {
//                    fs<<Variance[m1+m2*numbmodal]<<"     ";
//                }
//                fs<<"\n";
//            }

//        }
//        else{
//            for (int l=0; l<numbLeaves; l++) {
//                fs<<"The subclass "<<l<<" of general class "<<c<<" is of weight "<<LeavesVector[l]->GetNormWeight()<<" \n";
//                fs<<"The parameters of the distribution are : \n Mean \n";
//                float * Mean=LeavesVector[l]->GetMean();
//                float * Variance=LeavesVector[l]->GetVariance();
//                for (int m=0; m<numbmodal; m++) {
//                    fs<<Mean[m]<<"    ";
//                }
//                fs<<"\n Variance \n";
//                for (int m1=0; m1<numbmodal; m1++) {
//                    for (int m2=0; m2<numbmodal; m2++) {
//                        fs<<Variance[m1+m2*numbmodal]<<"     ";
//                    }
//                    fs<<"\n";
//                }
//            }
//        }
//    }
//    //fs<<"ghgh";
//    fs.close();
//    return;
//}
vector<TreeEM*> TreeEM::GetAllTreesFromLevel(int l){
    vector<TreeEM*> LevelElements;
    // First check if requested level is available
    if(l>this->GetNumberLevels()){
        cout<< "Requested level more than possible depth"<<endl;
        return LevelElements;
    }
    if(l==0){
        LevelElements.push_back(this);

    }
    else{
    int numbchild=this->GetNumberChildren();
    for(int c=0;c<numbchild;c++){
        vector<TreeEM *> LevelElements_child=this->GetChild(c)->GetAllTreesFromLevel(l-1);
        for(int c1=0;c1<LevelElements_child.size();c1++){
            LevelElements.push_back(LevelElements_child[c1]);
        }
    }
    }
    return LevelElements;
}

int TreeEM::GetNumberLevels(){
      int MaxLevel=0;
      int numbchild=this->GetNumberChildren();
//      cout<< "considered numbchild is "<<numbchild<<endl;
      for(int c=0;c<numbchild;c++){
//          cout << "the child looked is "<< c<<endl;
          int Level=(this->GetChild(c))->GetNumberLevels();
//          cout<< "Number of levels of child "<<c<<" is "<<Level<<endl;
          if (Level>MaxLevel){
              MaxLevel=Level;
          }
      }
//      cout<<"Number of levels is "<<MaxLevel+1<<" for "<<this<<endl;
      return MaxLevel+1;
    }

int TreeEM::GetLevel(){
    int Level=0;
    if(!(this->IsRoot()||this->WhatTypeOfTreeComponent()==INITIALNODE)){
        return 1+this->GetParent()->GetLevel();
    }
    return Level;
}

int * TreeEM::GetHierarchy(){
    int numbLevel=this->GetLevel();
    int * Hierarchy=NULL;
    if(numbLevel>0){
    Hierarchy=new int[numbLevel];
    int * Hierarchy_bef=this->GetParent()->GetHierarchy();
    int Level_bef=this->GetParent()->GetLevel();
    for(int l=0;l<Level_bef;l++){
        Hierarchy[l]=Hierarchy_bef[l];
    }
    Hierarchy[Level_bef]=this->GetParent()->FindIndex(this);
    if(Hierarchy_bef!=NULL){
        delete [] Hierarchy_bef;
        Hierarchy_bef=NULL;
    }
    }
    else{
        Hierarchy=NULL;
    }
    return Hierarchy;
}

int TreeEM::FindIndex(TreeEM* ChildToFind){
    int numbchild=this->GetNumberChildren();
    for(int c=0;c<numbchild;c++){
        if(this->GetChild(c)==ChildToFind){
            return c;
        }
    }
    return -1;
}

int TreeEM::FindGeneralClass(){
    // Check if Root or initial node
    if(this->IsRoot()){
        return -1;
    }
    else{
        if (this->GetParent()->IsRoot()){
            return this->GetParent()->FindIndex(this);
        }
        else{
            return this->GetParent()->FindGeneralClass();
        }
    }
}

TreeEM* TreeEM::FindGeneralClassPriors(){
    // Check if truly has upwards some priors
    if(this->GetPriors()==NULL){
        return NULL;
    }
    else{
        if(this->GetPriorsDirect()!=NULL){
            return this;
        }
        else{
            return this->GetParent()->FindGeneralClassPriors();
        }
    }
}

TreeEM* TreeEM::FindMainNode(){
    // if at Root returns NULL
    int * Hierarchy=this->GetHierarchy();
    if(Hierarchy==NULL){ // => we are considering the root
        return NULL;
    }
    else if(this->GetFlagOutliers()!=1){ // no outlier model or outlier model 2 not trilayered, then main class at first level
        return this->FindRoot()->GetChild(Hierarchy[0]);
    }
    else{ // outlier model case
        if(this->GetLevel()==1 && this!=this->GetNodeOutlier()) { // level 1 above main classes
            return NULL;
        }
        else if(this==this->GetNodeOutlier()){ // level 1 at node outlier
            return this;
        }
        else{ // under main node follow hierarchy : main node at level 2
            return this->FindRoot()->GetChild(Hierarchy[0])->GetChild(Hierarchy[1]);
        }
    }
}


void TreeEM::SaveTreeInTextFile(string filenameOut,SEG_PARAMETERS * segment_param){
    std::ofstream fs(filenameOut.c_str());
    cout<< filenameOut<< " to save tree in text file"<<endl;
    if(!fs.is_open())
    {
        std::cout<<"Cannot open the output file."<<std::endl;
        return;
    }
    int numbLevels=this->GetNumberLevels();
    int numbchild=this->GetNumberChildren();
    int numbmodal=this->GetNumberModalities();
    fs << "Level 0 1"<<endl;
    fs << "Images "<<numbmodal<<endl;
    for(int m=0;m<numbmodal;m++){
        string filename=segment_param->filename_input[m];
        string str(filename);
        fs<<filename<<endl;
    }
    if(segment_param->flag_mask && this->GetMask()!=NULL){
        fs<<"Mask "<<this->GetMask()->iname<<endl;
    }
    if(segment_param->flag_MRFOut){
        fs << "MRFImage "<<segment_param->filename_MRFOut<<endl;
    }
    if(segment_param->flag_GMatrixIn){
        fs<< "GMatrixIn "<<segment_param->filename_GMatrix<<endl;
    }
    if(segment_param->flag_GMatrixPost){
        fs<< "GMatrixPost "<<segment_param->filename_GMatrixPost<<endl;
    }
    if(segment_param->flag_Outliers){
        fs<<"OutliersMode "<<segment_param->OutliersMod<<endl;
        fs<<"OutliersWeight "<<segment_param->OutliersWeight<<endl;
    }
    else{
        fs<<"NoOutlier"<<endl;
    }
    fs<<"AtlasWeight "<<segment_param->AtlasWeight<<endl;
    fs<<"AtlasSmoothing "<<segment_param->AtlasSmoothing<<endl;
    fs<<"KernelSize "<<segment_param->KernelSize<<endl;
    if(segment_param->bias_order > 0){
        fs<<"BF "<<segment_param->bias_order<<endl;
        int numbBF=(int)((BForder+1)*(BForder+2)*(BForder+3))/6;
        float * BFCoeffs_toprint=this->GetBFCoeffs();
        for(int m=0;m<numbmodal;m++){
            for(int j=0;j<numbBF;j++){
                fs<<BFCoeffs_toprint[m*numbBF+j]<<"     ";
            }
            fs<<endl;
        }
    }

    if(segment_param->flag_DP){
        float * DPChildrenToPrint=this->GetDPChildren();
        if(segment_param->flag_DistClassInd){
        for(int c=0;c<numbchild;c++){

            fs<<"CP "<<c<<endl;
            for(int i=0;i<MaxSupport;i++){
                fs<<DPChildrenToPrint[i+c*MaxSupport]<< "   ";
            }
            fs<<endl;
        }
        }
    }
    fs<<"IF "<<this->GetIndFactor()<<endl;
    fs<<endl;
    int pa=0;
    for(int l=1;l<numbLevels;l++){
        vector<TreeEM*> LevelsElement=this->GetAllTreesFromLevel(l);
        int numbElements=LevelsElement.size();
        fs<<"Level "<<l<<" "<<numbElements<<endl;
        fs<<endl;
        for(int c=0;c<numbElements;c++){
            TreeEM * ElementToPrint=LevelsElement[c];
            int * Hierarchy=ElementToPrint->GetHierarchy();
            fs<<"Class ";
            for(int i=0;i<l;i++){
                fs<<Hierarchy[i]<<" ";
            }
            fs<<endl;
            if(Hierarchy!=NULL){
                delete [] Hierarchy;
                Hierarchy=NULL;
            }
            if(ElementToPrint->GetPriorsDirect()!=NULL){
                fs<<"Prior "<<ElementToPrint->GetPriorsDirect()->iname<<endl;
            }
            if(ElementToPrint->GetPriorsAdaptedDirect()!=NULL && pa<segment_param->filename_PriorsAdaptedOut.size()){
                fs<<"AdaptedPriors"<<segment_param->filename_PriorsAdaptedOut[pa]<<endl;
                pa++;
            }
            fs << "Weight "<<ElementToPrint->GetNormWeight()<<endl;
            int DistributionType=ElementToPrint->GetDistributionType();
            switch(DistributionType){
            case 0 :
                fs<<"Mixture"<<endl;
                break;
            case 2:
                fs<<"UniformDist"<<endl;
                break;
            default:{
                fs<<"Leaf "<<DistributionType<<endl;
                float * ValueToPrint=ElementToPrint->GetParametersValue();
                fs<<"Mean ";
                for(int m=0;m<numbmodal;m++){
                    fs<<ValueToPrint[m]<<" ";
                }
                fs<<endl;
                fs<<"Variance ";
                for(int m=0;m<numbmodal*numbmodal;m++){
                    fs<<ValueToPrint[numbmodal+m]<<" ";
                }
                fs<<endl;
            }
            }
            fs<<endl;
        }
    }
    fs.close();
    return;
}

void TreeEM::ModifyCountFiles(vector<string> CountFiles){
    if(this->GetNumberChildren()!=CountFiles.size()){
        cout<<"Pb in the number of count files and general classes"<<endl;
        return;
    }
    else{
        int numbchild=this->GetNumberChildren();
        for(int c=0;c<numbchild;c++){
            int numbLeaves=(this->GetChild(c)->GetNumberAllLeaves()>0?this->GetChild(c)->GetNumberAllLeaves():1);

            ifstream CountFile(CountFiles[c].c_str());
            ofstream TempCount("/Users/Carole/Documents/PhD/temporaryCount.txt");
            if(!CountFile || !TempCount){
                cout<<"Pb in opening files"<<endl;
            }
            int valueFile;
            int line=0;
            while(CountFile >> valueFile){
                if(line==numbLeaves-1){
                    valueFile+=1;
                }
                TempCount << valueFile<<"\n";
                line++;
            }
            rename("/Users/Carole/Documents/PhD/temporaryCount.txt",CountFiles[c].c_str());
        }
    }
}

float * TreeEM::GaussianBlurring(float * CountHistogram, float gauss_std, vector<int> dim){
    int numbDim=dim.size();
    int kernelsizemin=(int)floorf(gauss_std*KernelSize);
    int kernelsizemax=(int)ceilf(gauss_std*KernelSize);
    int kernelsize=0;
    (kernelsizemin/2)*2==kernelsizemin?kernelsize=kernelsizemax:kernelsize=kernelsizemin;
    // To do filtering, kernel must be at least of size 2.
    kernelsize=kernelsize<3?3:kernelsize;

    // Construction of the kernel
    float * Kernel=new float[kernelsize];
    int kernelradius=kernelsize/2;
    for(int i =0;i<kernelsize;i++){
        Kernel[i]=1.0/(sqrt(2.0*M_PI)*gauss_std)*expf(-0.5*pow((i-kernelradius)/gauss_std,2));
    }
    //Normalisation of the kernel
    float SumKernel=0;
    for(int i=0;i<kernelsize;i++){
        SumKernel+=Kernel[i];
    }
    if(SumKernel>0){
    for(int i=0;i<kernelsize;i++){
        Kernel[i]/=SumKernel;
    }
    }

    // Construction of the Shift array
    int * Shift=new int[numbDim];
    for(int c=0;c<numbDim;c++){
        if(c>0){
        Shift[c]=Shift[c-1]*dim[c-1];
        }
        else{
            Shift[c]=1;
        }
    }

    // Copying the initial array to blur
    int SizeHistogram=Shift[numbDim-1]*dim[numbDim-1];
    float * Gaussian_tmp=new float[SizeHistogram];
    float * Gaussian_tmp2=new float[SizeHistogram];
    int CountNonZero=0;
    int CountZeroB=0;
    for(int i=0;i<SizeHistogram;i++){
        Gaussian_tmp[i]=CountHistogram[i];
        if(Gaussian_tmp[i]>0){
            CountNonZero++;
        }
    }

    for(int c=0;c<numbDim;c++){ // To do the blurring in each direction
        for(int j=0;j<SizeHistogram/(dim[c]);j++){
            int i=0;
//            int RemainDim[numbDim-1];
//            int IndRemainDim[numbDim-1];
//            int RemainShift[numbDim-1];
//            int FinalShift[numbDim-1];
            int * RemainDim=new int[numbDim-1];
            int * IndRemainDim=new int[numbDim-1];
            int * RemainShift=new int[numbDim-1];
            int * FinalShift=new int[numbDim-1];
            for(int i=0;i<numbDim-1;i++){
                RemainShift[i]=1;
                FinalShift[i]=1;
            }
            for(int d=0;d<numbDim;d++){
                if(d!=c){
                    RemainDim[i]=dim[d];
                    IndRemainDim[i]=d;
                    i++;
                }
            }
            for(int i=0;i<numbDim-1;i++){
                if(i>0){
                    RemainShift[i]=RemainShift[i-1]*RemainDim[i-1];
                }
            }
            int t=j;
            int tmp=0;
            int * Index=new int[numbDim-1];
            for(int i=numbDim-2;i>=0;i--){
                tmp=t/RemainShift[i];
                t=t-tmp*RemainShift[i];
                Index[i]=tmp;
            }
            for(int i=0;i<numbDim-1;i++){
                if(c<IndRemainDim[i]){
                    FinalShift[i]=dim[c]*RemainShift[i];
                }
                else{
                    FinalShift[i]=RemainShift[i];
                }
            }
            int FinalIndex=0;
            for(int i=0;i<numbDim-1;i++){
                FinalIndex+=Index[i]*FinalShift[i];
            }
            for(int i=0;i<dim[c];i++){
                float tmpkernel=0;
                 int FinalIndexB=FinalIndex+i*Shift[c];
                for(int k=0;k<kernelsize;k++){
                    if(i-kernelsize/2+k>=0 && i-kernelsize/2+k<dim[c]){
                        tmpkernel+=Kernel[k]*Gaussian_tmp[FinalIndexB+(-kernelsize/2+k)*Shift[c]];
                    }
                }
                Gaussian_tmp2[FinalIndexB]=tmpkernel;
            }
            delete [] Index;
            Index=NULL;
            if(RemainShift!=NULL){
                delete [] RemainShift;
                RemainShift=NULL;
            }
            if(RemainDim!=NULL){
                delete [] RemainDim;
                RemainDim=NULL;
            }
            if(IndRemainDim!=NULL){
                delete [] IndRemainDim;
                IndRemainDim=NULL;
            }
            if(FinalShift!=NULL){
                delete [] FinalShift;
                FinalShift=NULL;
            }
        }
        float sumNorm=0;

        for(int l=0;l<SizeHistogram;l++){
            if(Gaussian_tmp2[l]>0){
                CountZeroB++;
            }
            sumNorm+=Gaussian_tmp2[l];
            Gaussian_tmp[l]=Gaussian_tmp2[l];
        }
//        for(int l=0;l<SizeHistogram;l++){
//            Gaussian_tmp[l]/=sumNorm;
//        }

    }
    delete [] Gaussian_tmp2;
    Gaussian_tmp2=NULL;
    delete [] Kernel;
    Kernel=NULL;
    delete [] Shift;
    Shift=NULL;
    return Gaussian_tmp;
}


/******************** METHODS RELATED TO MRF CALCULATION *******************/

float * TreeEM::GetNormRespShifted(int dim, int dir){
    // First make sure that the input values are ok or change them accordingly
    if(dir!=0){
    dir=dir>0?1:-1;
    }
    dim=dim>2?2:dim;
    dim=dim<0?0:dim;

    // Initialisation of the needed values
    int Dimensions[3];
    Dimensions[0]=this->GetDataImage()->nx;
    Dimensions[1]=this->GetDataImage()->ny;
    Dimensions[2]=this->GetDataImage()->nz;
    int Shift[3];
    Shift[0]=1;
    Shift[1]=this->GetDataImage()->nx;
    Shift[2]=this->GetDataImage()->nx*this->GetDataImage()->ny;
    int numbx=this->GetDataImage()->nx;
    int numby=this->GetDataImage()->ny;
    int numbz=this->GetDataImage()->nz;
    int indexShift=Shift[dim]*dir;
    int numbmodal=this->GetNumberModalities();
    int numel=this->GetNumberElements();
    int nvox=this->GetDataImage()->nvox;
    int numelmasked=this->GetNumberMaskedElements();
    int newIndex=0;
    int oldIndex=0;
    // Allocation of memory for shifted matrix
    float * NormRespShifted= new float[numel];//{0};
    int * L2S_PTR=this->GetL2S();
    float * NormRespPointedNC=this->GetNormResp();
    float * NormRespPointed = new float[numel];
        L2S_PTR=this->GetL2S();
        float * NormRespPointedNC_PTR = NormRespPointedNC;
        for (int i=0; i<numel; i++,L2S_PTR++) {

                NormRespPointed[i]=0;
                if(*L2S_PTR>=0){
                    NormRespPointed[i]=*NormRespPointedNC_PTR;
                    NormRespPointedNC_PTR++;
                }
                NormRespShifted[i]=0;
        }

        for (int i=0; i<numbx; i++) {
            for (int j=0; j<numby; j++) {
                for (int k=0; k<numbz; k++) {
                    newIndex=i+j*Shift[1]+k*Shift[2]+indexShift;
                    oldIndex=i+j*Shift[1]+k*Shift[2];
                    if (newIndex>0 && newIndex<numel) {
                        NormRespShifted[newIndex]=NormRespPointed[oldIndex];
                    }
                }
            }
        }

//    delete [] Dimensions;
//    Dimensions=NULL;
//    delete [] Shift;
//    Shift = NULL;
    delete [] NormRespPointed;
    NormRespPointed = NULL;
    return NormRespShifted;
}

float * TreeEM::GetSumNeighboursNormResp(){
    // First initialisation of the final NormResp obtained
    int numel=this->GetNumberElements();
    int numelmasked=this->GetNumberMaskedElements();
    float * SumNeighboursNormResp=new float[numelmasked];
    for(int i=0;i<numelmasked;i++){
        SumNeighboursNormResp[i]=0;
    }

    // Assuming that it is a 3D image
    for(int dim = 0;dim<3;dim++){
        for(int dir = -1;dir<=1;dir=dir+2){
            // Calculation of the corresponding NormRespShifted
            float * TmpShiftedNormResp=this->GetNormRespShifted(dim,dir);
//            SaveTmpResult(TmpShiftedNormResp,"/Users/Carole/Documents/PhD/TestShifted.nii.gz");
            float * SumNeighboursNormResp_PTR=SumNeighboursNormResp;
            int * L2S_PTR=this->GetL2S();
            // Then Summed in the final result
            for(int i=0;i<numel;i++,L2S_PTR++){
                if (*L2S_PTR>=0){
                    *SumNeighboursNormResp_PTR=TmpShiftedNormResp[i];
                    SumNeighboursNormResp_PTR++;
                }

            }
            // Clearing memory allocated temporarily
            delete [] TmpShiftedNormResp;
            TmpShiftedNormResp=NULL;
        }
    }
//SaveTmpResultMasked(SumNeighboursNormResp,"/Users/Carole/Documents/PhD/TestShifted.nii.gz");
    return SumNeighboursNormResp;
}

float * TreeEM::GMulSumNeighbNormRespExp(){
    // Initialisation of result
    int numelmasked=this->GetNumberMaskedElements();
    int numbLeaves=this->GetNumberAllLeaves();
    vector<TreeEM *> LeavesVector=this->GetAllLeaves();
    float * LeavesSumNeighbour = new float[numelmasked * numbLeaves];

    // Filling LeavesSumNeighbour
    for(int l=0;l<numbLeaves;l++){
        float * tmpLeavesSumNeighbour = LeavesVector[l]->GetSumNeighboursNormResp();
        for(int i=0;i<numelmasked;i++){
            LeavesSumNeighbour[i+l*numelmasked]=tmpLeavesSumNeighbour[i];
        }


        delete [] tmpLeavesSumNeighbour;
        tmpLeavesSumNeighbour=NULL;

    }



   // Initialising before multiplication

    float * Gmatrix=this->GetGMatrix();
    float * ResultMul=new float[numelmasked * numbLeaves];
    for(int i=0;i<numelmasked*numbLeaves;i++){
        ResultMul[i]=0;
    }

    for(int l=0;l<numbLeaves;l++){
            for(int i=0;i<numelmasked;i++){
                for(int l2=0;l2<numbLeaves;l2++){
                    ResultMul[i+l*numelmasked]+=Gmatrix[l*numbLeaves+l2]*LeavesSumNeighbour[i+l2*numelmasked];
            }
        }
    }
    int CountNan=0;
    for(int i=0;i<numelmasked;i++){
        for(int l=0;l<numbLeaves;l++){
            if(LeavesSumNeighbour[i+l*numelmasked]!=LeavesSumNeighbour[i+l*numelmasked]){
                CountNan++;
            }
        }
    }
    delete [] LeavesSumNeighbour;
    LeavesSumNeighbour=NULL;
 CountNan=0;

    for(int l=0;l<numbLeaves;l++){
        for(int i=0;i<numelmasked;i++){
            ResultMul[i+l*numelmasked]=exp(-ResultMul[i+l*numelmasked]);
//            if(ResultMul[i+l*numelmasked]!=ResultMul[i+l*numelmasked]){
//                CountNan++;
//            }
        }
    }

    return ResultMul;

}

//float * TreeEM::SumExpUmrf(){
//    // Check if you are at root level to make this request
//    if(!this->IsRoot()){
//        cout<<"Sum Umrf can only be asked at root"<<endl;
//        return NULL;
//    }
//    int numbLeaves=this->GetNumberAllLeaves();
//    int numelmasked=this->GetNumberMaskedElements();
//    vector<TreeEM *> LeavesVector=this->GetAllLeaves();
//    float * SumExpUmrfResult=new float[numelmasked];
//    for(int i=0;i<numelmasked;i++){
//        SumExpUmrfResult[i]=0;
//    }

//    for(int l=0;l<numbLeaves;l++){
//        float * PartExpUmrfResult=LeavesVector[l]->GetMRF();
//        for(int i=0;i<numelmasked;i++){
//            SumExpUmrfResult[i]+=PartExpUmrfResult[i];
//        }
//    }

//    return SumExpUmrfResult;
//}

void TreeEM::UpdateMRF(){

    float * tmpGMulSumNeigh=this->FindRoot()->GMulSumNeighbNormRespExp();
    int numelmasked=this->GetNumberMaskedElements();
    int numbLeaves=this->FindRoot()->GetNumberAllLeaves();
    vector<TreeEM*> LeavesVector=this->FindRoot()->GetAllLeaves();
    float * SumMRF= new float[numelmasked];
    for(int i=0;i<numelmasked;i++){
        SumMRF[i]=0;
    }
    // Creating the normalisation factor for each active voxel
    int CountNegativePB=0;
    int CountNan=0;
    for(int l=0;l<numbLeaves;l++){
    for(int i=0;i<numelmasked;i++){
        if(tmpGMulSumNeigh[i+l*numelmasked]<0){
            CountNegativePB++;
        }
        if(tmpGMulSumNeigh[i+l*numelmasked]!=tmpGMulSumNeigh[i+l*numelmasked]){
            CountNan++;
        }
        SumMRF[i]+=tmpGMulSumNeigh[i+l*numelmasked];
    }
    }
//    cout<<"Number of negative pb in sum MRF are "<<CountNegativePB<<" and Nan "<<CountNan<<endl;

    // For each leave, creating the MRF value, normalising it and set it
    for(int l=0;l<numbLeaves;l++){
        float * MRFToUpdate= new float[numelmasked];
        int CountMRFNormPb=0;
        for(int i=0;i<numelmasked;i++){
            if(SumMRF[i]>-10E-6){
            MRFToUpdate[i]=tmpGMulSumNeigh[i+l*numelmasked]/SumMRF[i];
            }
            else{
                CountMRFNormPb++;
//                cout<<"Pb with normalisation factor of MRF"<<endl;
                MRFToUpdate[i]=1.0/numbLeaves;
            }
        }
        if(CountMRFNormPb>0){
            cout<< "MRF Norm Pb is "<<CountMRFNormPb<<" at leave "<<l<<endl;
        }
        LeavesVector[l]->SetMRF(MRFToUpdate);
        delete [] MRFToUpdate;
        MRFToUpdate=NULL;
    }

    // Clearing memory
    delete [] SumMRF;
    SumMRF=NULL;
    delete [] tmpGMulSumNeigh;
    tmpGMulSumNeigh=NULL;
}

void TreeEM::SetMRF(float * MRFToUpdate){
    // Check if MRF to update is on a leaf
    if(!this->IsLeaf()){
//        cout<<"MRF is only set on leaves"<<endl;
        if(this->MRF!=NULL){
            delete [] this->MRF;
            this->MRF=NULL;
        }
        return;
    }
    // prepare by deleting already present MRF and creating memory for new one
    if(this->GetMRF()!=NULL){
        delete [] this->GetMRF();
        this->MRF=NULL;
    }
    if(MRFToUpdate!=NULL){
    int numelmasked=this->GetNumberMaskedElements();
    this->MRF=new float[numelmasked];
    for(int i=0;i<numelmasked;i++){
        this->MRF[i]=MRFToUpdate[i];
    }
    }
    return;
}

float * TreeEM::GetMRF(){
    if(this->IsLeaf()){
        return this->MRF;
    }
    else{
        if(this->MRF!=NULL){
            cout<<"There should not be any MRF at this stage"<<endl;
            this->SetMRF(NULL);
        }
        return NULL;
    }
}

float * TreeEM::CopyMRF(){
    float *MRFToCopy=this->GetMRF();
    if(MRFToCopy==NULL){
        return NULL;
    }
    else {
        int numelmasked=this->GetNumberMaskedElements();
        float * MRFCopied = new float[numelmasked];
        for(int i=0;i<numelmasked;i++){
            MRFCopied[i]=MRFToCopy[i];
        }
        return MRFCopied;
    }
}

float * TreeEM::PrepareGMatrixFromFile(string GMatrixFilename){
    // No preparation of GMatrix if not at root
    if(!this->IsRoot()){
        cout<<"No preparation of GMatrix if not at root"<<endl;
        return NULL;
    }

    // First reading and taking all needed information from file
    ifstream text (GMatrixFilename.c_str());
//    text.open(TreeTextFile,ios::in);
    if(!text){
        std::cout<<"could not open the text file properly !"<<endl;
              return NULL;
    }
    else{ // Checked before that we are at the root so no worry on this side
        int numbLeaves = this->GetNumberAllLeaves();
        std::string line;
//        int numbchild = this->GetNumberChildren();
        vector<TreeEM*> NodesVector=this->GetGeneralClassesVector();
        int numbchild=NodesVector.size();
        if(this->GetPriorsNodeVector().size()==0){
            numbchild=0;
        }
//        if(this->FlagOutliers>0){ // Case we are considering an outlier model and already checked beforehand that there are some pb
//            numbchild=this->GetPriorsNodeVector().size();
//            if(numbchild==0){ // means that no priors is put on the "normal classes"
//                numbchild=this->GetChild(0)->GetNumberChildren(); // In this case the number of anatomical classes must be at second level
//            }
//        }
//        else{
//            numbchild=this->GetNumberChildren();
//        }
        vector< vector<int> > NeighborhoodSummary;
        float * GValues = new float[3];
        GValues[0]=1;
        GValues[1]=0;
        GValues[2]=0;
        int numbClasses=0;
        while(getline(text,line)){
            istringstream in(line);
            std:: string type;
            in >> type;
            if (type == "GClasses"){
//                int numbClasses=0;
                in >> numbClasses;
                if(numbClasses!=numbchild && this->GetFlagOutliers()==0){
                    cout<< "Incompatibility of general classes number with Gfile"<<endl;
                    return NULL;
                }
            }
            if(type == "Neighborhood"){
                int CurrentClass=0;
                int NeighborhoodClass=0;
                vector<int> NeighboursCount;
                in >> CurrentClass;
                if(CurrentClass >= numbchild){
                    cout<< "Impossible class neighborhood offered"<<endl;
                    return NULL;
                }
                while(in >> NeighborhoodClass){
//                    in >> NeighborhoodClass;
                    if(NeighborhoodClass >= numbchild){
                        cout<< "Impossible class neighborhood offered"<<endl;
                        return NULL;
                    }
                    NeighboursCount.push_back(NeighborhoodClass);
                }
                NeighborhoodSummary.push_back(NeighboursCount);
            }
            if(type == "GValue"){
                int Level=0;
                float Value=0;
                in >> Level;
                in >> Value;
                // Checking if Value between 0 and 1 and thresholding if necessary
//                Value=Value>1?1:Value;
                Value=Value<0?0:Value;
                if(Level<=2){
                GValues[Level]=Value;
                }
            }
        }
        //
//        cout<<"GValues are ";
        for(int i=0;i<=2;i++){
//            cout<<GValues[i]<<" "; // At some point need to check that values ordered in GMatrix
        }
//        cout<<endl;
        // Checking if obtained informations coherent with preexistent data
        if(NeighborhoodSummary.size()!=numbClasses){
            cout<< "Incompatibility between number of classes and neighborhood information"<<endl;
            return NULL;
        }

        // Construction of GMatrix
        float * GMatrixResult = new float[numbLeaves * numbLeaves];
        for (int i=0;i<numbLeaves*numbLeaves;i++){
            GMatrixResult[i]=GValues[0];
        }
        if(numbchild==0){ // case where numbchild != numbClasses in particular when there are no priors and we cannot say anything about neighborhood properties
            for(int l1=0;l1<numbLeaves;l1++){
                for(int l2=0;l2<numbLeaves;l2++){
                    if(l1!=l2){
                        GMatrixResult[l1*numbLeaves+l2]=GValues[2];
                    }
                }
            }
            if (GValues!=NULL){
                delete [] GValues;
                GValues=NULL;
            }
            for(int c=0;c<numbchild;c++){
                NeighborhoodSummary[c].clear();
            }
            return GMatrixResult;
        }
        else{
        int * NumberLeavesPerChild = new int[numbClasses];
        int * IndexChangeChild=new int [numbClasses+1];
        for(int i=0;i<numbClasses+1;i++){
            IndexChangeChild[i]=0;
        }
//        vector<TreeEM*> NodesVector; // give the vector of nodes vector not counting the outlier node if exists
//        if(this->FlagOutliers==1){
//            NodesVector=this->GetPriorsNodeVector();
//            if (NodesVector.size()==0){ // case no priors are used
//                NodesVector=this->GetChild(0)->GetChildren();
//            }
//        }
//        else if(this->FlagOutliers==0){
//               NodesVector=this->GetChildren();
//            }
        // Once the node vector defined, obtention of the index changes from child to child.
        for(int c=0;c<numbClasses;c++){

            NumberLeavesPerChild[c]=NodesVector[c]->GetNumberAllLeaves();
            // If no Leaves detected means that the general class is no mixture but a leaf in itself
            NumberLeavesPerChild[c]=NumberLeavesPerChild[c]==0?1:NumberLeavesPerChild[c];

            IndexChangeChild[c+1]=NumberLeavesPerChild[c]+IndexChangeChild[c];

        }
        // First filling blocks of second level
        for(int l1=0;l1<numbLeaves;l1++){
            for(int l2=0;l2<numbLeaves;l2++){
                // Determining at which general class each leaf corresponds
                int Child1=0;
                int Child2=0;
                int test1=0;
                int test2=0;
                for(int c=0;c<numbchild;c++){
                    test1=IndexChangeChild[c];
                    test2=IndexChangeChild[c+1];
                    if (l1>=test1 && l1 <test2){
                        Child1=c;
                        break;
                    }
                }
                if(l1>=test2 && Child1==0){ // case where no child found and l1 correspond to out of the normal children => outlier
                    Child1=numbClasses;
                }
                test1=0;
                test2=0;
                for(int c=0;c<numbchild;c++){
                    test1=IndexChangeChild[c];
                    test2=IndexChangeChild[c+1];
                    if (l2>=test1 && l2 <test2){
                        Child2=c;
                        break;
                    }
                }
                if(l2>=test2 && Child2==0){ // same as with l1 => outlier
                    Child2=numbClasses;
                }
                // 1 on the diagonal
                if(l2==l1){
                    GMatrixResult[l1+l2*numbLeaves]=0;
                }
                // if leaves considered belong to same general class
                else if(Child1==Child2){
                    GMatrixResult[l1+l2*numbLeaves]=GValues[2];
                }
                else if(Child1>=numbClasses || Child2>=numbClasses){ // case where one of the leaves is considered as outlier : can be neighbour to everything so GValue[2] chosen
                    GMatrixResult[l1+l2*numbLeaves]=GValues[2];
                }
                else if(this->GetPriorsVector().size()>0){
                // if leaves considered belong to neighbouring general class
                int NumberOfNeighboursChild1=NeighborhoodSummary[Child1].size();
                bool NeighC1C2=0;
                for(int c=0;c<NumberOfNeighboursChild1;c++){
                    if(NeighborhoodSummary[Child1][c]==Child2){
                        NeighC1C2=1;
                        break;
                    }
                }
                 if(NeighC1C2){
                    GMatrixResult[l1+l2*numbLeaves]=GValues[1];
                }
                }

            }
        }
        //cleaning memory for GValues
        if (GValues!=NULL){
            delete [] GValues;
            GValues=NULL;
        }
        for(int c=0;c<numbchild;c++){
            NeighborhoodSummary[c].clear();
        }
        if(NumberLeavesPerChild !=NULL){
            delete [] NumberLeavesPerChild;
            NumberLeavesPerChild=NULL;
        }
        if(IndexChangeChild !=NULL){
            delete [] IndexChangeChild;
            IndexChangeChild=NULL;
        }
        NeighborhoodSummary.clear();
        return GMatrixResult;
    }
    }
}

void TreeEM::SetGMatrix(float * GMatrixToSet){
    // Gmatrix only set not NULL at Root
    if(!this->IsRoot()){
        cout<<"GMatrix can only be set at root"<<endl;
        if(this->GetGMatrixDirect()!=NULL){
            delete [] this->GetGMatrixDirect();
            this->GMatrix=NULL;
        }
        return;
    }
    // delete current GMatrix if there is one existing
    if(this->GetGMatrixDirect()!=NULL){
//        delete [] this->GetGMatrixDirect();
        delete [] this->GMatrix;
        this->GMatrix=NULL;
    }
    // allocate memory for new GMatrix
    int numbLeaves=this->GetNumberAllLeaves();
    this->GMatrix=new float [numbLeaves*numbLeaves];
    for(int l=0;l<numbLeaves*numbLeaves;l++){
        this->GMatrix[l]=GMatrixToSet[l];
    }
    return;
}

float * TreeEM::GetGMatrix(){
    if(!this->IsRoot()){
        return this->GetParent()->GetGMatrix();
    }
    else{
        return this->GMatrix;
    }
}

float * TreeEM::GetGMatrixDirect(){
    return this->GMatrix;
}

float * TreeEM::CopyGMatrix(){
    float * GMatrixToCopy=this->GetGMatrixDirect();
    if(GMatrixToCopy==NULL){
        return NULL;
    }

        int numbLeaves=this->FindRoot()->GetNumberAllLeaves();
        float * CopiedGMatrix=new float[numbLeaves * numbLeaves];
        for(int i=0;i<numbLeaves*numbLeaves;i++){
            CopiedGMatrix[i]=GMatrixToCopy[i];
        }

    return CopiedGMatrix;
}

void TreeEM::ClearMRFModel(){
    if(this->GetMRF()!=NULL){
        if(!this->IsLeaf()){
            delete [] this->GetMRF();
            this->SetMRF(NULL);
        }
    }
    // recursive part
    int numbchild=this->GetNumberChildren();
    for(int c=0;c<numbchild;c++){
        this->GetChild(c)->ClearMRFModel();
    }
}

/*** METHODS FOR OPTIMIZATION OF G MATRIX IN MRF ****/

// Return as an array of integers the result of the hard segmentation on the leaves
int * TreeEM::MakeHardSegLeaves(){
    int numbLeaves=this->GetNumberAllLeaves();
    int numelmasked=this->GetNumberMaskedElements();
    vector<TreeEM*> LeavesVector=this->GetAllLeaves();
    vector<float *> NormRespVector;
    for(int l=0;l<numbLeaves;l++){
        NormRespVector.push_back(LeavesVector[l]->GetNormResp());
    }
    int * HardSegLeavesResult=new int[numelmasked];
    for(int i=0;i<numelmasked;i++){
        HardSegLeavesResult[i]=0;
        float maxNormResp=-1;
        int maxClass=0;
        for(int l=0;l<numbLeaves;l++){
            if(NormRespVector[l][i]>maxNormResp){
                maxNormResp=NormRespVector[l][i];
                maxClass=l;
            }
        }
        HardSegLeavesResult[i]=maxClass;
    }
    return HardSegLeavesResult;
}

// Return the HardSeg on the leaves shifted according to dimension or direction Size is these of image (numel)
int * TreeEM::HardSegLeavesShifted(int dim, int dir){
    // First make sure that the input values are ok or change them accordingly
    if(dir!=0){
    dir=dir>0?1:-1;
    }
    dim=dim>2?2:dim;
    dim=dim<0?0:dim;

    // Initialisation of the needed values
    int Dimensions[3];
    Dimensions[0]=this->GetDataImage()->nx;
    Dimensions[1]=this->GetDataImage()->ny;
    Dimensions[2]=this->GetDataImage()->nz;
    int Shift[3];
    Shift[0]=1;
    Shift[1]=this->GetDataImage()->nx;
    Shift[2]=this->GetDataImage()->nx*this->GetDataImage()->ny;
    int numbx=this->GetDataImage()->nx;
    int numby=this->GetDataImage()->ny;
    int numbz=this->GetDataImage()->nz;
    int indexShift=Shift[dim]*dir;
    int numel=this->GetNumberElements();
    int newIndex=0;
    int oldIndex=0;
    // Allocation of memory for shifted matrix
    int * HardSegLeavesShifted= new int[numel];//{0};
    int * L2S_PTR=this->GetL2S();
    int * HardSegLeavesPointedNC=this->MakeHardSegLeaves();
    int * HardSegLeavesPointed = new int[numel];
        L2S_PTR=this->GetL2S();
        int * HardSegLeavesPointedNC_PTR = HardSegLeavesPointedNC;
        for (int i=0; i<numel; i++,L2S_PTR++) {

                HardSegLeavesPointed[i]=0;
                if(*L2S_PTR>=0){
                    HardSegLeavesPointed[i]=*HardSegLeavesPointedNC_PTR;
                    HardSegLeavesPointedNC_PTR++;
                }
                HardSegLeavesShifted[i]=0;
        }

        for (int i=0; i<numbx; i++) {
            for (int j=0; j<numby; j++) {
                for (int k=0; k<numbz; k++) {
                    newIndex=i+j*Shift[1]+k*Shift[2]+indexShift;
                    oldIndex=i+j*Shift[1]+k*Shift[2];
                    if (newIndex>0 && newIndex<numel) {
                        HardSegLeavesShifted[newIndex]=HardSegLeavesPointed[oldIndex];
                    }
                }
            }
        }
// Clearing memory before returning result
//    delete [] Dimensions;
//    Dimensions=NULL;
//    delete [] Shift;
//    Shift = NULL;
    delete [] HardSegLeavesPointed;
    HardSegLeavesPointed = NULL;
    delete [] HardSegLeavesPointedNC;
    HardSegLeavesPointedNC=NULL;
    return HardSegLeavesShifted;
}

// Summation of the result of the HardSeg on the 6 next neighbours of size numelmasked * numbLeaves (gives the gi vector)
int * TreeEM::SumHardSegNeighboursLeaves(){
    int numbLeaves=this->GetNumberAllLeaves();
    int numelmasked=this->GetNumberMaskedElements();
    int * SumHardSegLeavesNeighbours=new int[numelmasked * numbLeaves];
    for(int i=0;i<numelmasked*numbLeaves;i++){
        SumHardSegLeavesNeighbours[i]=0;
    }

    for(int dim =0;dim<3;dim++){
        for (int dir =-1;dir<=1;dir=dir+2){
            int * tmpHardSegLeaves=this->HardSegLeavesShifted(dim,dir);
            for(int i=0;i<numelmasked;i++){
                SumHardSegLeavesNeighbours[i+tmpHardSegLeaves[i]*numelmasked]++;
            }
            delete [] tmpHardSegLeaves;
            tmpHardSegLeaves = NULL;
        }

    }
    return SumHardSegLeavesNeighbours;
}

float * TreeEM::MakeHistogramHardSegLeaves(){
    int numbLeaves=this->GetNumberAllLeaves();
    int numelmasked=this->GetNumberMaskedElements();
    int * HardSegLeaves=this->MakeHardSegLeaves();
    int * SumHardSegLeavesNeighbours=this->SumHardSegNeighboursLeaves();
    int SizeHist=pow_int(numbLeaves,7);
    float * HistogramFinal=new float [SizeHist];
    int * Shift = new int[7];
    for(int i=0;i<7;i++){
        Shift[i]=pow_int(numbLeaves,i);
    }
    for(int i=0;i<numelmasked;i++){
        int Index=HardSegLeaves[i];
        for(int l=0;l<numbLeaves;l++){
            Index+=l*Shift[SumHardSegLeavesNeighbours[i]+1];
        }
        HistogramFinal[Index]++;
    }
    for(int i=0;i<pow_int(numbLeaves,7);i++){
        HistogramFinal[i]/=(float)numelmasked;
    }
    return HistogramFinal;
}

// Using current Softseg, fill the histogram of the possible configurations for each voxel and its 6 neighbours. The only non zero elements are those such that index correspond to sum of 6
float * TreeEM::HistogramSoftSegMRF(){
    int numbLeaves=this->GetNumberAllLeaves();
//    bool test=this->IsNormRespValidLeaves();
    // Obtention of the possible combinations for the segmentation
    int * TabSize = new int[7];
    for(int i=0;i<7;i++){
        TabSize[i]=numbLeaves;
    }
    int * Combination=CombinationBis(TabSize,7);
//    int * Shift = new int[7];
//    for(int i=0;i<7;i++){
//        Shift[i]=pow_int(numbLeaves,i);
//    }
    int * Shift = new int [numbLeaves+1];
    Shift[0]=1;
    Shift[1]=numbLeaves;
    for(int i=2;i<=numbLeaves;i++){
        Shift[i]=Shift[i-1]*7;
    }

    int numelmasked=this->GetNumberMaskedElements();
    // Storage of the shifted NormResp for all leaves
    vector<TreeEM*> LeavesVector=this->GetAllLeaves();
    vector<float *> NormRespLeaves;
    for(int l=0;l<numbLeaves;l++){
        for(int dim=0;dim<3;dim++){
            for(int dir =-1;dir<=1;dir=dir+2){
               NormRespLeaves.push_back(LeavesVector[l]->GetNormRespShifted(dim,dir));
            }
        }
    }

//    SaveTmpResult(NormRespLeaves[0],"/Users/Carole/Documents/PhD/NormRespLeavesShifted.nii.gz");
//    SaveTmpResultMasked(LeavesVector[0]->GetNormResp(),"/Users/Carole/Documents/PhD/NormRespNonShifted.nii.gz");

    // Allocation and initialisation of the Histogram
    int SizeHist=pow_int(7,numbLeaves)*numbLeaves;
    PrecisionTYPE * HistogramFinal = new PrecisionTYPE[SizeHist];
    for(int i=0;i<SizeHist;i++){
        HistogramFinal[i]=0;
    }
int numberCombi=pow_int(numbLeaves,7);
    // Convert Combination into gi vector
    int * GiVectorConverted = new int[numbLeaves*numberCombi];
    for(int i=0;i<numbLeaves*numberCombi;i++){
        GiVectorConverted[i]=0;
    }
//    int CountImpossibleCombination=0;
//    int CountPbIndexRecovery=0;

    for(int i=0;i<numberCombi;i++){
        int Ind=7*i;
        // First element in each session of combination corresponds to z and not g therefore omitted when reconstructing g vector
        // Ultimately in gi as a vector count of neighbours allocated to a special leaf for all leaves (Norm1 of each gi must be 6)
        for(int j=1;j<7;j++){
            GiVectorConverted[i*numbLeaves+Combination[Ind+j]]++;
        }
//        // Check if sum of GiVectorConverted is equal to 6
        int CheckSum=0;
        for(int l=0;l<numbLeaves;l++){
            CheckSum+=GiVectorConverted[i*numbLeaves+l];
        }


        // Find corresponding index in the histogram
//        int Index=Combination[i*numbLeaves];
        int Index=Combination[Ind];
        for(int l=0;l<numbLeaves;l++){
            Index+=GiVectorConverted[i*numbLeaves+l]*Shift[l+1];
//            Index+=l*Shift[GiVectorConverted[i*numbLeaves+l]+1];
        }


//        // Check if Index corresponds to proper sum
//        int tmp=Index;
//        int tmp2=0;
//        int sumTry=0;
//        for(int l=numbLeaves;l>0;l--){
//            tmp2=tmp/Shift[l];
//            sumTry+=tmp2;
//            tmp=tmp-tmp2*Shift[l];
//        }
//        if(sumTry!=6){
//            cout<< "Pb with index recovery at index "<<Ind<<" with value "<<sumTry<<endl;
//            CountPbIndexRecovery++;
//        }

        // Determine corresponding value to add
        int * S2L_PTR=this->GetS2L();

        for(int e=0;e<numelmasked;e++){
            PrecisionTYPE Value=LeavesVector[Combination[Ind]]->GetNormResp()[e]; // Value of norm resp corresponding to ze
            for(int j=1;j<7;j++){
                Value*=(PrecisionTYPE)NormRespLeaves[Combination[Ind+j]*6+j-1][S2L_PTR[e]];
//                // Check sum of normresp equal to 1 for all shifted leaves at all masked voxels
//                float testSumNormRespLeaves=0;
//                for(int l=0;l<numbLeaves;l++){
//                    testSumNormRespLeaves+=NormRespLeaves[l*6+j-1][S2L_PTR[e]];
//                }
//                if(testSumNormRespLeaves>1+10E-6){
//                    cout<<"SumNormRespLeaves shifted more than 1 for j "<<j<<" at index "<<e<<" with value "<<testSumNormRespLeaves<<endl;
//                }
//                if(testSumNormRespLeaves<1-10E-6){
//                    cout<<"SumNormRespLeaves shifted less than 1 for j "<<j<<" at index "<<e<<" with value "<<testSumNormRespLeaves<<endl;
//                }
//                if(Value==0){
////                    cout<< "Impossible combination at index "<<e<<endl;
//                    CountImpossibleCombination++;
//                }
            }
            if(Value/numelmasked>1){
                cout<<"Pb with value more than 1"<<endl;
            }
            HistogramFinal[Index]+=(PrecisionTYPE)Value;
        }
    }
    //Normalising Histogram
    float * HistogramReturn=new float[SizeHist];
    for(int i=0;i<SizeHist;i++){
        HistogramReturn[i]=HistogramFinal[i]/(float)numelmasked;
    }
    delete [] GiVectorConverted;
    GiVectorConverted=NULL;
    delete [] HistogramFinal;
    HistogramFinal=NULL;
    //Checking sum over Histogram Final
    PrecisionTYPE SumCheck=0;
    for(int i=0;i<SizeHist;i++){
        SumCheck+=(PrecisionTYPE)HistogramReturn[i];
        if (HistogramReturn[i]!=0){
            cout << "Histogram non zero at "<<i<<" with value "<<HistogramReturn[i]<<endl;
        }
    }
    for(int i=0;i<SizeHist;i++){
        HistogramReturn[i]/=SumCheck;
    }
    // Clearing memory
    for(int i=0;i<6*numbLeaves;i++){
        delete[] NormRespLeaves[i];
        NormRespLeaves[i]=NULL;
    }
    NormRespLeaves.clear();
    delete[] Shift;
    Shift=NULL;
    delete [] TabSize;
    TabSize=NULL;
    delete[] Combination;
    Combination=NULL;
    return HistogramReturn;
}

// return as a float array the log ratio described in equation 7 of VanLeemput 98
// When 0 value, put at -1E32 as default value. Will be discarded farther in the process
float * TreeEM::MRFOptMakeLogRatio(){
    int numbLeaves=this->GetNumberAllLeaves();
//    int SizeHist=pow_int(numbLeaves,7);
    int SizeHist=numbLeaves*pow_int(7,numbLeaves);
    float * HistogramFinal=this->HistogramSoftSegMRF();


    // check on non zero values of HistogramFinal
    int CountNonZero=0;
    float SumHistogram=0;
    for(int i=0;i<SizeHist;i++){
        if(HistogramFinal[i]>0){
            cout<<"Histogram non zero at "<< i << " with value "<<HistogramFinal[i]<<endl;
            CountNonZero++;
            SumHistogram+=HistogramFinal[i];
        }
    }

    // Initialisation of the result
    int SizeLogRatio=(int)(SizeHist*(numbLeaves-1))/2;
    float * LogRatio=new float[SizeLogRatio];

    for(int i=0;i<SizeLogRatio;i++){
        LogRatio[i]=0;
    }

    int SizeGCombi=SizeHist/numbLeaves;
    int CountPbLogRatio=0;
    for(int j=0;j<SizeGCombi;j++){
        int i1=0;
        int i2=0;

        for(int i=0;i<numbLeaves*(numbLeaves-1)/2;i++){
            i2++;
            if(i2==numbLeaves){
                i1++;
                i2=i1+1;
            }
                if(HistogramFinal[i2+j*numbLeaves]>1E-6 && HistogramFinal[i1+j*numbLeaves]>1E-6){
                LogRatio[i+j*(numbLeaves*(numbLeaves-1))/2]=logf(HistogramFinal[i1+j*numbLeaves]/HistogramFinal[i2+j*numbLeaves]);
                }
                else{
                    CountPbLogRatio++;
                    LogRatio[i+j*(numbLeaves*(numbLeaves-1))/2]=-10E32;
                }
            }
        }
    //Clearing memory before returning result
    delete [] HistogramFinal;
    HistogramFinal=NULL;
    cout<<"LogRatio calculated"<<endl;
    return LogRatio;

}

// return the integer array containing all the differences of int vectors obtained as possible configurations of one voxel and its 6 neighbours
int * TreeEM::MRFOptAMatrix(){
    int numbLeaves=this->GetNumberAllLeaves();
//    int SizeHist=pow_int(numbLeaves,7);
    int SizeHist=numbLeaves*pow_int(7,numbLeaves);
    int SizeAMatrix=(SizeHist*(numbLeaves-1)*pow_int(numbLeaves,2))/2;
//    int SizeAMatrix=pow_int(numbLeaves,6)*(numbLeaves)*(numbLeaves-1)/2*pow_int(numbLeaves,2);
    int * AMatrixResult=new int[SizeAMatrix];
    for(int i=0;i<SizeAMatrix;i++){
        AMatrixResult[i]=0;
    }

//    int * ShiftGOnly=new int[6];
//    for(int i=0;i<6;i++){
//        ShiftGOnly[i]=pow_int(numbLeaves,i);
//    }

    int * ShiftGOnly=new int[numbLeaves];
    for(int i=0;i<numbLeaves;i++){
        ShiftGOnly[i]=pow_int(7,i);
    }
    int SizeGCombi=SizeHist/numbLeaves;
    int * VMatrix=new int[SizeHist * numbLeaves*numbLeaves];
    for(int i=0;i<SizeHist*numbLeaves*numbLeaves;i++){
        VMatrix[i]=0;
    }
    int CountPbAMatrix=0;
    for(int j=0;j<SizeGCombi;j++){
        int i1=0;
        int i2=0;
//        int Index[numbLeaves];
        int * Index=new int[numbLeaves];
        int tmp=j;
        for(int l=numbLeaves-1;l>=0;l--){
            Index[l]=tmp/ShiftGOnly[l];
            tmp=tmp-Index[l]*ShiftGOnly[l];
        }
        // Conversion of Index into gi Vector
//        int Gi[numbLeaves];
        int * Gi=new int[numbLeaves];
        for(int l=0;l<numbLeaves;l++){
            Gi[l]=0;
        }
int CheckGiSum=0;
        for(int l=0;l<numbLeaves;l++){
//            Gi[Index[l]]++;
            Gi[l]=Index[l];
            CheckGiSum+=Gi[l];
        }

        for(int l=0;l<numbLeaves;l++){
//            int Zi[numbLeaves];
            int * Zi=new int[numbLeaves];
            for(int kl=0;kl<numbLeaves;kl++){
                Zi[kl]=0;
            }
            Zi[l]=1;

            for(int l1=0;l1<numbLeaves;l1++){
                for(int l2=0;l2<numbLeaves;l2++){
                    VMatrix[j*numbLeaves*numbLeaves*numbLeaves+l*numbLeaves*numbLeaves+l1*numbLeaves+l2]=Zi[l1]*(VMatrix[j*numbLeaves*numbLeaves*numbLeaves+l*numbLeaves*numbLeaves+l1*numbLeaves+l2]+Gi[l2]);
                }
            }
            if(Zi!=NULL){
                delete [] Zi;
                Zi=NULL;
            }
        }
        if(Gi!=NULL){
            delete [] Gi;
            Gi=NULL;
        }




//        for(int i=0;i<numbLeaves*(numbLeaves-1)/2;i++){
//            i2++;
//            if(i2==numbLeaves){
//                i1++;
//                i2=i1+1;
//            }
//            for(int l=0;l<numbLeaves;l++){
//                AMatrixResult[j*(numbLeaves)*(numbLeaves-1)/2+i1*numbLeaves+l]=AMatrixResult[j*(numbLeaves)*(numbLeaves-1)/2+i1*numbLeaves+l]+Gi[l];
//                AMatrixResult[j*(numbLeaves)*(numbLeaves-1)/2+i2*numbLeaves+l]=AMatrixResult[j*(numbLeaves)*(numbLeaves-1)/2+i2*numbLeaves+l]-Gi[l];
//            }
//        }
    }

//    // printing V result
//    cout<< "Printing VResult"<<endl;
//    int SizeLogRatio=SizeHist*(numbLeaves-1)/2;
//    int LengthMatrixG=pow_int(numbLeaves,2);
//    for(int i=0;i<SizeHist;i++){
//        for(int j=0;j<LengthMatrixG;j++){
//            cout<<VMatrix[i*LengthMatrixG+j]<<"  ";
//        }
//        cout<<endl;
//    }

    for(int j=0;j<SizeGCombi;j++){
        int i1=0;
        int i2=0;
                for(int i=0;i<numbLeaves*(numbLeaves-1)/2;i++){
                    i2++;
                    if(i2==numbLeaves){
                        i1++;
                        i2=i1+1;
                    }
                    for(int l=0;l<numbLeaves*numbLeaves;l++){
                        AMatrixResult[j*numbLeaves*numbLeaves*(numbLeaves)*(numbLeaves-1)/2+i*numbLeaves*numbLeaves+l]=AMatrixResult[j*numbLeaves*numbLeaves*(numbLeaves)*(numbLeaves-1)/2+i*numbLeaves*numbLeaves+l]+VMatrix[j*numbLeaves*numbLeaves*numbLeaves+i1*numbLeaves*numbLeaves+l];
                        AMatrixResult[j*numbLeaves*numbLeaves*(numbLeaves)*(numbLeaves-1)/2+i*numbLeaves*numbLeaves+l]=AMatrixResult[j*numbLeaves*numbLeaves*(numbLeaves)*(numbLeaves-1)/2+i*numbLeaves*numbLeaves+l]-VMatrix[j*numbLeaves*numbLeaves*numbLeaves+i2*numbLeaves*numbLeaves+l];
                    }
                }
    }
//    // printing A result
//    cout<< "Printing AResult"<<endl;
//    for(int i=0;i<SizeLogRatio;i++){
//        for(int j=0;j<LengthMatrixG;j++){
//            cout<<AMatrixResult[i*LengthMatrixG+j]<<"  ";
//        }
//        cout<<endl;
//    }

    delete [] ShiftGOnly;
    ShiftGOnly=NULL;
    delete [] VMatrix;
    VMatrix=NULL;
    cout<<"AMatrixCalculated"<<endl;
    return AMatrixResult;
}

// Solving the LS problem and returning the obtained result that will form the G matrix for the 6 neighbours
float * TreeEM::MRFOptSolveLS(){
    float * LogRatio=this->MRFOptMakeLogRatio();
    int * AMatrix=this->MRFOptAMatrix();
    int numbLeaves=this->GetNumberAllLeaves();
//    int SizeLogRatio=pow_int(numbLeaves,6)*(numbLeaves)*(numbLeaves-1)/2;
    int SizeLogRatio=pow_int(7,numbLeaves)*(numbLeaves)*(numbLeaves-1)/2;
    int SizeAMatrix=SizeLogRatio*pow_int(numbLeaves,2);
    // find indices not to use in the solution
    int * AuthorisedIndices=new int[SizeLogRatio];
    int SizeNonValid=0;
//    float minLogRatio=10E32;
    for(int i=0;i<SizeLogRatio;i++){
        AuthorisedIndices[i]=0;
//        if(LogRatio[i]<minLogRatio){
//            minLogRatio=LogRatio[i];
//        }
        if(LogRatio[i]<=-9E32){
            AuthorisedIndices[i]=-1;
            SizeNonValid++;
        }
    }

    // Obtain the reduced Matrices to use instead
    float * LogRatioReduced=new float[SizeLogRatio-SizeNonValid];
    int SizeLogRatioReduced=SizeLogRatio-SizeNonValid;
    int SizeAMatrixReduced=SizeAMatrix-SizeNonValid*pow_int(numbLeaves,2);
    float * AMatrixReduced=new float[SizeAMatrixReduced];
    int LengthMatrixG=pow_int(numbLeaves,2);
    int j=0;
//    int CountNanAMatrix=0;
    int CountNanLogRatio=0;
    for(int i=0;i<SizeLogRatio;i++){
        if(AuthorisedIndices[i]>=0){
            cout<< i<<"  ";
            LogRatioReduced[j]=LogRatio[i];
            if(LogRatioReduced[j]!=LogRatioReduced[j]){
                CountNanLogRatio++;
            }
            for(int l=0;l<LengthMatrixG;l++){
                AMatrixReduced[j*LengthMatrixG+l]=AMatrix[i*LengthMatrixG+l];
            }
            j++;
        }
    }
    cout<<endl;
    delete [] AuthorisedIndices;
    AuthorisedIndices=NULL;
    delete [] LogRatio;
    LogRatio=NULL;
    delete [] AMatrix;
    AMatrix=NULL;
    // Using SVD to solve LS linear problem


//    SVD SVDforMRF=SVD(AMatrixReduced,SizeLogRatio-SizeNonValid,LengthMatrixG);
//    float * VforMRF=SVDforMRF.getV();
//    float * UforMRF=SVDforMRF.getU();
//    float * SforMRF=SVDforMRF.getS();


//    int SizeSUT[4];
////    SizeSUT[0]=SizeLogRatio-SizeNonValid;
//    SizeSUT[0]=SizeSUT[1]=SizeSUT[2]=LengthMatrixG;
//    SizeSUT[3]=SizeLogRatio-SizeNonValid;
//    int un=(SizeLogRatio-SizeNonValid)>=LengthMatrixG?LengthMatrixG:(SizeLogRatio-SizeNonValid);
//    float * UT=TransposeMatrix(UforMRF,SizeLogRatio-SizeNonValid,un);


//    float * SUT=ProductMatrix(SforMRF,UT,SizeSUT);
//    delete[]UT;
//    UT=NULL;

//    int SizeVSUT[4];
//    SizeVSUT[0]=SizeVSUT[1]=SizeVSUT[2]=LengthMatrixG;
//    SizeVSUT[3]=SizeLogRatio-SizeNonValid;

//    float * VSUT=ProductMatrix(VforMRF,SUT,SizeVSUT);
//    delete [] SUT;
//    SUT=NULL;

//    int SizeFinal[4];
//    SizeFinal[0]=LengthMatrixG;
//    SizeFinal[1]=SizeFinal[2]=SizeLogRatio-SizeNonValid;
//    SizeFinal[3]=1;

//    float * FinalThetaResult=ProductMatrix(VSUT,LogRatioReduced,SizeFinal);
//    delete [] VSUT;
//    VSUT=NULL;


    float * AMatrixReducedT=TransposeMatrix(AMatrixReduced,SizeLogRatioReduced,LengthMatrixG);
//    // printing AT
//    for(int i=0;i<SizeLogRatioReduced;i++){
//        for(int j=0;j<LengthMatrixG;j++){
//            cout<<AMatrixReducedT[j*SizeLogRatioReduced+i]<<"  ";
//        }
//        cout<<endl;
//    }
    int SizeATA[4]={LengthMatrixG,SizeLogRatioReduced,SizeLogRatioReduced,LengthMatrixG};
    float * ATA=ProductMatrix(AMatrixReducedT,AMatrixReduced,SizeATA);
//    //printing ATA
//    for (int i=0;i<LengthMatrixG;i++){
//        for(int j=0;j<LengthMatrixG;j++){
//            cout<<ATA[i*LengthMatrixG+j]<<"  ";
//        }
//        cout<<endl;
//    }

    int CountNanATA=0;
    for(int i=0;i<LengthMatrixG*LengthMatrixG;i++){
        if(ATA[i]!=ATA[i]){
            CountNanATA++;
        }
    }
    cout<<CountNanATA<<" nan in the calculation of ATA"<<endl;
    invertMatrix(ATA,LengthMatrixG);
    int CountNanInvertATA=0;
    for(int i=0;i<LengthMatrixG*LengthMatrixG;i++){
        if(ATA[i]!=ATA[i]){
            CountNanInvertATA++;
        }
    }
//    cout<<"Printing inv ATA"<<endl;
//    for (int i=0;i<LengthMatrixG;i++){
//        for(int j=0;j<LengthMatrixG;j++){
//            cout<<ATA[i*LengthMatrixG+j]<<"  ";
//        }
//        cout<<endl;
//    }
    cout<<CountNanInvertATA<<" nan in the inversion of ATA"<<endl;
    int SizeATLR[4]={LengthMatrixG,SizeLogRatioReduced,SizeLogRatioReduced,1};
    float * ATLR=ProductMatrix(AMatrixReducedT,LogRatioReduced,SizeATLR);
    delete [] AMatrixReduced;
    AMatrixReduced=NULL;
    delete [] AMatrixReducedT;
    AMatrixReducedT=NULL;
    int SizeinvATAATLR[4]={LengthMatrixG,LengthMatrixG,LengthMatrixG,1};
    float * FinalResultBis=ProductMatrix(ATA,ATLR,SizeinvATAATLR);

    delete [] AMatrixReduced;
    AMatrixReduced=NULL;
    delete [] LogRatioReduced;
    LogRatioReduced=NULL;

    return FinalResultBis;
}

void TreeEM::SaveTmpResultMasked(float * ResultToSave,string filename){
            nifti_image * Result=nifti_copy_nim_info(this->GetDataImage());
            Result->data = (void *) calloc(Result->nvox, sizeof(float));
            float * Result_PTRtmp=static_cast<float *>(Result->data);
            float * ToSave_PTR=ResultToSave;
            int * L2S_PTR =this->GetL2S();
            int numel=this->GetNumberElements();
            for (int i=0; i<numel; i++, Result_PTRtmp++,L2S_PTR++) {
                *Result_PTRtmp=0;
                if(*L2S_PTR>=0){
                *Result_PTRtmp=*ToSave_PTR;
                    ToSave_PTR++;
                }
            }
            nifti_set_filenames(Result, filename.c_str(), 0, 0);
            nifti_image_write(Result);
            nifti_image_free(Result);
}

void TreeEM::SaveTmpResult(float * ResultToSave,string filename){
            nifti_image * Result=nifti_copy_nim_info(this->GetDataImage());
            Result->data = (void *) calloc(Result->nvox, sizeof(float));
            float * Result_PTRtmp=static_cast<float *>(Result->data);
            float * ToSave_PTR=ResultToSave;
            int numel=this->GetNumberElements();
            for (int i=0; i<numel; i++, Result_PTRtmp++,ToSave_PTR++) {
                *Result_PTRtmp=0;
                *Result_PTRtmp=*ToSave_PTR;

            }
            nifti_set_filenames(Result, filename.c_str(), 0, 0);
            nifti_image_write(Result);
            nifti_image_free(Result);
}

void TreeEM::SavePriorsAdapted(SEG_PARAMETERS* segment_param){
//    if(segment_param->AtlasWeight==0){
//        return;
//    }
    vector<float *> PriorsAdaptedVector=this->GetPriorsAdaptedVector();
    int numbGclasses=PriorsAdaptedVector.size();

    for(int c=0;c<numbGclasses;c++){
        nifti_image * Result=nifti_copy_nim_info(this->GetDataImage());
        Result->dim[4]=1;
        nifti_update_dims_from_array(Result);
        Result->data = (void *) calloc(Result->nvox, sizeof(float));
        int numel=Result->nx*Result->ny*Result->nz;
        float * PriorsAdapted_PTR=PriorsAdaptedVector[c];
        if(PriorsAdapted_PTR!=NULL){
        float * Result_PTRtmp=static_cast<float *>(Result->data);
        for(int i=0;i<numel;i++,Result_PTRtmp++,PriorsAdapted_PTR++){
            *Result_PTRtmp=*PriorsAdapted_PTR;
        }
        }
        nifti_set_filenames(Result, segment_param->filename_PriorsAdaptedOut[c].c_str(), 0, 0);
        nifti_image_write(Result);
        nifti_image_free(Result);
    }


}

void TreeEM::SaveMRFImage(SEG_PARAMETERS * segment_param){
    vector<TreeEM*> LeavesVector=this->GetAllLeaves();
    int numbLeaves=LeavesVector.size();
    if(LeavesVector[0]->GetMRF()==NULL){ // if no MRF, nothing to save, change segment_param so that nothing will appear in txt file
        segment_param->flag_MRFOut=0;
        return;
    }
    nifti_image * Result=nifti_copy_nim_info(this->GetDataImage());
    Result->dim[0]=4;
    Result->dim[4]=numbLeaves;
    nifti_update_dims_from_array(Result);
    Result->data = (void *) calloc(Result->nvox, sizeof(float));
    float *Result_PTR=static_cast<float*>(Result->data);
    int numel=this->GetNumberElements();
    for(int l=0;l<numbLeaves;l++){
        float *Result_PTRtmp=&Result_PTR[l*numel];
        int * L2S_PTR=this->GetL2S();
        float * MRFToSave_PTR=LeavesVector[l]->GetMRF();
        for(int i=0;i<numel;i++,Result_PTRtmp++,L2S_PTR++){
            *Result_PTRtmp=0;
            if(*L2S_PTR>=0){
                *Result_PTRtmp=*MRFToSave_PTR;
                MRFToSave_PTR++;
            }
        }
    }
    nifti_set_filenames(Result, segment_param->filename_MRFOut.c_str(), 0, 0);
    nifti_image_write(Result);
    nifti_image_free(Result);

}
















