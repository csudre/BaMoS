#pragma once

//#ifndef _Seg_InitAndRun_h
//#define _Seg_InitAndRun_h
//#endif

//
//  TreeEM.h
//  TreeSeg
//
//  Created by Carole Sudre on 08/04/2013.
//  Copyright (c) 2013 Carole Sudre. All rights reserved.
//

//#ifndef TreeSeg_TreeEM_h
//#define TreeSeg_TreeEM_h

#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>

#include "nifti1.h"
#include "nifti1_io.h"
#include "_seg_matrix.h"
#include "_seg_tools.h"
#include "_seg_EM.h"
#include "_SVDCov.h"
#include "_EigenvaluesCov.h"



using namespace std;

class TreeEM;

#define PrecisionTYPE double
#define MaxNumbLeaves 30
#define MaxSupport 15
#define MaxNumbModal 5
#define MaxBForder 5
#define WeightThreshold 10E-6
#define MaxNumbBF ((MaxBForder+1)*(MaxBForder+2)*(MaxBForder+3))/6


typedef struct SEG_PARAMETERS{
    int  bias_order;
    PrecisionTYPE Bias_threshold;
    PrecisionTYPE ConvThreshold;
    PrecisionTYPE DPGauss_std;
    PrecisionTYPE AtlasWeight;
    PrecisionTYPE AtlasSmoothing;
    int KernelSize;
    float quantMin;
    float quantMax;
    float OutliersWeight;
    int OutliersMod;
    int  maxIteration;
    int maxIterationBF;
    int minIteration;
    int maxNumbLeaves;
    int numbDP;
    int numbCount;
    int choiceInitSplit;
    int uniformTypeChange;
//    int  verbose_level;
    int numb_classes;
    int numbmodal;
    int smoothing_order;
    PrecisionTYPE BWfactor;
    bool flag_Outliers;
    bool flag_manual_priors;
    bool flag_input;
    bool flag_mask;
    bool flag_out;
    bool flag_Bias;
    bool flag_bc_out;
    bool flag_bf_out;
    bool flag_data_out;
    bool flag_EMfirst_out;
    bool flag_BiASM;
    bool flag_DP;
    bool flag_DistClassInd;
    bool flag_Count;
    bool flag_Countmod;
    bool flag_intxt;
    bool flag_CEM;
    bool flag_CommonChanges;
    bool flag_MRF;
    bool flag_MRFPost;
    bool flag_MRFIn;
    bool flag_GMatrix;
    bool flag_GMatrixIn;
    bool flag_GMatrixPost;
    bool flag_PriorsAdaptedOut;
    bool flag_MRFOut;
    int SMOrder;
    int PriorsKept;
    vector<string> filename_input;
    vector<string> filename_out;
    vector<string> filename_DP;
    vector<string> filename_Count;
    string filename_datacorrected;
    string filename_intxt;
    string filename_correction;
    string filename_mask;
    string filename_datatxt;
    string filename_EMfirst;
    string filename_GMatrixPost;
    string filename_GMatrix;
    vector<string> filename_priors;
    vector<string> filename_PriorsAdaptedOut;
    string filename_MRFOut;

//    SEG_PARAMETERS(){
//        filename_datacorrected=NULL;
//        filename_correction=NULL;
//        filename_mask=NULL;
//    }

//    ~SEG_PARAMETERS(){
//        if(filename_datacorrected!=NULL){
//            delete[] filename_datacorrected;
//            filename_datacorrected=NULL;
//        }
//        if(filename_correction!=NULL){
//            delete[] filename_correction;
//            filename_correction=NULL;
//        }
//        if(filename_mask!=NULL){
//            delete[] filename_mask;
//            filename_mask=NULL;
//        }
//    }

  }SEG_PARAMETERS;

typedef struct Parameters{
    int DistributionType;
    int SizeParameters;
    float * ValueParameters;
    Parameters(){
        DistributionType=0;
        SizeParameters=0;
        ValueParameters=NULL;
    }
~Parameters(){
    if(ValueParameters!=NULL){
        delete [] ValueParameters;
        ValueParameters=NULL;
    }
}
} Parameters;

typedef struct GenKLD{
    float KLDTot;
//    float * KLD;
    GenKLD(){
        KLDTot=-1;
//        KLD=NULL;
    }
    ~GenKLD(){
        KLDTot=-1;
//        if(KLD!=NULL){
//            delete[] KLD;
//            KLD=NULL;
//        }
    }

    GenKLD* CopyGenKLD(){
        GenKLD * CopiedGenKLD=new GenKLD();
        CopiedGenKLD->KLDTot=this->KLDTot;
//        CopiedGenKLD->KLD=new float[MaxNumbModal];
//        if(this->KLD!=NULL){
//            for(int m=0;m<MaxNumbModal;m++){
//                CopiedGenKLD->KLD[m]=this->KLD[m];
//            }
//        }
//        else{
//            CopiedGenKLD->KLD=NULL;
//        }
        return CopiedGenKLD;
    }
}KLD;

typedef struct MergeKLD{
    TreeEM * Parent;
    int Child1;
    int Child2;
    float KLDTot;
//    float * KLD;

    MergeKLD(){
        Parent=NULL;
        Child1=-1;
        Child2=-1;
        KLDTot=-1;
//        KLD=NULL;
//        KLD=new float[MaxNumbModal];
//        for(int m=0;m<MaxNumbModal;m++){
//            KLD[m]=0;
//        }
    }

    ~MergeKLD(){
        Parent=NULL;
        Child1=-1;
        Child2=-1;
        KLDTot=-1;
//        if(KLD!=NULL){
//            delete [] KLD;
//            KLD=NULL;
//        }
    }

    MergeKLD * CopyMergeKLD(){
        MergeKLD * CopiedMergeKLD=new MergeKLD();
        CopiedMergeKLD->KLDTot=this->KLDTot;
        CopiedMergeKLD->Parent=this->Parent;
        CopiedMergeKLD->Child1=this->Child1;
        CopiedMergeKLD->Child2=this->Child2;
//        int numbmodal=CopiedSplitKLD->Parent->GetNumberModalities();
//        CopiedMergeKLD->KLD=new float[MaxNumbModal];
//        if(this->KLD!=NULL){
//        for(int m=0;m<MaxNumbModal;m++){
//            CopiedMergeKLD->KLD[m]=this->KLD[m];
//        }
//        }
//        else{
//            CopiedMergeKLD->KLD=NULL;
//        }
        return CopiedMergeKLD;
    }
}MergeKLD;

typedef struct SplitKLD{
    TreeEM * Parent;
    int ChildToSplit;
    int TypeChange; // 0 for change to Gaussian, 1 for classical split, 2 for change from 1 uniform to 2 gaussians
    float KLDTot;
    float * KLD;

    SplitKLD(){
        Parent=NULL;
        ChildToSplit=-1;
        KLDTot=-1;
        TypeChange=1;
//        KLD=NULL;
//        KLD=new float[MaxNumbModal];
//        for(int m=0;m<MaxNumbModal;m++){
//            KLD[m]=0;
//        }
    }

    ~SplitKLD(){
        Parent=NULL;
//        if(KLD!=NULL){
//            delete [] KLD;
//            KLD=NULL;
//        }
    }

    SplitKLD * CopySplitKLD(){
        SplitKLD * CopiedSplitKLD=new SplitKLD();
        CopiedSplitKLD->KLDTot=this->KLDTot;
        CopiedSplitKLD->Parent=this->Parent;
        CopiedSplitKLD->ChildToSplit=this->ChildToSplit;
        CopiedSplitKLD->TypeChange=this->TypeChange;
//        int numbmodal=CopiedSplitKLD->Parent->GetNumberModalities();
//        CopiedSplitKLD->KLD=new float[MaxNumbModal];
//        if(this->KLD!=NULL){
//        for(int m=0;m<MaxNumbModal;m++){
//            CopiedSplitKLD->KLD[m]=this->KLD[m];
//        }
//        }
//        else{
//            CopiedSplitKLD->KLD=NULL;
//        }
        return CopiedSplitKLD;
    }

} SplitKLD;

class TreeEM{
    TreeEM * Parent; // NULL if Root of the Tree
    vector<TreeEM *> Children; // Children of the tree : NULL if leaf
    float * DPChildren;// DP for number of subclasses per child only stored at root
    bool FlagDistClassInd;
    int FlagOutliers; // Only need to be stored at root.
    nifti_image * DataImage;
    nifti_image * Mask;
    nifti_image * Priors;
    float * PriorsAdapted;
    int * S2L;
    int * L2S;
//    PrecisionTYPE * NonNormResp; // Non normalized responsability
    float * NormResp; // Normalised responsability
//    PrecisionTYPE * Distribution;// Probabilistic distribution
    float NonNormWeight; // Non normalised weight
    float NormWeight; // Normalised weight
    Parameters * ParametersDistribution; // Parameters corresponding to the distribution ( non null pointer for the value only if leaf of the tree and therefore non mixture distribution)
    bool * SplitCheck;
    bool * MergeCheck;
    float * BFCoeffs;
    //PrecisionTYPE * BFBasisFunctions;
//    PrecisionTYPE * BFCorrection;
    float * DataBFCorrected;
    float IndFactor;
    int NumberMaskedElements;
    int * HardSeg;
    float * MRF;
    float * GMatrix;

    public :

    float * LogGaussBlur(float * GaussianBlur,int TotalSize);
    inline float pow_int(const float x,
                                    int exp);
    void invertMatrix(float * MatrixToInvert,int sizeMatrix);
    float determinant(float * data,int size);
    float * ProductMatrix(float * MatrixL, float * MatrixR, int * Size);
    float * TransposeMatrix(float * MatrixToTranspose, int SizeM,int SizeN);
    float Determinant_lib(float **a,int n);
    // Constructors, copy and destructor
    TreeEM();
    TreeEM(nifti_image * Data);

    TreeEM * CopyTree(TreeEM * ParentNew);
    nifti_image * CopyImage();
    nifti_image * CopyMask();
    nifti_image * CopyPriors();
    float * CopyPriorsAdapted();
    int * CopyS2L();
    int * CopyL2S();
//    PrecisionTYPE * CopyDistribution();
//    PrecisionTYPE * CopyNonNormResp();
    float * CopyNormResp();
    Parameters * CopyParameters();
    bool * CopyMergeCheck();
    bool * CopySplitCheck();
    float * CopyBFCorrection();
    float * CopyDataBFCorrected();
    float * CopyBFCoeffs();
    float CopyIndFactor();
    float * CopyDPChildren();
    int * CopyHardSeg();
    ~TreeEM();
    void MakeEmpty();

    // Get methods
    TreeEM * GetParent();
    vector<TreeEM * >GetChildren();
    int GetNumberChildren();
    int GetNumberGeneralClasses();
    TreeEM * GetChild(int c);
    nifti_image * GetDataDirect();
    nifti_image * GetDataImage();
    float GetMaxDataModal(int m);
    float GetMinDataModal(int m);
    float * GetDataModalMasked(int m);
    float * GetDataShifted(int dim,int dir);
    float * GetDataShiftedCorrected(int dim, int dir);
    float * GetDataCorrelation(int dim);
    float GetDataCorrelationMod(int dim, int* modalChoice);
    float MakeIndFactorMod(int dim);
    float MakeIndFactorTotMod();
    float MakeIndFactor(int dim);
    float MakeIndFactorTot();
    int GetNumberElements();
    int GetNumberModalities();
    int GetNumberMaskedElements();
    int MakeNumberMaskedElements();
//    PrecisionTYPE * GetDistribution();
    int GetNumbbins();
nifti_image * GetMask();
    nifti_image * GetMaskDirect();
    nifti_image * GetPriors();
    nifti_image * GetPriorsDirect();
    float * GetPriorsAdaptedDirect();
    float * GetPriorsAdapted();
int * GetS2L();
    int * GetS2LDirect();
int * GetL2S();
    int * GetL2SDirect();
    int * GetHardSeg();
    int * GetHardSegDirect();
    void MakeHardSeg();
    void MakeL2S();
    void MakeS2L();
//PrecisionTYPE * GetNonNormResp();
float * GetNormResp();
float  GetNonNormWeight();
float  GetNormWeight();
    float GetPartNormWeight();
   float GetPartNormWeightAbovePriors();
Parameters * GetParameters();
int GetDistributionType();
    float * GetMean();
    float * GetVariance();
    float * GetMeanDirect();
    float * GetDiagVarianceDirect();
    float * GetVarianceDirect();
    float * GetParametersValue();
    int GetSizeParameters();
    vector<nifti_image*> GetPriorsVector();
    vector<TreeEM*> GetPriorsNodeVector();
    vector<TreeEM*> GetGeneralClassesVector();
    vector<TreeEM*> GetMainNodesVector();
    TreeEM * GetNodeOutlier();
    vector<float *> GetPriorsAdaptedVector();
    void GetLeaves(vector<TreeEM*> & LeavesVector);
    vector<TreeEM *> GetAllLeaves();
    vector<TreeEM *> GetAllDirectLeaves();
    int GetNumberDirectLeaves();
    bool * GetSplitCheck();
    bool * GetMergeCheck();
    void GetNumberLeaves(int & numbleaves);
    int GetNumberAllLeaves();
    int GetNumberFreeParameters();
    float GetLogLikelihood();
    float GetLogLikelihoodCEM(int c);
    float * GetBFCoeffs();
    float * GetBFCoeffsDirect();
    //PrecisionTYPE * GetBFBasisFunctionsDirect();
    //PrecisionTYPE * GetBFBasisFunctions();
    float * GetBFCorrection();
    float * GetBFCorrectionDirect();
    float * GetDataBFCorrected();
    float * GetDataBFCorrectedDirect();
float GetIndFactor();
float * GetDPChildren();
float * GetDPChildrenDirect();
bool GetFlagDistClassInd();
int GetFlagOutliers();
    enum TreeComponent {
        ROOT = 0,
        BRANCH=1,
        LEAF = 2,
        INITIALNODE=3
    };
    enum PartialResultType {
        DISTRIBUTION=0,
        NONNORMRESP=1,
        NORMRESP=2
    };




    float * GetPartialResult(PartialResultType ResultType,SEG_PARAMETERS * segment_param);

    // Set methods
    void SetParent(TreeEM * Parent);
    void ReinitialiseTree();
    void SetChildren(vector<TreeEM *> Children);
    void AddChild(TreeEM * ChildAdded);
    void SetData(nifti_image * data);
    void NormaliseDataImage();
    void QuantilizeDataImage(SEG_PARAMETERS * segment_param);
    void ConvertDataImageToFloat();
//    void SetDistribution(float * Distribution);
    void SetMask(nifti_image * MaskImage);
    void SetNumberMaskedElements(int NumberMaskedElementsInput);
    void MakeMaskBinary();
    void SetPriors(nifti_image * Priors);
    void MakePriorsFloat();
    void MakePriorsProbabilityType();
    void SetPriorsAdapted(float * PriorsAdapted);
    void MakePriorsAdaptedProbabilityType();
    void SetS2L(int * S2L);
    void SetL2S(int * L2S);
//    void SetNonNormResp(PrecisionTYPE * NonNormResp);
    void SetNormResp(float * NormResp);
    void InitialiseNormRespWithPriors();
    void InitialiseNormResp();
    void SetNonNormWeight(float NonNormWeight);
    void SetNormWeight(float  NormWeight);
    void SetParameters(Parameters * ParametersToSet);
    void SetDistributionType(int DistributionType);
    void SetBFCoeffs(float * BFCoeffs);
    void SetBFCorrection(float * BFCorrection);
    void SetDataBFCorrected(float * DataBFCorrectedInput);
    void MakeParametersMixture();
    void SetIndFactor(float IndFactorInput);
    void SetDPChildren(float* DPChildrenInput);
    void SetFlagDistClassInd(bool flag_DistClassInd);
    void SetFlagOutliers(int flag_Outliers);
    void SetHardSeg(int * HardSegInput);


// EM related methods
    void UpdateParameters();
    void UpdateParametersCEM(int c);
    void UpdateHardSeg();
    void UpdateGaussianParameters();
    void UpdateGaussianParametersCEM(int c);
    void UpdateBFCoeffs();
    void UpdateBFCorrection();
    void UpdateDataBFCorrected();
    void UpdateGaussianMean();
    void UpdateGaussianMeanCEM(int c);
    void UpdateGaussianVariance();
    void UpdateGaussianVarianceCEM(int c);
    void UpdatePriorsAdapted(SEG_PARAMETERS * segment_param);
//    void UpdateDistribution();

//    void SumToDistributionWeightedDistributionOfChild(int c);

//    void ReputDistributionToZero();
    void UpdateGaussianDistribution();

    void UpdateNonNormResp(SEG_PARAMETERS * segment_param);
    void UpdateNonNormRespCEM();
    void UpdateNormRespRoot();

//    void ReputNonNormRespToZero();
//    void CalculatingNonNormResp();

//    void SumToNonNormRespWeightedNonNormRespOfChild(int c);

//    void MultiplyCurrentNonNormRespByNormWeight();

    void UpdateNormResp();
    void MakeNonNormWeightedSum(float * NonNormSum,SEG_PARAMETERS * segment_param);
    void MakeNonNormWeightedSumCEM(float * NonNormSum);
    void MakeWeightedDist(float * WeigtedDist);
    float * MakeUniformDistribution();
    float * MakeGaussianDistribution();
//    void NormaliseNonNormResp();

    void DivideNonNormRespByRootNonNormResp();

    void UpdateNonNormWeights();
    void UpdateNonNormWeightsCEM(int c);
    void UpdateNormWeights();
    float * NormalisedPriorsForCCL();
    float CompleteLogLikelihoodFinBis();
    float CompleteLogLikelihoodFin(SEG_PARAMETERS * segment_param);
    float CompleteLogLikelihood();
    void EMCompleteIteration(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration, SEG_PARAMETERS * segment_param);
    void CEMPartialIteration(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration, int c);
    void CEMCompleteIteration(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration,SEG_PARAMETERS * segment_param);
    void RunFullEM(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration, SEG_PARAMETERS * segment_param);
    void RunPartialCEM(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration, int c);
    void RunFullCEM(float & CompleteLogLikelihood, float & OldCompleteLogLikelihood, int & Iteration,SEG_PARAMETERS * segment_param);
    void ReplugParameters(TreeEM * TreeToReplug);
    TreeEM * RunFullBiASM(SEG_PARAMETERS * segment_param);
    TreeEM * RunFullBiASM_bis(SEG_PARAMETERS * segment_param);
    //Tree related methods
    void CreateAndAddChild();
    void CreateAndAddChildPriors(SEG_PARAMETERS * segment_param,nifti_image * Priors,int DistributionType);
    void CreateAndAddChildWeight(SEG_PARAMETERS *segment_param, float Weight, int DistributionType);
//int GetLevel();
    void CollapseOnlyChild();
    void CollapseOnlyChildTot();
    void CollapseChild(TreeEM * ChildToCollapse);
    void AllLeavesToChildrenLevel();
    int Index(TreeEM * ElementToFind);
    void FindAllAtCertainLevel(vector<TreeEM*> & ResultVector, int Level);
    void FindDepth( int & ResDepth);
    TreeEM * FindRoot();
    void DividePartNonNormResp(PrecisionTYPE * SumNonNormResp);
    void InitialiseBasicTree(SEG_PARAMETERS * segment_param);
    void ClearSMChecks();
    void PutAllLeavesToChildrenLevel();
    void ModifyNormWeightsForChildrenLevel();
    void ModifyNormWeightsForMainNodesChildrenLevel();
    void PutAllLeavesToMainNodesChildrenLevel();
    void RecPutAllLeavesToChildrenLevel();

// Tests on tree
bool IsEqual(TreeEM * TreeEqualityCheck);
bool IsLeaf();
    bool IsBranch();
    bool IsRoot();
void DoesBelongToCurrentTree(TreeEM * Element, bool & test);
int WhatTypeOfTreeComponent();
    bool IsTreeValid();
    bool AreWeightsValid();
    bool IsNormRespValid();
    bool AreNormRespValid();
    bool IsNormRespValidLeaves();
    bool IsNormRespValidGeneral();
    bool IsBasicTree();
    bool IsNumbChildOKWithDistType();
    bool IsDataImageNormalised();
    bool IsDataFloat();
    bool IsPriorsFloat();
    bool AreBFCoeffsDirectUseable();
    bool IsStructureSimilar(TreeEM * TreeTest);
    // Parameters structure

    int CalculateSizeParameters();
    int CalculateSizeParameters(int DistributionType);
    void CreateAllocateAndInitializeParameters(int DistributionType);
    bool CheckForValidityOfParametersStructure();
    bool CheckForSizeParametersValidity();

    bool ArePriorsNormalised();
    bool ArePriorsAdaptedNormalised();
    bool IsPriorsVectorFilledWithNULL();
    bool IsPriorsAdaptedVectorFilledWithNULL();
    bool IsOneOfThePriorsNULL();
    bool IsOneOfThePriorsAdaptedNULL();
    void NormalisePriors();
    void NormalisePriorsAdapted();
    bool IsPriorsProbabilityType();
    bool IsPriorsAdaptedProbabilityType();
    bool CheckForValidityOfParametersStructure(Parameters * ParametersToCheck);
    bool CheckForSizeParametersValidity( Parameters * ParametersToSet);
    bool IsMRFZero();

    /**************** METHODS FOR SM OPERATIONS *********************/
    // Functions needed for calculation of KLD

    vector<float *> GetDistHistogram();
    float * GetDistHistogramTotal();
    vector<float *> GetUniformDistributionHist();
    vector<float *> GetGaussianDistributionHist();
    float * GetUniformDistributionHistTotal();
    float * GetGaussianDistributionHistTotal();
    vector<float *> GetDataHistogram();
    float * GetDataHistogramTotal();

    float GetDistCompKLDTotal();
    float * GetDistCompKLD();

    float CriterionCalculation();
    float CriterionCalculationDPInd();
    float CriterionCalculationDPNonInd(int ChildModified);
    float CriterionCalculationSplitDP();
    float CriterionCalculationMergeDP();

    // Functions needed for merging operations
    GenKLD * GetGenKLD();
    MergeKLD * GetKLDMerge(int Child1,int Child2);
    vector<MergeKLD *> GetVectorKLDMergeChildren();
    bool IsMergeChecked();
    vector<MergeKLD *> GetVectorKLDMergeLeaves();
    MergeKLD * GetToMerge();
    void MergeOperation(int Child1,int Child2,SEG_PARAMETERS * segment_param);
    Parameters * ParametersForMerging(int Child1,int Child2);
    vector<MergeKLD*> GetMergeMoreVertical(int numbclasses);
    vector<int> OrderingMergingChildren();
    vector<MergeKLD*> OrderingMergingLeaves();
    vector<MergeKLD*> GetMergeOrderVertical();
    TreeEM * RunMergeOperation(MergeKLD * MergeTry, bool & AcceptanceDecision,SEG_PARAMETERS * segment_param);
    TreeEM * RunMergeOperation(vector<MergeKLD *> MergeTry, bool & AcceptanceDecision,SEG_PARAMETERS * segment_param);
    // Functions needed for splitting operations
    SplitKLD * GetKLDSplit(int ChildToSplit);
    vector<SplitKLD *> GetVectorKLDSplitChildren(SEG_PARAMETERS * segment_param);
    vector<SplitKLD *> GetVectorKLDSplitLeaves(SEG_PARAMETERS * segment_param);
    bool IsSplitChecked();
    SplitKLD * GetToSplit(SEG_PARAMETERS * segment_param);
    vector<int> OrderingSplittingChildren();
    vector<SplitKLD*> OrderingSplittingLeaves(SEG_PARAMETERS* segment_param);
    vector<SplitKLD*> GetSplitOrder();
    vector<SplitKLD*> GetSplitMoreVertical(int numbclasses,SEG_PARAMETERS* segment_param);
    int * Combination(int n,int k);
    int * CombinationBis(int * TabSize,int numbclasses);
    vector<SplitKLD*> GetToSplitMore(vector<int> ClassesToSplit,SEG_PARAMETERS * segment_param);
    vector<SplitKLD*> GetSplitOrderVertical(SEG_PARAMETERS * segment_param);
    inline int Factorial(int x);
    int NumbComb(int k, int n);
    void SplitOperation(SplitKLD * SplitTry, int choiceInit, SEG_PARAMETERS* segment_param);
    void SplitOperation(int ChildToSplit,int choiceInit,SEG_PARAMETERS* segment_param);
    Parameters ** ParametersDoubleForSplittingUniform(int choiceInit);
    Parameters ** ParametersForSplitting(int choiceInit);
    Parameters * ParametersForSplittingUniform3();
    float * MeanForSplitInitialisation(float MaxEigen, float * MaxVector);
    float * MeanForSplitInitialisationUniform(float MaxEigen, float* MaxVector);
    void SortingSVDEigen(float * SVDEigen, int * Index);
    float * VarianceForSplitInitialisation(float MaxEigen, float * MaxVector);
    float * VarianceForSplitInitialisationUniform(float MaxEigen, float * MaxVector);
    TreeEM * RunSplitOperation(SplitKLD * SplitTry, bool & AcceptanceDecision, SEG_PARAMETERS * segment_param);
    TreeEM * RunSplitOperation(vector<SplitKLD *> SplitTry, bool & AcceptanceDecision, SEG_PARAMETERS * segment_param);

    float * AdaptPriors(SEG_PARAMETERS * segment_param);
    nifti_image * BuildConstantPriors(float WeightInput);
    nifti_image * CreatePriorsFromAdaptedPriors();
    nifti_image * TransformNormRespIntoPriors(SEG_PARAMETERS * segment_param);
    // Methods for bias field correction
    float * MakeBasisFunctions();
    float * MakeWMatrixChildren();
    float * MakeWMatrixLeaves();
    float * MakeRMatrixChildren();
    float * MakeAtWAMatrixChildren();
    float * MakeAtWRVectorChildren();
    float * MakeInvAtWAMatrixChildren();
    float * MakeFinalBFCoeffsChildren();
    float * MakeBFCoeffsDirect();
    float * MakeBFCorrection();
    float * MakeDataBFCorrected();

    // Saving Results
    void SaveAllClasses(string filenameOut, SEG_PARAMETERS * segment_param);
    void SaveGeneralClasses(string filenameOut, SEG_PARAMETERS * segment_param);
    void SaveBFBasisFunctions(string filenameBF);
    void SaveBFBasisFunctions(int BFnumber, string filenameBF);
    void SaveBFCorrection(string filenameBF);
    void SaveBFCorrectedData(string filenameBF);
    void SavePriorsAdapted(SEG_PARAMETERS* segment_param, string filename);

    int GetLevel();
    int GetNumberLevels();
    vector<TreeEM*> GetAllTreesFromLevel(int l);
    int FindIndex(TreeEM * TreeToLook);
    int FindGeneralClass();
    TreeEM* FindGeneralClassPriors();
    int * GetHierarchy();
    TreeEM * FindMainNode();

    void SaveTreeInTextFile(string filenameOut,SEG_PARAMETERS * segment_param);
    void ModifyCountFiles(vector<string> CountFiles);

    float * GaussianBlurring(float * CountHistogram, float gauss_std, vector<int> dim);

    // MRF functions
    void ClearMRFModel();
    float * CopyGMatrix();
    float * CopyMRF();
    float * GetMRF();
    float * GetGMatrix();
    float * GetGMatrixDirect();
    void SetGMatrix(float * GMatrixToSet);
    void SetMRF(float * MRFToSet);
    void UpdateMRF();
    float * PrepareGMatrixFromFile(string GMatrixFilename);
    float * GMulSumNeighbNormRespExp();
    float * GetSumNeighboursNormResp();
    float * GetNormRespShifted(int dim, int dir);

    int * MakeHardSegLeaves();
    int * HardSegLeavesShifted(int dim,int dir);
    int * SumHardSegNeighboursLeaves();
    float * MakeHistogramHardSegLeaves();
    float * HistogramSoftSegMRF();
    float * MRFOptMakeLogRatio();
    float * MRFOptSolveLS();
    int * MRFOptAMatrix();

    void SaveTmpResultMasked(float * ResultToSave,string filename);
    void SaveTmpResult(float * ResultToSave,string filename);
    void SaveMRFImage(SEG_PARAMETERS * segment_param);
    void SavePriorsAdapted(SEG_PARAMETERS* segment_param);

    template <class NewTYPE>
    int seg_changeDatatype(nifti_image *image)
    {
        switch(image->datatype){
        case DT_BINARY:
            seg_changeDatatype1<NewTYPE,bool>(image);
            break;
        case NIFTI_TYPE_UINT8:
            seg_changeDatatype1<NewTYPE,unsigned char>(image);
            break;
        case NIFTI_TYPE_INT8:
            seg_changeDatatype1<NewTYPE,char>(image);
            break;
        case NIFTI_TYPE_UINT16:
            seg_changeDatatype1<NewTYPE,unsigned short>(image);
            break;
        case NIFTI_TYPE_INT16:
            seg_changeDatatype1<NewTYPE,short>(image);
            break;
        case NIFTI_TYPE_UINT32:
            seg_changeDatatype1<NewTYPE,unsigned int>(image);
            break;
        case NIFTI_TYPE_INT32:
            seg_changeDatatype1<NewTYPE,int>(image);
            break;
        case NIFTI_TYPE_FLOAT32:
            seg_changeDatatype1<NewTYPE,float>(image);
            break;
        case NIFTI_TYPE_FLOAT64:
            seg_changeDatatype1<NewTYPE,double>(image);
            break;
        default:
            fprintf(stderr,"[NiftyReg ERROR] seg_changeDatatype\tThe initial image data type (%d) is not supported\n",image->datatype);
            exit(1);
        }
        return 1;
    }

    template <class NewTYPE, class DTYPE>
    int seg_changeDatatype1(nifti_image *image)
    {
        // the initial array is saved and freeed
//        DTYPE *initialValue = (DTYPE *)malloc(image->nvox*sizeof(DTYPE));
        DTYPE * initialValue=new DTYPE[image->nvox];
        DTYPE *ImagePTR=static_cast<DTYPE *>(image->data);
//        memcpy(initialValue, image->data, image->nvox*sizeof(DTYPE));
        int numel=image->nvox;
//        int numel=image->nx*image->ny*image->nz;
        for(int i=0;i<numel;i++){
            initialValue[i]=ImagePTR[i];
        }

        // the new array is allocated and then filled
        if(sizeof(NewTYPE)==sizeof(unsigned char)) image->datatype = NIFTI_TYPE_UINT8;
        else if(sizeof(NewTYPE)==sizeof(float)) image->datatype = NIFTI_TYPE_FLOAT32;
        else if(sizeof(NewTYPE)==sizeof(double)) image->datatype = NIFTI_TYPE_FLOAT64;
        else{
            fprintf(stderr,"[NiftyReg ERROR] seg_changeDatatype\tOnly change to unsigned char, float or double are supported\n");
            exit(1);
        }
        free(image->data);
        image->nbyper = sizeof(NewTYPE);
        image->data = (void *)calloc(image->nvox,sizeof(NewTYPE));
        NewTYPE *dataPtr = static_cast<NewTYPE *>(image->data);
        if(sizeof(NewTYPE)==sizeof(unsigned char)){
            for(unsigned int i=0; i<image->nvox; i++) dataPtr[i] = (unsigned char)(round(initialValue[i]));
        }
        else{
            for(unsigned int i=0; i<image->nvox; i++) dataPtr[i] = (NewTYPE)(initialValue[i]);
        }
//        nifti_set_filenames(image,"/Users/Carole/Documents/PhD/TestImage.nii.gz",0,0);
//        nifti_image_write(image);

        free(initialValue);
        return 1;
    }

};



