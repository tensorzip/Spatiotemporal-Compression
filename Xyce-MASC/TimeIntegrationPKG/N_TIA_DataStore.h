//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
//   Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
//   NTESS, the U.S. Government retains certain rights in this software.
//
//   This file is part of the Xyce(TM) Parallel Electrical Simulator.
//
//   Xyce(TM) is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   Xyce(TM) is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with Xyce(TM).
//   If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//
// Purpose       : This file handles the class that defines the data arrays
//                 needed for the time integration algorithms.
//
// Special Notes :
//
// Creator       : Buddy Watts
//
// Creation Date : 6/1/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_DataStore_h
#define Xyce_N_TIA_DataStore_h

// ---------- Standard Includes ----------
#include <vector>

using std::vector;
// ----------   Xyce Includes   ----------
#include <N_LAS_fwd.h>
#include <N_TIA_fwd.h>
#include <N_LAS_EpetraMatrix.h>
#include <N_TIA_TwoLevelError.h>
#include <save_data.h>
#include <time.h>


//对于以下的代码，之后进行迁移
class bitmap{
public:
    int trim_m(vector<vector<double>>& m);
    bool compare_m(vector<vector<double>> A,vector<vector<double>> B);
    int get_m_size();
    int get_memo_cos();
    int init_bitmap(const vector<vector<double>> m);

    int construct_bitmap_0_simple();

    int construct_bitmap_0();
    int construct_bitmap_1();

    int convert_to_int_bitmap();
    int convert_bitmap_to_m(vector<vector<double>>& m);

    //xyce专有函数
    int convert_bitmap_to_Mptr(Xyce::Linear::Matrix* m);

    inline int fast_index(uint32_t bits,int type);
    int find_pos(int v, int &x,int &y);

    void set_compression_ration();
    vector<double> get_nzv();
    bool compare_bitmap(bitmap m);


    bitmap operator+(bitmap b_m);

private:
    vector<double> nzv;
    //用char来存储bitmap
    vector<vector<unsigned char> >bitmap_init;

    //对于小型矩阵，建立简单的bitmap_0
    vector<unsigned char> bitmap_0_simple;

    vector<uint32_t> bitmap_0;
    vector<uint32_t> bitmap_1;



    //用uint32_t来存储bitmap,初始化这两个vector则自动清空以上char类型的vector,以减少空间开销
//    vector<vector<uint32_t> >bitmap_0_int;
//    vector<vector<uint32_t>> bitmap_1_int;

    int compression_ratio0;
    int compression_ratio1;
    int m_size;
    int convert_sign;
};

//class save_data_for_sens{
//public:
//    bool compress_matrix();

//private:
//    std::vector<Xyce::Linear::Matrix*> Jac_v;
//};

namespace Xyce {
namespace TimeIntg {

//-----------------------------------------------------------------------------
// Class         : DataStore
// Purpose       : This is the class for defining data arrays needed in the
//                 time integration algorithms.
// Special Notes :
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
class DataStore
{
  public:
    DataStore(int max_order, const Linear::Builder &linear_builder);
    ~DataStore();

    void allocateSensitivityArrays(int nP, 
        bool includeTransientDirect,
        bool includeTransientAdjoint);

    void deleteSensitivityArrays();

    void allocateHBVectors();
    void allocateWaMPDEVectors();


    //此函数之后删除
    int bitmap_memo_cos(Xyce::Linear::Matrix* m_p);

  private:
    DataStore(const DataStore& rhs);
    DataStore &operator=(const DataStore& rhs);

  // DataStore Functions
  public:

    void updateSolDataArrays();
    bool updateStateDataArrays();
    void updateSolDataArraysAdjoint (int timeIndex);

    void setConstantHistory();
    void setConstantSensitivityHistory();
    void setZeroHistory();
    void setConstantHistoryAdjoint ();

    void setErrorWtVector(const TIAParams &tia_params, const std::vector<char> &     variable_type);
    double WRMS_errorNorm();

    bool equateTmpVectors ();
    bool usePreviousSolAsPredictor ();

    void stepLinearCombo ();

    double partialErrorNormSum();
    double partialQErrorNormSum();

    double partialSum_m1(int currentOrder);
    double partialSum_p1(int currentOrder, int maxOrder);
    double partialSum_q1();

    double delta_x_errorNorm_q1();
  
    bool getSolnVarData( const int & gid, std::vector<double> & varData );
    bool getStateVarData( const int & gid, std::vector<double> & varData );
    bool setSolnVarData( const int & gid, const std::vector<double> & varData );
    bool setStateVarData( const int & gid, const std::vector<double> & varData );
    bool getStoreVarData( const int & gid, std::vector<double> & varData );
    bool setStoreVarData( const int & gid, const std::vector<double> & varData );

    int getNumSolnVarData() const { return 11; }
    int getNumStateVarData() const { return 7; }
    int getNumStoreVarData() const { return 4; }
  
    bool setNextSolVectorPtr (Linear::Vector * solVecPtr);
    bool setNextSolVectorPtr (Linear::Vector & solVecPtr);// TT: added
    bool unsetNextSolVectorPtr ();
  
    bool resetAll(double absolute_error_tolerance, double relative_error_tolerance);

    bool resetFastTimeData ();
    //每次添加新矩阵时调用此函数
    bool construct_tensor(Xyce::Linear::Matrix* m,vector<Xyce::Linear::Matrix*> &tensor);
    bool decompress_tensor(vector<Xyce::Linear::Matrix*> &tensor);
    //位图存储
    bool construct_tensor_bits(Xyce::Linear::Matrix* m,vector<bitmap> &tensor);
    bool decompress_tensor_bits(vector<bitmap> &tensor);
    //判断两个矩阵是否相同（待改进）
    bool compare_matrix(Xyce::Linear::Matrix* A,Xyce::Linear::Matrix* B);

    // 评估压缩后的内存消耗，之后的版本会改进
    unsigned long evaluate_memo_cos(bool is_bitmap);
    //仿真实验
    void evaluate_memo_cos_virtual(Xyce::Linear::Matrix* m,bool is_bitmap);


    vector<vector<double>> convert_m_ptr(Xyce::Linear::Matrix* A);
  public:
    Xyce::Linear::Matrix* virtual_m;

    bitmap virtual_bitmap;

    unsigned long virtual_size;
    unsigned long virtual_size_0;
    //总步数
    int virtual_step;
    int virtual_num;

    std::vector<Xyce::Linear::Matrix*> Jac_v;

    std::vector<Xyce::Linear::Matrix*> dfdx_v;


    std::vector<Xyce::Linear::Matrix*> dqdx_v;

    std::vector<Linear::Vector *> rhs_v;

    std::vector<bitmap> Jac_bits;

    vector<bool> list;

    const Linear::Builder&        builder_;
 
    // limiter flag:
    bool limiterFlag;

    // TIA Arrays (pointers) for  Integration Solution Process:
    unsigned int maxOrder;
    unsigned int solutionSize;
    unsigned int stateSize;
    unsigned int storeSize;
    unsigned int leadCurrentSize;

    // temporary vectors:
    Linear::Vector * tmpSolVectorPtr;
    Linear::Vector * tmpStaVectorPtr;
    Linear::Vector * tmpStaDerivPtr;
    Linear::Vector * tmpStoVectorPtr;
    Linear::Vector * tmpLeadCurrentVectorPtr;
    Linear::Vector * tmpLeadDeltaVPtr;

    Linear::Vector * tmpLeadCurrentQDerivVectorPtr;

    // Predictors
    Linear::Vector * xn0Ptr;

    // Solutions:
    Linear::Vector * currSolutionPtr;
    Linear::Vector * lastSolutionPtr;
    Linear::Vector * nextSolutionPtr;

    // Used for pointer switching, related to nextSolPtrSwitched_
    Linear::Vector * savedNextSolutionPtr;

    // States:
    Linear::Vector * currStatePtr;
    Linear::Vector * lastStatePtr;
    Linear::Vector * nextStatePtr;

    // Storage:
    Linear::Vector * currStorePtr;
    Linear::Vector * lastStorePtr;
    Linear::Vector * nextStorePtr;
    
    // Lead current and power vectors 
    Linear::Vector * currLeadCurrentPtr;
    Linear::Vector * nextLeadCurrentPtr;
    
    Linear::Vector * currLeadDeltaVPtr;
    Linear::Vector * nextLeadDeltaVPtr;
    
    // for lead current calculations.  F component is
    // held in the store vector, Q component is here
    Linear::Vector * currLeadCurrentQPtr;
    Linear::Vector * nextLeadCurrentQPtr;    

    // Number of sensitivity parameters
    int numParams;

    // sensitivity vectors:
    Linear::MultiVector* sensRHSPtrVector;

    // adjoint sparse storage experiment
    Linear::FilteredMultiVector * sparseSensRHSMV;

    Linear::MultiVector* nextDfdpPtrVector;

    Linear::MultiVector* currDqdpPtrVector;
    Linear::MultiVector* nextDqdpPtrVector;

    Linear::MultiVector* nextDbdpPtrVector;

    Linear::MultiVector* currDXdpPtrVector;
    Linear::MultiVector* nextDXdpPtrVector;

    // adjoint sparse storage experiment
    std::vector<std::vector<int> > masterIndexVector;
    std::vector<int> masterIndexVectorSize;

    // matvecs:
    Linear::MultiVector* currDQdxDXdpPtrVector;
    Linear::MultiVector* lastDQdxDXdpPtrVector;
    Linear::MultiVector* nextDQdxDXdpPtrVector;

    Linear::MultiVector* currDFdxDXdpPtrVector;
    Linear::MultiVector* lastDFdxDXdpPtrVector;
    Linear::MultiVector* nextDFdxDXdpPtrVector;

    // Adjoint sensitivity solutions:
    Linear::Vector* nextLambdaPtr;
    Linear::Vector* currLambdaPtr;
    Linear::Vector* lastLambdaPtr;

    Linear::Vector* nextDQdxLambdaPtr;
    Linear::Vector* currDQdxLambdaPtr;
    Linear::Vector* lastDQdxLambdaPtr;

    Linear::Vector* nextDFdxLambdaPtr;
    Linear::Vector* currDFdxLambdaPtr;
    Linear::Vector* lastDFdxLambdaPtr;

    // Derivatives of States:
    Linear::Vector * currStateDerivPtr;
    Linear::Vector * nextStateDerivPtr;

    // Derivatives of Store for lead curent calculations
    Linear::Vector * currLeadCurrentQDerivPtr;
    Linear::Vector * nextLeadCurrentQDerivPtr;

    // Derivatives of dq/dp for sensitivity calculations
    Linear::MultiVector* nextDqdpDerivPtrVector;

    // Error Vectors
    Linear::Vector * errWtVecPtr;

    // Jacobian and RHS (pointers to objects in linear system)
    Linear::Matrix * JMatrixPtr;
    Linear::Vector * RHSVectorPtr;

    // NonLinear Solution Vectors
    Linear::Vector * newtonCorrectionPtr;
    Linear::Vector * qNewtonCorrectionPtr;

    // Mask for error norms (to allow some equations not to take part in
    // weighted norms)
    Linear::Vector * deviceErrorWeightMask_;

    // 2-level information:
    std::vector<TwoLevelError> innerErrorInfoVec;

    // new-DAE data (originally from the new-DAE derrived class)
    // Error Vectors
    Linear::Vector * qErrWtVecPtr;

    // DAE formulation vectors
    Linear::Vector * daeQVectorPtr;
    Linear::Vector * daeFVectorPtr;
    Linear::Vector * daeBVectorPtr;

    // DAE formulation matrices
    Linear::Matrix * dQdxMatrixPtr;
    Linear::Matrix * dFdxMatrixPtr;

    int get_dQdx_size();
    int get_dFdx_size();

    // HB temporary Matvec storage vectors
    Linear::Vector * dQdxVecVectorPtr;
    Linear::Vector * dFdxVecVectorPtr;

    // voltage limiting vectors
    Linear::Vector * dFdxdVpVectorPtr;
    Linear::Vector * dQdxdVpVectorPtr;

    // History arrays
    std::vector<Linear::Vector*> xHistory;
    std::vector<Linear::Vector*> qHistory;
    std::vector<Linear::Vector*> sHistory;    // state history
    std::vector<Linear::Vector*> stoHistory;  // store history
    std::vector<Linear::Vector*> leadCurrentHistory;  // history for lead current Q component.
    std::vector<Linear::Vector*> leadCurrentQHistory;  // history for lead current Q component.
    std::vector<Linear::Vector*> leadCurrentQDerivHistory;  // history for lead current dQ/dt component.
    std::vector<Linear::Vector*> leadDeltaVHistory;  // history for junction voltage
    
    // sensitivity histories
    std::vector< Linear::MultiVector* > dbdpHistory;
    std::vector< Linear::MultiVector* > dfdpHistory;
    std::vector< Linear::MultiVector* > dqdpHistory;

    // histories used for transient adjoints:
    int itAdjointIndex;
    bool adjointDcop;
    std::vector< int > orderHistory;
    std::vector< double > dtHistory;
    std::vector< double > timeHistory;
    std::vector< Linear::Vector *> solutionHistory;
    std::vector< Linear::Vector *> stateHistory;
    std::vector< Linear::Vector *> storeHistory;

   //temporary function for detect size

    int get_solutionHistory_size();

    // outer loop over parameter list, inner is the history
    std::vector< Linear::MultiVector * > functionSensitivityHistory;

    // adjoint sparse storage experiment
    std::vector< Linear::FilteredMultiVector* > sparseFunctionSensitivityHistory;

    long get_sparseFunctionSensitivity_size();

    // adjoint tmp matrix storage
    Linear::Matrix * tmpMatrixPtr;

    // Predictors
    Linear::Vector * qn0Ptr;

    // Step-size selection temporary vectors for two-level Newton
    Linear::Vector * delta_x;

    // Temporary vectors for WaMPDE interpolation
    Linear::Vector * tmpXn0APtr;
    Linear::Vector * tmpXn0BPtr;

    // These are for MPDE fast time scale points
    std::vector<double> timeSteps;
    std::vector<bool> timeStepsBreakpointFlag;
    std::vector<Linear::Vector*> fastTimeSolutionVec;
    std::vector<Linear::Vector*> fastTimeStateVec;
    std::vector<Linear::Vector*> fastTimeQVec;
    std::vector<Linear::Vector*> fastTimeStoreVec;

    std::vector<double> objectiveVec_;
    std::vector<double> dOdpVec_;
    std::vector<double> dOdpAdjVec_;
    std::vector<double> scaled_dOdpVec_;
    std::vector<double> scaled_dOdpAdjVec_;
    std::vector<double> paramOrigVals_; 

  private:
    bool nextSolPtrSwitched_;

    std::vector<int> indexIVars;
    std::vector<int> indexVVars;
    std::vector<int> indexMaskedVars;

    double absErrTol_, relErrTol_;
    double solsMaxValue;

    // Vectors for new LTE strategies
    Linear::Vector * maxSolutionPtr;
    Linear::Vector * relSolutionPtr;
    int index; 

    bool allocateSensitivityArraysComplete_;
    bool includeTransientAdjoint_;
    bool includeTransientDirect_;
};

} // namespace TimeIntg
} // namespace Xyce



#endif // Xyce_N_TIA_DataStore_h

