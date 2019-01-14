#ifndef FTL_HITS_TREE
#define FTL_HITS_TREE

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeBase.h"

using namespace std;

//---Define the TTree branches
#define DYNAMIC_TREE_NAME FTLElectronsTree

#define DATA_TABLE                              \
  DATA(int, event)                              \
  DATA(int, lumi)                               \
  DATA(int, run)                                \
  DATA(float, genVtx_x)                         \
  DATA(float, genVtx_y)                         \
  DATA(float, genVtx_z)                         \
  DATA(float, genVtx_t)                   

#define DATA_CLASS_TABLE                                                \
  DATA(vector<int>,   electrons_idx)                                    \
  DATA(vector<float>, electrons_pt)                                     \
  DATA(vector<float>, electrons_eta)                                    \
  DATA(vector<float>, electrons_phi)                                    \
  DATA(vector<int>,   electrons_isEB)                                   \
  DATA(vector<int>,   electrons_isEE)                                   \
  DATA(vector<float>, electrons_fbrem)                                  \
  DATA(vector<float>, electrons_dEtaIn)                                 \
  DATA(vector<float>, electrons_dPhiIn)                                 \
  DATA(vector<float>, electrons_eop)                                    \
  DATA(vector<float>, electrons_ecalEnergy)                             \
  DATA(vector<float>, electrons_sigmaEtaEta)                            \
  DATA(vector<float>, electrons_sigmaIetaIeta)                          \
  DATA(vector<float>, electrons_sigmaIphiIphi)                          \
  DATA(vector<float>, electrons_e1x5)                                   \
  DATA(vector<float>, electrons_e2x5Max)                                \
  DATA(vector<float>, electrons_e5x5)                                   \
  DATA(vector<float>, electrons_r9)                                     \
  DATA(vector<float>, electrons_hcalOverEcal)                           \
  DATA(vector<float>, electrons_mva)                                    \
  DATA(vector<float>, electrons_eta_atBTL)                              \
  DATA(vector<float>, electrons_phi_atBTL)                              \
  DATA(vector<float>, electrons_local_x_atBTL)                          \
  DATA(vector<float>, electrons_local_y_atBTL)                          \
  DATA(vector<float>, electrons_local_z_atBTL)                          \
  DATA(vector<float>, electrons_global_R_atBTL)                         \
  DATA(vector<float>, electrons_eta_atETL)                              \
  DATA(vector<float>, electrons_phi_atETL)                              \
  DATA(vector<float>, electrons_local_x_atETL)                          \
  DATA(vector<float>, electrons_local_y_atETL)                          \
  DATA(vector<float>, electrons_local_z_atETL)                          \
  DATA(vector<float>, electrons_global_R_atETL)                         \
  DATA(vector<float>, electrons_x)                                      \
  DATA(vector<float>, electrons_y)                                      \
  DATA(vector<float>, electrons_z)                                      \
  DATA(vector<float>, electrons_t)                                      \
  DATA(vector<float>, electrons_energy)                                 \
  DATA(vector<int>,   electrons_hasMTD)                                 \
  DATA(vector<float>, electrons_mcMatch_genPdgId)                       \
  DATA(vector<float>, electrons_mcMatch_genPt)                          \
  DATA(vector<float>, electrons_mcMatch_genEta)                         \
  DATA(vector<float>, electrons_mcMatch_genPhi)                         \
  DATA(vector<float>, electrons_mcMatch_genVtx_x)                       \
  DATA(vector<float>, electrons_mcMatch_genVtx_y)                       \
  DATA(vector<float>, electrons_mcMatch_genVtx_z)                       \
  DATA(vector<float>, electrons_mcMatch_genVtx_t)                       \
  DATA(vector<float>, electrons_mcMatch_DR)                             \
  DATA(vector<int>,   matchedSimHits_n)                                 \
  DATA(vector<int>,   matchedRecHits_n)                                 \
  DATA(vector<int>,   matchedClusters_n)                                \
  DATA(vector, matchedSimHits_idx,           <vector<int> >)            \
  DATA(vector, matchedSimHits_det,           <vector<int> >)            \
  DATA(vector, matchedSimHits_energy,        <vector<float> >)          \
  DATA(vector, matchedSimHits_energyCorr,    <vector<float> >)          \
  DATA(vector, matchedSimHits_time,          <vector<float> >)          \
  DATA(vector, matchedSimHits_rr,            <vector<int> >)            \
  DATA(vector, matchedSimHits_module,        <vector<int> >)            \
  DATA(vector, matchedSimHits_modType,       <vector<int> >)            \
  DATA(vector, matchedSimHits_crystal,       <vector<int> >)            \
  DATA(vector, matchedSimHits_ieta,          <vector<int> >)            \
  DATA(vector, matchedSimHits_iphi,          <vector<int> >)            \
  DATA(vector, matchedSimHits_entry_local_x, <vector<float> >)          \
  DATA(vector, matchedSimHits_entry_local_y, <vector<float> >)          \
  DATA(vector, matchedSimHits_entry_local_z, <vector<float> >)          \
  DATA(vector, matchedSimHits_entry_global_R,<vector<float> >)          \
  DATA(vector, matchedSimHits_exit_local_x,  <vector<float> >)          \
  DATA(vector, matchedSimHits_exit_local_y,  <vector<float> >)          \
  DATA(vector, matchedSimHits_exit_local_z,  <vector<float> >)          \
  DATA(vector, matchedSimHits_exit_global_R, <vector<float> >)          \
  DATA(vector, matchedSimHits_electron_Deta,    <vector<float> >)       \
  DATA(vector, matchedSimHits_electron_Dphi,    <vector<float> >)       \
  DATA(vector, matchedSimHits_electron_DR,      <vector<float> >)       \
  DATA(vector, matchedSimHits_electron_Dz,      <vector<float> >)       \
  DATA(vector, matchedSimHits_electron_RDphi,   <vector<float> >)       \
  DATA(vector, matchedSimHits_electron_dist,    <vector<float> >)       \
  DATA(vector, matchedRecHits_idx,           <vector<int> >)            \
  DATA(vector, matchedRecHits_det,           <vector<int> >)            \
  DATA(vector, matchedRecHits_energy,        <vector<float> >)          \
  DATA(vector, matchedRecHits_energyCorr,    <vector<float> >)          \
  DATA(vector, matchedRecHits_time,          <vector<float> >)          \
  DATA(vector, matchedRecHits_rr,            <vector<int> >)            \
  DATA(vector, matchedRecHits_module,        <vector<int> >)            \
  DATA(vector, matchedRecHits_modType,       <vector<int> >)            \
  DATA(vector, matchedRecHits_crystal,       <vector<int> >)            \
  DATA(vector, matchedRecHits_ieta,          <vector<int> >)            \
  DATA(vector, matchedRecHits_iphi,          <vector<int> >)            \
  DATA(vector, matchedRecHits_local_x,       <vector<float> >)          \
  DATA(vector, matchedRecHits_local_y,       <vector<float> >)          \
  DATA(vector, matchedRecHits_local_z,       <vector<float> >)          \
  DATA(vector, matchedRecHits_global_R,      <vector<float> >)          \
  DATA(vector, matchedRecHits_electron_Deta,    <vector<float> >)       \
  DATA(vector, matchedRecHits_electron_Dphi,    <vector<float> >)       \
  DATA(vector, matchedRecHits_electron_DR,      <vector<float> >)       \
  DATA(vector, matchedRecHits_electron_Dz,      <vector<float> >)       \
  DATA(vector, matchedRecHits_electron_RDphi,   <vector<float> >)       \
  DATA(vector, matchedRecHits_electron_dist,    <vector<float> >)       \
  DATA(vector, matchedRecHits_sietaieta,     <vector<float> >)          \
  DATA(vector, matchedRecHits_siphiiphi,     <vector<float> >)          \
  DATA(vector, matchedClusters_idx,           <vector<int> >)           \
  DATA(vector, matchedClusters_det,           <vector<int> >)           \
  DATA(vector, matchedClusters_energy,        <vector<float> >)         \
  DATA(vector, matchedClusters_energyCorr,    <vector<float> >)         \
  DATA(vector, matchedClusters_time,          <vector<float> >)         \
  DATA(vector, matchedClusters_rr,            <vector<int> >)           \
  DATA(vector, matchedClusters_module,        <vector<int> >)           \
  DATA(vector, matchedClusters_modType,       <vector<int> >)           \
  DATA(vector, matchedClusters_crystal,       <vector<int> >)           \
  DATA(vector, matchedClusters_ieta,          <vector<int> >)           \
  DATA(vector, matchedClusters_iphi,          <vector<int> >)           \
  DATA(vector, matchedClusters_size,          <vector<int> >)           \
  DATA(vector, matchedClusters_size_x,        <vector<int> >)           \
  DATA(vector, matchedClusters_size_y,        <vector<int> >)           \
  DATA(vector, matchedClusters_local_x,       <vector<float> >)         \
  DATA(vector, matchedClusters_local_y,       <vector<float> >)         \
  DATA(vector, matchedClusters_local_z,       <vector<float> >)         \
  DATA(vector, matchedClusters_global_R,      <vector<float> >)         \
  DATA(vector, matchedClusters_electron_Deta,    <vector<float> >)      \
  DATA(vector, matchedClusters_electron_Dphi,    <vector<float> >)      \
  DATA(vector, matchedClusters_electron_DR,      <vector<float> >)      \
  DATA(vector, matchedClusters_electron_Dz,      <vector<float> >)      \
  DATA(vector, matchedClusters_electron_RDphi,   <vector<float> >)      \
  DATA(vector, matchedClusters_electron_dist,    <vector<float> >)          

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
