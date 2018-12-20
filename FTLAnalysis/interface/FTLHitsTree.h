#ifndef FTL_HITS_TREE
#define FTL_HITS_TREE

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeBase.h"

using namespace std;

//---Define the TTree branches
#define DYNAMIC_TREE_NAME FTLHitsTree

#define DATA_TABLE                       \
  DATA(int, event)                       \
  DATA(int, lumi)                        \
  DATA(int, run)                         \
  DATA(int, simHits_n)                   \
  DATA(int, recHits_n)                   \
  
#define DATA_CLASS_TABLE                                                \
  DATA(vector<int>,   track_idx)                                        \
  DATA(vector<float>, track_pt)                                         \
  DATA(vector<float>, track_eta)                                        \
  DATA(vector<float>, track_phi)                                        \
  DATA(vector<float>, track_eta_atBTL)                                  \
  DATA(vector<float>, track_phi_atBTL)                                  \
  DATA(vector<float>, track_energy)                                     \
  DATA(vector<float>, track_normalizedChi2)                             \
  DATA(vector<int>,   track_numberOfValidHits)                          \
  DATA(vector<int>,   track_numberOfLostHits)                           \
  DATA(vector<int>,   track_isHighPurity)                               \
  DATA(vector<float>, track_mcMatch_genPdgId)                           \
  DATA(vector<float>, track_mcMatch_genPt)                              \
  DATA(vector<float>, track_mcMatch_genEta)                             \
  DATA(vector<float>, track_mcMatch_genPhi)                             \
  DATA(vector<float>, track_mcMatch_DR)                                 \
  DATA(vector<float>, simHits_energy)                                   \
  DATA(vector<float>, simHits_time)                                     \
  DATA(vector<int>,   simHits_rr)                                       \
  DATA(vector<int>,   simHits_module)                                   \
  DATA(vector<int>,   simHits_modType)                                  \
  DATA(vector<int>,   simHits_crystal)                                  \
  DATA(vector<int>,   simHits_ieta)                                     \
  DATA(vector<int>,   simHits_iphi)                                     \
  DATA(vector<float>, simHits_entry_local_x)                            \
  DATA(vector<float>, simHits_entry_local_y)                            \
  DATA(vector<float>, simHits_entry_local_z)                            \
  DATA(vector<float>, simHits_entry_global_R)                           \
  DATA(vector<float>, simHits_exit_local_x)                             \
  DATA(vector<float>, simHits_exit_local_y)                             \
  DATA(vector<float>, simHits_exit_local_z)                             \
  DATA(vector<float>, simHits_exit_global_R)                            \
  DATA(vector<float>, recHits_energy)                                   \
  DATA(vector<float>, recHits_time)                                     \
  DATA(vector<int>,   recHits_rr)                                       \
  DATA(vector<int>,   recHits_module)                                   \
  DATA(vector<int>,   recHits_modType)                                  \
  DATA(vector<int>,   recHits_crystal)                                  \
  DATA(vector<int>,   recHits_ieta)                                     \
  DATA(vector<int>,   recHits_iphi)                                     \
  DATA(vector<float>, recHits_local_x)                                  \
  DATA(vector<float>, recHits_local_y)                                  \
  DATA(vector<float>, recHits_local_z)                                  \
  DATA(vector<float>, recHits_global_R)                                 \
  DATA(vector<int>,   matchedSimHits_n)                                 \
  DATA(vector<int>,   matchedRecHits_n)                                 \
  DATA(vector, matchedSimHits_idx,           <vector<int> >)            \
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
  DATA(vector, matchedSimHits_track_Deta,    <vector<float> >)          \
  DATA(vector, matchedSimHits_track_Dphi,    <vector<float> >)          \
  DATA(vector, matchedSimHits_track_DR,      <vector<float> >)          \
  DATA(vector, matchedSimHits_track_Dz,      <vector<float> >)          \
  DATA(vector, matchedSimHits_track_RDphi,   <vector<float> >)          \
  DATA(vector, matchedSimHits_track_dist,    <vector<float> >)          \
  DATA(vector, matchedRecHits_idx,           <vector<int> >)            \
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
  DATA(vector, matchedRecHits_track_Deta,    <vector<float> >)          \
  DATA(vector, matchedRecHits_track_Dphi,    <vector<float> >)          \
  DATA(vector, matchedRecHits_track_DR,      <vector<float> >)          \
  DATA(vector, matchedRecHits_track_Dz,      <vector<float> >)          \
  DATA(vector, matchedRecHits_track_RDphi,   <vector<float> >)          \
  DATA(vector, matchedRecHits_track_dist,    <vector<float> >)          \
  DATA(vector, matchedRecHits_sietaieta,     <vector<float> >)          \
  DATA(vector, matchedRecHits_siphiiphi,     <vector<float> >)

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
