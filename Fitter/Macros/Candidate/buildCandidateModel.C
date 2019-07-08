#ifndef Candidate_buildCandidateModel_C
#define Candidate_buildCandidateModel_C


#include "RooWorkspace.h"
#include "RooArgSet.h"

#include <iostream>
#include <string>

#include "addModel.C"
#include "../Utilities/rooDataUtils.h"


bool buildCandidateModel(RooWorkspace& ws, GlobalInfo&  info, const std::string& chg)
{
  // Import the Candidate Models to the local workspace
  if (!addModel(ws, info, chg)) { return false; }
  //
  // Set Fixed parameters to constant (clean up)
  setFixedVarsToContantVars(ws);
  //
  // save the initial values of the model we've just created
  saveSnapshot(ws, "initialParameters", info.Par.at("dsName"+chg));
  //
  return true;
};


#endif // #ifndef Candidate_buildCandidateMassModel_C
