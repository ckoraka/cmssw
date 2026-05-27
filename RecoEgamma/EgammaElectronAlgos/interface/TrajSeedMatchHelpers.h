#ifndef RecoEgamma_EgammaElectronAlgos_TrajSeedMatchHelpers_h
#define RecoEgamma_EgammaElectronAlgos_TrajSeedMatchHelpers_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/TrajSeedMatchingCutsPortable.h"

#include <vector>
#include <array>

namespace egamma {

  inline edm::ParameterSetDescription makeMatchingCutsPSetDescription() {
    edm::ParameterSetDescription desc;

    edm::ParameterSetDescription cutsDesc;
    cutsDesc.add<std::vector<double> >("dPhiMaxHighEt",       {0.003});
    cutsDesc.add<std::vector<double> >("dPhiMaxHighEtThres",  {0.0});
    cutsDesc.add<std::vector<double> >("dPhiMaxLowEtGrad",    {0.0});
    cutsDesc.add<std::vector<double> >("dRZMaxHighEt",        {0.005});
    cutsDesc.add<std::vector<double> >("dRZMaxHighEtThres",   {30});
    cutsDesc.add<std::vector<double> >("dRZMaxLowEtGrad",     {-0.002});
    cutsDesc.add<std::vector<double> >("etaBins",             {});

    edm::ParameterSet defaults;
    defaults.addParameter<std::vector<double> >("dPhiMaxHighEt",      std::vector<double>{0.003});
    defaults.addParameter<std::vector<double> >("dPhiMaxHighEtThres", std::vector<double>{0.0});
    defaults.addParameter<std::vector<double> >("dPhiMaxLowEtGrad",   std::vector<double>{0.0});
    defaults.addParameter<std::vector<double> >("dRZMaxHighEt",       std::vector<double>{0.005});
    defaults.addParameter<std::vector<double> >("dRZMaxHighEtThres",  std::vector<double>{30});
    defaults.addParameter<std::vector<double> >("dRZMaxLowEtGrad",    std::vector<double>{-0.002});
    defaults.addParameter<std::vector<double> >("etaBins",            std::vector<double>{});

    desc.addVPSet("matchingCuts", cutsDesc, std::vector<edm::ParameterSet>{defaults, defaults, defaults});
    return desc;
  }

}


#endif
