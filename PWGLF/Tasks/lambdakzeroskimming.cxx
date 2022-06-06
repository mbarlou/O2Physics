// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

// Example V0 skimming task
// ========================
//
// This code loops over FullV0s
// produces Filtered V0s

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/StrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "MathUtils/Utils.h"
#include "LFDerived.h"
#include <TLorentzVector.h>
#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::math_utils::detail;
using std::array;

#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

namespace o2::aod
{

namespace lfv0
{
DECLARE_SOA_INDEX_COLUMN(LFCollision, lfCollision);
// DECLARE_SOA_COLUMN(LFCollision, lfCollision);
DECLARE_SOA_COLUMN(Spec, spec, int8_t);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float pt, float eta) { return pt * TMath::CosH(eta); });
// add the V0 quantities

} // namespace lfv0

/*DECLARE_SOA_TABLE(LFCollisions, "AOD", "LFCOLLISION", o2::soa::Index<>,
                  o2::aod::collision::PosZ);
using LFCollision = LFCollisions::iterator;*/

DECLARE_SOA_TABLE(LFV0s, "AOD", "LFV0s", o2::soa::Index<>,
                  lfv0::Spec,
                  lfv0::Pt, lfv0::Eta, lfv0::Phi,
                  lfv0::P<lfv0::Pt, lfv0::Eta>);
using LFV0 = LFV0s::iterator;

} // namespace o2::aod

struct lambdakzeroSkimmingTask {

  Produces<aod::LFCollisions> outputFilteredCollisions;
  Produces<aod::LFV0s> outputFilteredV0s;

  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", .1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", .1, "DCA Pos To PV"};
  Configurable<float> v0radius{"v0radius", 5.0, "v0radius"};
  Configurable<float> rapidity{"rapidity", 0.5, "rapidity"};
  Configurable<float> v0eta{"v0eta", 0.8, "v0eta"};
  Configurable<int> saveDcaHist{"saveDcaHist", 0, "saveDcaHist"};

  Configurable<float> k0s_dcanegtopv{"k0s_dcanegtopv", .1, "k0s DCA Neg To PV"};
  Configurable<float> k0s_dcapostopv{"k0s_dcapostopv", .1, "K0s DCA Pos To PV"};

  Configurable<float> lam_dcanegtopv{"lam_dcanegtopv", 0.25, "Lambdas DCA Neg To PV"};
  Configurable<float> lam_dcapostopv{"lam_dcapostopv", .1, "Lambdas DCA Pos To PV"};

  HistogramRegistry histos{"Histos", {}, QAObject};

  void init(o2::framework::InitContext&)
  {

    histos.add("events", "", kTH1F, {{10, 0., 10., "Number of events"}});

    histos.add("lam_pt", "", kTH1F, {{40, 0., 10., "lam Pt"}});
    histos.add("lam_eta", "", kTH1F, {{100, -1., 1., "lam Eta"}});
    histos.add("lam_phi", "", kTH1F, {{100, 0., 180., "lam Phi"}});
    histos.add("lam_mass", "", kTH1F, {{3000, 0.0, 3.0, "lam Inv. Mass"}});

    histos.add("K0s_pt", "", kTH1F, {{40, 0., 10., "K0s Pt"}});
    histos.add("K0s_eta", "", kTH1F, {{100, -1., 1., "K0s Eta"}});
    histos.add("K0s_phi", "", kTH1F, {{100, 0., 180., "K0s Phi"}});
    histos.add("K0s_mass", "", kTH1F, {{3000, 0.0, 3.0, "K0s Inv. Mass"}});

    histos.add("V0_pt_unfiltered", "", kTH1F, {{40, 0., 10., "V0 pt unf"}});
    histos.add("V0_eta_unfiltered", "", kTH1F, {{100, -1., 1., "V0 eta unf"}});
    histos.add("V0_phi_unfiltered", "", kTH1F, {{100, 0., 180., "V0 phi unf"}});
    histos.add("K0s_mass_unfiltered", "", kTH1F, {{3000, 0.0, 3.0, "K0s Inv. Mass unf"}});
    histos.add("Lam_mass_unfiltered", "", kTH1F, {{3000, 0.0, 3.0, "Lambdas Inv. Mass unf"}});
  }
  static constexpr float defaultLifetimeCuts[1][2] = {{25., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

  Filter V0Filter = nabs(aod::v0data::dcapostopv) > dcapostopv&& nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters < dcav0dau;
  // more V0 selection creteria will be added here
  using TracksWithPID = soa::Join<aod::FullTracks, aod::pidTOFPi, aod::pidTPCPi>;
  void process(soa::Join<aod::Collisions, aod::EvSels, aod::CentV0Ms>::iterator const& collision,
               //  o2::aod::V0 const&,
               soa::Filtered<aod::V0Datas> const& fullV0s,
               TracksWithPID const&)
  // soa::aod::PosTrack const& posTrack, soa::aod::NegTrack const& negTrack)
  {
    histos.fill(HIST("events"), 1.);

    /* if (!collision.alias()[kINT7]) {
       return;
     }
     if (!collision.sel7()) {
        return;
      }*/

    histos.fill(HIST("events"), 2.);
    outputFilteredCollisions(collision.posZ());
    for (auto& v0 : fullV0s) {

      histos.fill(HIST("V0_pt_unfiltered"), v0.pt());
      histos.fill(HIST("V0_eta_unfiltered"), v0.eta());
      histos.fill(HIST("V0_phi_unfiltered"), v0.phi());
      histos.fill(HIST("K0s_mass_unfiltered"), v0.mK0Short());
      histos.fill(HIST("Lam_mass_unfiltered"), v0.mLambda());

      // const auto& pos = v0.negTrack_as<TracksWithPID>();
      // const auto& neg = v0.posTrack_as<TracksWithPID>();

      if (v0.dcaV0daughters() < dcav0dau) {
        if (TMath::Abs(v0.eta()) < v0eta) {
          if (v0.v0radius() > v0radius && v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa) {
            if (TMath::Abs(v0.yLambda()) < rapidity) {
              if (v0.dcanegtopv() < lam_dcanegtopv && v0.dcapostopv() < lam_dcapostopv) {
                if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(kLambda0) < lifetimecut->get("lifetimecutLambda")) {

                  // add V0 selection creteria
                  outputFilteredV0s(0, v0.pt(), v0.eta(), v0.phi());

                  histos.fill(HIST("lam_pt"), v0.pt());
                  histos.fill(HIST("lam_eta"), v0.eta());
                  histos.fill(HIST("lam_phi"), v0.phi());
                  histos.fill(HIST("lam_mass"), v0.mLambda());
                }
              }
            }

            if (TMath::Abs(v0.yK0Short()) < rapidity) {
              if (v0.dcanegtopv() < k0s_dcanegtopv && v0.dcapostopv() < k0s_dcapostopv) {
                if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(kK0Short) < lifetimecut->get("lifetimecutK0S")) {
                  // if (v0.qtarm(collision.pxpos(), collision.pypos(), collision.pzpos(), collision.pxneg(), collision.pyneg(), collision.pzneg()) > 2. * TMath::Abs(v0.alpha(collision.pxpos(), collision.pypos(), collision.pzpos(), collision.pxneg(), collision.pyneg(), collision.pzneg()))) {
                  //  if (v0.qtarm() > 2.) {
                  // add V0 selection creteria
                  outputFilteredV0s(1, v0.pt(), v0.eta(), v0.phi());
                  histos.fill(HIST("K0s_pt"), v0.pt());
                  histos.fill(HIST("K0s_eta"), v0.eta());
                  histos.fill(HIST("K0s_phi"), v0.phi());
                  histos.fill(HIST("K0s_mass"), v0.mK0Short());
                  // }
                }
              }
            }
          }
        }
      }
    }
  }

  //	TFile* DataRawYield_Var = new TFile(Form("ParameterNo[%d].root", cut),"recreate");
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdakzeroSkimmingTask>(cfgc)};
}
