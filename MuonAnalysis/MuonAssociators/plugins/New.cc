// -*- C++ -*-
// //
// // Package: MuonMCClassifAndAna
// // Class: MuonMCClassifAndAna
// //
// /**\class MuonMCClassifAndAna MuonMCClassifAndAna.cc MuonAnalysis/MuonAssociators/src/MuonMCClassifAndAna.cc
//
//

/*

CLASSIFICATION: For each RECO Muon, match to SIM particle, and then:
- If the SIM is not a Muon, label as Punchthrough (1)
- If the SIM is a Muon, then look at it's provenance.
A) the SIM muon is also a GEN muon, whose parent is NOT A HADRON AND NOT A TAU
-> classify as "primary" (4).
B) the SIM muon is also a GEN muon, whose parent is HEAVY FLAVOURED HADRON OR A TAU
-> classify as "heavy flavour" (3)
C) classify as "light flavour/decay" (2)
In any case, if the TP is not preferentially matched back to the same RECO muon,
label as Ghost (flip the classification)
FLAVOUR:
- for non-muons: 0
- for primary muons: 13
- for non primary muons: flavour of the mother: abs(pdgId) of heaviest quark, or 15 for tau
*/


//
//// Original Author: Nov 16 16:12 (lxplus231.cern.ch)
//// Created: Sun Nov 16 16:14:09 CET 2008
//// $Id: MuonMCClassifAndAna.cc,v 1.9 2013/06/24 12:53:19 speer Exp $
////
////
//// system include files
#include <memory>
#include <set>


// user include files

 
