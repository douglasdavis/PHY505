//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: RunAction.hh 74265 2013-10-02 14:41:20Z gcosmo $
// 
/// \file RunAction.hh
/// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include <vector>

class G4Run;

/// Run action class
///
/// It accumulates statistic and computes dispersion of the energy deposit 
/// and track lengths of charged particles with use of analysis tools:
/// H1D histograms are created in BeginOfRunAction() for the following 
/// physics quantities:
/// - Edep in absorber
/// - Edep in gap
/// - Track length in absorber
/// - Track length in gap
/// The same values are also saved in the ntuple.
/// The histograms and ntuple are saved in the output file in a format
/// accoring to a selected technology in Analysis.hh.
///
/// In EndOfRunAction(), the accumulated statistic and computed 
/// dispersion is printed.
///

class RunAction : public G4UserRunAction {

public:

  RunAction();
  virtual ~RunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void   EndOfRunAction(const G4Run*);

  void UpdateTree(const G4double l1e, const G4double l1t,
		  const G4double l2e, const G4double l2t,
		  const G4double l3e, const G4double l3t,
		  const G4double l4e, const G4double l4t,
		  const G4double l5e, const G4double l5t,
		  const G4double l6e, const G4double l6t,
		  const G4double l7e, const G4double l7t,
		  const G4double l8e, const G4double l8t);

  
private:

  TFile *fOutFile;
  TTree *fOutTree;

  G4double fl1e;
  G4double fl1t;
  G4double fl2e;
  G4double fl2t;
  G4double fl3e;
  G4double fl3t;
  G4double fl4e;
  G4double fl4t;
  G4double fl5e;
  G4double fl5t;
  G4double fl6e;
  G4double fl6t;
  G4double fl7e;
  G4double fl7t;
  G4double fl8e;
  G4double fl8t;

};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

