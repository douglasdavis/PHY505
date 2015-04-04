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
// $Id: RunAction.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
//#include "Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
  : G4UserRunAction()
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     
  G4RunManager::GetRunManager()->SetPrintProgress(1000);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  fOutFile = new TFile("zout.root","RECREATE");
  fOutTree = new TTree("RunTree","RunTree");

  fOutTree->Branch("L1_Edep",&fl1e,"L1_Edep/D");
  fOutTree->Branch("L1_Tlen",&fl1t,"L1_Tlen/D");
  fOutTree->Branch("L2_Edep",&fl2e,"L2_Edep/D");
  fOutTree->Branch("L2_Tlen",&fl2t,"L2_Tlen/D");
  fOutTree->Branch("L3_Edep",&fl3e,"L3_Edep/D");
  fOutTree->Branch("L3_Tlen",&fl3t,"L3_Tlen/D");
  fOutTree->Branch("L4_Edep",&fl4e,"L4_Edep/D");
  fOutTree->Branch("L4_Tlen",&fl4t,"L4_Tlen/D");  
  fOutTree->Branch("L5_Edep",&fl5e,"L5_Edep/D");
  fOutTree->Branch("L5_Tlen",&fl5t,"L5_Tlen/D");
  fOutTree->Branch("L6_Edep",&fl6e,"L6_Edep/D");
  fOutTree->Branch("L6_Tlen",&fl6t,"L6_Tlen/D");
  fOutTree->Branch("L7_Edep",&fl7e,"L7_Edep/D");
  fOutTree->Branch("L7_Tlen",&fl7t,"L7_Tlen/D");
  fOutTree->Branch("L8_Edep",&fl8e,"L8_Edep/D");
  fOutTree->Branch("L8_Tlen",&fl8t,"L8_Tlen/D");  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::UpdateTree(const G4double l1e, const G4double l1t,
			   const G4double l2e, const G4double l2t,
			   const G4double l3e, const G4double l3t,
			   const G4double l4e, const G4double l4t,
			   const G4double l5e, const G4double l5t,
			   const G4double l6e, const G4double l6t,
			   const G4double l7e, const G4double l7t,
			   const G4double l8e, const G4double l8t)
			   
{
  fl1e = l1e;
  fl1t = l1t;

  fl2e = l2e;
  fl2t = l2t;

  fl3e = l3e;
  fl3t = l3t;

  fl4e = l4e;
  fl4t = l4t;

  fl5e = l5e;
  fl5t = l5t;

  fl6e = l6e;
  fl6t = l6t;

  fl7e = l7e;
  fl7t = l7t;

  fl8e = l8e;
  fl8t = l8t;

  fOutTree->Fill();
}
    
void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  fOutTree->Write();
  fOutFile->Close(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
