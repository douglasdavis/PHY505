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
// $Id: EventAction.cc 75604 2013-11-04 13:17:26Z gcosmo $
// 
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "RunAction.hh"
#include "CalorimeterSD.hh"
#include "CalorHit.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
  : G4UserEventAction(),
    fLay1HCID(-1),
    fLay2HCID(-1),
    fLay3HCID(-1),
    fLay4HCID(-1),
    fLay5HCID(-1),
    fLay6HCID(-1),
    fLay7HCID(-1),
    fLay8HCID(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorHitsCollection* 
EventAction::GetHitsCollection(G4int hcID,
			       const G4Event* event) const
{
  CalorHitsCollection* hitsCollection = static_cast<CalorHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("EventAction::GetHitsCollection()","MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void EventAction::PrintEventStatistics(G4double lay1Edep, G4double lay1TrackLength,
				       G4double lay2Edep, G4double lay2TrackLength) const
{
  // print event statistics
  G4cout
    << "   Lay1rber: total energy: " 
    << std::setw(7) << G4BestUnit(lay1Edep, "Energy")
    << "       total track length: " 
    << std::setw(7) << G4BestUnit(lay1TrackLength, "Length")
    << G4endl
    << "        Lay2: total energy: " 
    << std::setw(7) << G4BestUnit(lay2Edep, "Energy")
    << "       total track length: " 
    << std::setw(7) << G4BestUnit(lay2TrackLength, "Length")
    << G4endl;
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{  
  // Get hits collections IDs (only once)
  if ( fLay1HCID == -1 ) {
    fLay1HCID = G4SDManager::GetSDMpointer()->GetCollectionID("Lay1HitsCollection");
    fLay2HCID = G4SDManager::GetSDMpointer()->GetCollectionID("Lay2HitsCollection");
    fLay3HCID = G4SDManager::GetSDMpointer()->GetCollectionID("Lay3HitsCollection");
    fLay4HCID = G4SDManager::GetSDMpointer()->GetCollectionID("Lay4HitsCollection");
    fLay5HCID = G4SDManager::GetSDMpointer()->GetCollectionID("Lay5HitsCollection");
    fLay6HCID = G4SDManager::GetSDMpointer()->GetCollectionID("Lay6HitsCollection");
    fLay7HCID = G4SDManager::GetSDMpointer()->GetCollectionID("Lay7HitsCollection");
    fLay8HCID = G4SDManager::GetSDMpointer()->GetCollectionID("Lay8HitsCollection");
  }

  // Get hits collections
  CalorHitsCollection* lay1HC = GetHitsCollection(fLay1HCID, event);
  CalorHitsCollection* lay2HC = GetHitsCollection(fLay2HCID, event);
  CalorHitsCollection* lay3HC = GetHitsCollection(fLay3HCID, event);
  CalorHitsCollection* lay4HC = GetHitsCollection(fLay4HCID, event);
  CalorHitsCollection* lay5HC = GetHitsCollection(fLay5HCID, event);
  CalorHitsCollection* lay6HC = GetHitsCollection(fLay6HCID, event);
  CalorHitsCollection* lay7HC = GetHitsCollection(fLay7HCID, event);
  CalorHitsCollection* lay8HC = GetHitsCollection(fLay8HCID, event);

  // Get hit with total values
  CalorHit* lay1Hit = (*lay1HC)[lay1HC->entries()-1];
  CalorHit* lay2Hit = (*lay2HC)[lay2HC->entries()-1];
  CalorHit* lay3Hit = (*lay3HC)[lay3HC->entries()-1];
  CalorHit* lay4Hit = (*lay4HC)[lay4HC->entries()-1];
  CalorHit* lay5Hit = (*lay5HC)[lay5HC->entries()-1];
  CalorHit* lay6Hit = (*lay6HC)[lay6HC->entries()-1];
  CalorHit* lay7Hit = (*lay7HC)[lay7HC->entries()-1];
  CalorHit* lay8Hit = (*lay8HC)[lay8HC->entries()-1];
 
  // Print per event (modulo n)
  //
  G4int eventID = event->GetEventID();
  G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;     
  }  

  RunAction* runac = (RunAction*)(G4RunManager::GetRunManager()->GetUserRunAction());
  runac->UpdateTree(lay1Hit->GetEdep(),lay1Hit->GetTrackLength(),
		    lay2Hit->GetEdep(),lay2Hit->GetTrackLength(),
		    lay3Hit->GetEdep(),lay3Hit->GetTrackLength(),
		    lay4Hit->GetEdep(),lay4Hit->GetTrackLength(),
		    lay5Hit->GetEdep(),lay5Hit->GetTrackLength(),
		    lay6Hit->GetEdep(),lay6Hit->GetTrackLength(),
		    lay7Hit->GetEdep(),lay7Hit->GetTrackLength(),
		    lay8Hit->GetEdep(),lay8Hit->GetTrackLength());
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
