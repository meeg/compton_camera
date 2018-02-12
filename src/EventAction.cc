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
// $Id: EventAction.cc 75214 2013-10-29 16:04:42Z gcosmo $
//
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "TrackerHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "Analysis.hh"

#include <unordered_map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* detector)
: G4UserEventAction()
{fDetector = detector;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  // get number of stored trajectories

  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  // periodic printing

  /*
  G4int eventID = event->GetEventID();
  if ( eventID < 1000 || eventID % 1000 == 0) {
    G4cout << ">>> Event: " << eventID  << G4endl;
    if ( trajectoryContainer ) {
      G4cout << "    " << n_trajectories
             << " trajectories stored in this event." << G4endl;
    }
    G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);
    G4cout << "    "  
           << hc->GetSize() << " hits stored in this event" << G4endl;
  }
  */

  G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);

  //sum hits
  std::unordered_map<G4int,std::vector<TrackerHit*>*> **hit_map_array = new std::unordered_map<G4int,std::vector<TrackerHit*>*>*[fDetector->GetNbOfChambers()];
  for (int i=0; i<fDetector->GetNbOfChambers(); i++) {
      hit_map_array[i] = new std::unordered_map<G4int,std::vector<TrackerHit*>*>;
  }
  for (int i=0;i<hc->GetSize();i++) {
      TrackerHit *hit = (TrackerHit*)hc->GetHit(i);
      std::unordered_map<G4int,std::vector<TrackerHit*>*> *hit_map = hit_map_array[hit->GetChamberNb()];
      std::vector<TrackerHit*> *hit_vector;
      if (hit_map->find(hit->GetTrackID()) != hit_map->end()) {
          hit_vector = hit_map->find(hit->GetTrackID())->second;
      } else {
          hit_vector = new std::vector<TrackerHit*>();
          hit_map->emplace(hit->GetTrackID(),hit_vector);
      }
      hit_vector->push_back(hit);
  }

  fChamberNbVec.clear();
  fEdepVec.clear();
  fPosXVec.clear();
  fPosYVec.clear();
  fPosZVec.clear();

  int nhits = 0;

  for (int i=0; i<fDetector->GetNbOfChambers(); i++) {
      for (std::unordered_map<G4int,std::vector<TrackerHit*>*>::iterator it(hit_map_array[i]->begin());it != hit_map_array[i]->end(); it++) {
          std::vector<TrackerHit*>* hit_vector = it->second;
          G4int trackID = hit_vector->front()->GetTrackID();
          G4int chamberNb = hit_vector->front()->GetChamberNb();
          G4double edep = 0.0;
          G4ThreeVector weightedPos = G4ThreeVector(0.0,0.0,0.0);
          for (std::vector<TrackerHit*>::iterator hit_it(hit_vector->begin()); hit_it!=hit_vector->end(); hit_it++) {
              edep += (*hit_it)->GetEdep();
              weightedPos += ((*hit_it)->GetEdep()) * ((*hit_it)->GetPos());
          }
          TrackerHit *combinedHit = new TrackerHit();
          combinedHit->SetTrackID(trackID);
          combinedHit->SetChamberNb(chamberNb);
          combinedHit->SetEdep(edep);
          combinedHit->SetPos(weightedPos/edep);

          nhits++;
          fChamberNbVec.push_back(chamberNb);
          fEdepVec.push_back(edep);
          fPosXVec.push_back(weightedPos.x()/edep);
          fPosYVec.push_back(weightedPos.y()/edep);
          fPosZVec.push_back(weightedPos.z()/edep);

          combinedHit->Print();
      }
  }
  if (nhits>0) {
  auto analysisManager = G4AnalysisManager::Instance();
  //printf("%d\n",analysisManager->GetNofNtuples());
  //printf("%d\n",analysisManager->GetFirstNtupleColumnId());
  analysisManager->FillNtupleIColumn(0, nhits);
  analysisManager->AddNtupleRow();
  }

  for (int i=0; i<fDetector->GetNbOfChambers(); i++) {
      for (std::unordered_map<G4int,std::vector<TrackerHit*>*>::iterator it(hit_map_array[i]->begin());it != hit_map_array[i]->end(); it++) delete (it->second);
      delete hit_map_array[i];
  }
  delete hit_map_array;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
