#include "detector.hh"

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{}

MySensitiveDetector::~MySensitiveDetector()
{}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{
	G4Track *track =aStep->GetTrack();
	
	track->SetTrackStatus(fStopAndKill);
	
	G4StepPoint *preStepPoint = aStep->GetPreStepPoint();  //ENTER DETECTOR
	G4StepPoint *postStepPoint = aStep->GetPostStepPoint(); //LEAVE DETECTOR
	
	G4ThreeVector posPhoton = preStepPoint->GetPosition();	
		
	//G4cout <<posPhoton<<G4endl;
	
	G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
	
	
	G4AnalysisManager *man = G4AnalysisManager::Instance();
	man->FillNtupleIColumn(0,evt);
	man->FillNtupleDColumn(1,posPhoton[0]);
	man->FillNtupleDColumn(2,posPhoton[1]);
	man->FillNtupleDColumn(3,posPhoton[2]);
	man->AddNtupleRow(0);
	
}


