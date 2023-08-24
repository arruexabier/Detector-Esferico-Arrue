#include "generator.hh"

MyPrimaryGenerator::MyPrimaryGenerator()
{
	
	fParticleGun = new G4ParticleGun(1);
	
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName="opticalphoton";
	G4ParticleDefinition *particle = particleTable->FindParticle("opticalphoton");
	
	fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-0.015*m));
	fParticleGun->SetParticleEnergy(7.293187824072908*eV);
	fParticleGun->SetParticleDefinition(particle);	
	fParticleGun->SetParticlePolarization(G4RandomDirection());	

}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
	delete fParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{

    fParticleGun->SetParticleMomentumDirection(G4RandomDirection());

	fParticleGun->GeneratePrimaryVertex(anEvent);
}
