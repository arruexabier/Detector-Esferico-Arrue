#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4OpticalPhoton.hh"
#include "G4RandomDirection.hh"
#include "Randomize.hh"
#include "G4GenericMessenger.hh"
#include "G4RandomDirection.hh"

class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
	MyPrimaryGenerator();
	~MyPrimaryGenerator();
	
	virtual void GeneratePrimaries(G4Event*);
	
private:
	G4ParticleGun *fParticleGun;
};

#endif
