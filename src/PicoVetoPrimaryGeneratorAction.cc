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
/// \file PicoVeto/src/PicoVetoPrimaryGeneratorAction.cc
/// \brief Implementation of the PicoVetoPrimaryGeneratorAction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PicoVetoPrimaryGeneratorAction.hh"
#include "PicoVetoPrimaryGeneratorMessenger.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>

#include "PicoVetoSensitiveDetector.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PicoVetoPrimaryGeneratorAction::PicoVetoPrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(), 
   fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  //create a messenger for this class
  fGunMessenger = new PicoVetoPrimaryGeneratorMessenger(this);

  
  //default kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
   G4ParticleDefinition* particle = particleTable->FindParticle("e-");
  //  G4ParticleDefinition* particle = particleTable->FindParticle("mu-"); //muon
  // G4ParticleDefinition* particle = particleTable->FindParticle("opticalphoton");

  
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleTime(0.0*ns);

    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PicoVetoPrimaryGeneratorAction::~PicoVetoPrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PicoVetoPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  /*
  // monoenergetic beam in -y direction : 
  //fParticleGun->SetParticlePosition(G4ThreeVector(0.0*cm,1000*cm,0.0*cm));
  fParticleGun->SetParticlePosition(G4ThreeVector(0.0*cm,220*cm,0.0*cm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,-1,0));
  fParticleGun->SetParticleEnergy(1*GeV);
  */
  
  // fParticleGun->SetParticlePosition(G4ThreeVector(100.0*cm,150*cm,0.0*cm));
  // fParticleGun->SetParticlePosition(G4ThreeVector(50.0*cm,280*cm,0.0*cm));
  // fParticleGun->SetParticlePosition(G4ThreeVector(2.0*m,3.7*m,0.0*m));
    fParticleGun->SetParticlePosition(G4ThreeVector(2.*m,2.7*m,0.0*m));

  // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.5,-1,0));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,-1,0));
  //  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-0.5,1,0));
  
  fParticleGun->SetParticleEnergy(0.5*GeV);
 



  // -------------------------------------------------------------------------------
  // SNOLAB muon spectrum (from DOI: 10.1103/PhysRevD.73.053004):
  
   // initial position
  G4double randomP = G4UniformRand(); // entre 0 et 1
  if(randomP==0) randomP = 1E-33;
  G4double randomPhi = G4UniformRand()*360*deg; // entre 0 and 360 deg 
  G4double angle = -sqrt( (0.267447*0.267447*log(4.10762/randomP))/0.5 )+1.44933+0.00024299999999999322;// cosTheta distribution fit 2
  G4double theta = G4UniformRand()*90*deg*angle; // between 0 and 90 deg with prob towards 0 from angle...
  G4double rayon = 4.0;
  G4double posX = rayon * cos(randomPhi) * sin(acos(angle));
  G4double posZ = rayon * sin(randomPhi) * sin(acos(angle));
  G4double posY = rayon * angle - 2.0;

    
  // initial momentum towards center of bottom of tank : 
  G4double momX = -posX;
  G4double momY = -posY-2;
  G4double momZ = -posZ;
   

  // initial energy : 
  G4double Energie = 0.;
  G4double b = 0.4;
  G4double gamma = 3.77;
  G4double epsilon = 693.;
  G4double Anorm = (10000000000000./7.812132)*22.;
  G4double h = 6.011;
  G4double randomE = G4UniformRand()+0.00001; // entre 0.00001 et 1
  Energie = pow((Anorm*exp(-b*h*(gamma-1.))/randomE),(1./gamma)) - epsilon*(1-exp(-b*h));


  
  // log position of all particle source in a .txt file : 
  std::ofstream log("sourcePOS.txt", std::ios_base::app | std::ios_base::out);
  log << posX << "\t" << posY << "\t" << posZ << "\n";
  // G4cout << G4endl; 
  // G4cout <<posX << "\t" << posY << "\t" << posZ  << G4endl; 
  
  // log energy of all particle source in a .txt file : 
  std::ofstream logE("sourceENE.txt", std::ios_base::app | std::ios_base::out);
  logE << Energie << "\n";
  // G4cout << G4endl; 
  // G4cout <<posX << "\t" << posY << "\t" << posZ  << G4endl; 


  

  /*
  // Set Source Data : 
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(momX,momY,momZ));
  fParticleGun->SetParticlePosition(G4ThreeVector(posX*m,posY*m,posZ*m));
  fParticleGun->SetParticleEnergy(Energie*GeV);
  */


  
   //
   fParticleGun->GeneratePrimaryVertex(anEvent);


  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PicoVetoPrimaryGeneratorAction::SetOptPhotonPolar()
{
 G4double angle = G4UniformRand() * 360.0*deg;
 SetOptPhotonPolar(angle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PicoVetoPrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
 if (fParticleGun->GetParticleDefinition()->GetParticleName()!="opticalphoton")
   {
     G4cout << "--> warning from PrimaryGeneratorAction::SetOptPhotonPolar() :"
               "the particleGun is not an opticalphoton" << G4endl;
     return;
   }

 G4ThreeVector normal (1., 0., 0.);
 G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
 G4ThreeVector product = normal.cross(kphoton);
 G4double modul2       = product*product;
 
 G4ThreeVector e_perpend (0., 0., 1.);
 if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product;
 G4ThreeVector e_paralle    = e_perpend.cross(kphoton);
 
 G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
 fParticleGun->SetParticlePolarization(polar);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
