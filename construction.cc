#include "construction.hh"

MyDetectorConstruction::MyDetectorConstruction()
{
	fMessenger = new G4GenericMessenger(this,"/config/","Configuration of Detector");
	
	fMessenger->DeclareProperty("DetPos",DetPos,"Position of Detector");
	fMessenger->DeclareProperty("CaseA",CaseA,"CaseA");
	fMessenger->DeclareProperty("CaseB",CaseB,"CaseB");
	fMessenger->DeclareProperty("CaseC",CaseC,"CaseC");
	fMessenger->DeclareProperty("CaseD",CaseD,"CaseD");
	fMessenger->DeclareProperty("PMT",PMTR,"Detector PmtR11410");
	fMessenger->DeclareProperty("zGun",zGun,"Posicion z del Emisor");
	fMessenger->DeclareProperty("Gas",Gas,"Añadir Xenon gas en pared");
	fMessenger->DeclareProperty("transparency",transparency,"Transparencia de Mesh");
	fMessenger->DeclareProperty("Patron1",Patron1,"Patron1");
	fMessenger->DeclareProperty("Patron2",Patron2,"Patron2");
	fMessenger->DeclareProperty("Patron3",Patron3,"Patron3");
	
	DetPos=-0*cm;
	CaseA = true;  //TEFLON EN ESFERA
	CaseB = true;  //TEFLON EN TAPA
	CaseC = true;  //TEFLON EN TUBO INTERIOR
	CaseD = true;  //TEFLON EN TUBO EXTERIOR

	PMTR=false;

	transparency=1.0;
	zGun=-0.015*m;
	
	Gas=false;
	
	Patron1=false;
	Patron2=false;
	Patron3=false;
	
	DefineMaterials();

		
}

MyDetectorConstruction::~MyDetectorConstruction()
{}




G4double XenonRefractiveIndex(G4double energy, G4double density)
{
  G4double P[3] = {71.23, 77.75, 1384.89}; // [eV^3 cm3 / mole]
  G4double E[3] = {8.4, 8.81, 13.2};       // [eV]

  energy = energy / eV;
  G4double virial = 0.;

  for (G4int i=0; i<3; i++)
  virial = virial + P[i] / (energy*energy - E[i]*E[i]);

  density = density / g * cm3;

  G4double mol_density = density / 131.29;
  G4double alpha = virial * mol_density;

  G4double n2 = (1. - 2*alpha) / (1. + alpha);
  if (n2 < 1.) {
    n2 = 1.;
  }

  return sqrt(n2);
}

G4double GXeDensity(G4double pressure)
{
  // Computes Xe (gas) density at T = 293 K
  // Values are taken from the reference file nexus/data/gxe_density_table.txt
  // (which, in turn, is downloaded from https://webbook.nist.gov/chemistry/fluid).
  // We assume a linear interpolation between any pair of values in the database.

  G4double density;
  const G4int n_pressures = 6;
  G4double data[n_pressures][2] = {{  1.0 * bar,   5.419 * kg/m3},
                                   {  5.0 * bar,  27.721 * kg/m3},
                                   { 10.0 * bar,  57.160 * kg/m3},
                                   { 13.5 * bar,  78.949 * kg/m3},
                                   { 20.0 * bar, 122.510 * kg/m3},
                                   { 30.0 * bar, 199.920 * kg/m3}};
  G4bool found = false;

  for (G4int i=0; i<n_pressures-1; ++i) {
    if  (pressure >= data[i][0] && pressure < data[i+1][0]) {
      G4double x1 = data[i][0];
      G4double x2 = data[i+1][0];
      G4double y1 = data[i][1];
      G4double y2 = data[i+1][1];
      density = y1 + (y2-y1)*(pressure-x1)/(x2-x1);
      found = true;
      break;
    }
  }

  if (!found) {
    if (pressure == data[n_pressures-1][0]) {
      density = data[n_pressures-1][1];
    }
    else {
      throw "Unknown xenon density for this pressure!";
    }
  }

  return density;
}




void MyDetectorConstruction::DefineMaterials()
{
	G4NistManager *nist = G4NistManager::Instance();
	
	
	//XENON
	lXe = nist->FindOrBuildMaterial("G4_lXe");
	
	G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();
	
	constexpr G4double optPhotMinE_ = 6 * eV;
	constexpr G4double optPhotMaxE_ = 9  * eV;
	constexpr G4double noAbsLength_ = 1.e8  * m;
	
	//REFRACTIVE INDEX
	const G4int ri_entries = 200;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;
	
	std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }
	
	G4double density = 2.953 * g/cm3;
    std::vector<G4double> ri_index;
    for (G4int i=0; i<ri_entries; i++) {
      ri_index.push_back(XenonRefractiveIndex(ri_energy[i], density));
    }
	mptWorld->AddProperty("RINDEX",ri_energy,ri_index);
	
	//ABSORTION LENGTH
	std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> abs_length = {noAbsLength_, noAbsLength_};
    mptWorld->AddProperty("ABSLENGTH", abs_energy, abs_length);
		
	//RAYLEIGH
	std::vector<G4double> rayleigh_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> rayleigh_length = {36.*cm, 36.*cm};
    mptWorld->AddProperty("RAYLEIGH", rayleigh_energy, rayleigh_length);
    
	lXe->SetMaterialPropertiesTable(mptWorld);
	
	
	
	

	
	//STEEL
	Steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");

	G4MaterialPropertiesTable *mptSteel = new G4MaterialPropertiesTable();
	
	G4double energy[2]={optPhotMinE_,optPhotMaxE_};
	G4double rindexSteel[2]={1.29,1.29};
	mptSteel->AddProperty("RINDEX",energy,rindexSteel,2);
	
	std::vector<G4double> ENERGIESSteel = {optPhotMinE_,optPhotMaxE_};
	std::vector<G4double> REFLECTIVITYSteel = {0.98, 0.98};
	mptSteel->AddProperty("REFLECTIVITY", ENERGIESSteel, REFLECTIVITYSteel);
	Steel->SetMaterialPropertiesTable(mptSteel);
	
	

	//TEFLON
	Teflon = nist->FindOrBuildMaterial("G4_TEFLON");

	G4MaterialPropertiesTable *mptTeflon = new G4MaterialPropertiesTable();
	// REFLECTIVITY IN LXE (from https://link.springer.com/content/pdf/10.1140/epjc/s10052-020-7800-6.pdf)
	std::vector<G4double> ENERGIES = {optPhotMinE_,optPhotMaxE_};
	std::vector<G4double> REFLECTIVITY = {0.98, 0.98};
	mptTeflon->AddProperty("REFLECTIVITY", ENERGIES, REFLECTIVITY);

	// REFLEXION BEHAVIOR
	// Specular reflection about the normal to a microfacet.
	// Such a vector is chosen according to a gaussian distribution with
	// sigma = SigmaAlhpa (in rad) and centered in the average normal.
	std::vector<G4double> specularlobe  = {0., 0.};
	// specular reflection about the average normal
	std::vector<G4double> specularspike = {0., 0.};
	// 180 degrees reflection.
	std::vector<G4double> backscatter   = {0., 0.};
	// 1 - the sum of these three last parameters is the percentage of Lambertian reflection

	mptTeflon->AddProperty("SPECULARLOBECONSTANT", ENERGIES, specularlobe);
	mptTeflon->AddProperty("SPECULARSPIKECONSTANT",ENERGIES, specularspike);
	mptTeflon->AddProperty("BACKSCATTERCONSTANT",  ENERGIES, backscatter);

	Teflon->SetMaterialPropertiesTable(mptTeflon);
	
	
	
	
	
	
	//KOVAR
	Kovar = nist->FindOrBuildMaterial("G4_Kovar");

	G4Element* Fe = nist->FindOrBuildElement("Fe");
	G4Element* Ni = nist->FindOrBuildElement("Ni");
	G4Element* Co = nist->FindOrBuildElement("Co");

	Kovar = new G4Material("Kovar", 8.35*g/cm3, 3);
	Kovar->AddElement(Fe, 54);
	Kovar->AddElement(Ni, 29);
	Kovar->AddElement(Co, 17);
	
	Kovar->SetMaterialPropertiesTable(new G4MaterialPropertiesTable());
	
	
	
	
	//FUSED SILICA
	silica = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	
	G4MaterialPropertiesTable* mptsilica = new G4MaterialPropertiesTable();
    // REFRACTIVE INDEX
    // The range is chosen to be up to ~10.7 eV because Sellmeier's equation
    // for fused silica is valid only in that range

    // The following values for the refractive index have been calculated
    // using Sellmeier's equation:
    //    n^2 - 1 = B_1 * \lambda^2 / (\lambda^2 - C_1) +
    //            + B_2 * \lambda^2 / (\lambda^2 - C_2) +
    //            + B_3 * \lambda^2 / (\lambda^2 - C_3),
    // with wavelength \lambda in micrometers and
    //    B_1 = 4.73E-1, B_2 = 6.31E-1, B_3 = 9.06E-1
    //    C_1 = 1.30E-2, C_2 = 4.13E-3, C_3 = 9.88E+1.
    G4double B_1 = 4.73e-1;
    G4double B_2 = 6.31e-1;
    G4double B_3 = 9.06e-1;
    G4double C_1 = 1.30e-2;
    G4double C_2 = 4.13e-3;
    G4double C_3 = 9.88e+1;

    std::vector<G4double> rIndex;
    for (int i=0; i<ri_entries; i++) {
      G4double lambda = CLHEP::h_Planck*CLHEP::c_light/ri_energy[i]*1000; // in micron
      G4double n2 = 1 + B_1*pow(lambda,2)/(pow(lambda,2)-C_1)
        + B_2*pow(lambda,2)/(pow(lambda,2)-C_2)
        + B_3*pow(lambda,2)/(pow(lambda,2)-C_3);
      rIndex.push_back(sqrt(n2));
      // G4cout << "* FusedSilica rIndex:  " << std::setw(5) << ri_energy[i]/eV
      //       << " eV -> " << rIndex[i] << G4endl;
    }
    mptsilica->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energysilica = {
      optPhotMinE_,  6.46499 * eV,
      6.54000 * eV,  6.59490 * eV,  6.64000 * eV,  6.72714 * eV,
      6.73828 * eV,  6.75000 * eV,  6.82104 * eV,  6.86000 * eV,
      6.88000 * eV,  6.89000 * eV,  7.00000 * eV,  7.01000 * eV,
      7.01797 * eV,  7.05000 * eV,  7.08000 * eV,  7.08482 * eV,
      7.30000 * eV,  7.36000 * eV,  7.40000 * eV,  7.48000 * eV,
      7.52000 * eV,  7.58000 * eV,  7.67440 * eV,  7.76000 * eV,
      7.89000 * eV,  7.93000 * eV,  8.00000 * eV,
      optPhotMaxE_
    };

    std::vector<G4double> absLength = {
      noAbsLength_, noAbsLength_,
      200.0 * cm,   200.0 * cm,  90.0 * cm,  45.0 * cm,
      45.0 * cm,    30.0 * cm,  24.0 * cm,  21.0 * cm,
      20.0 * cm,    19.0 * cm,  16.0 * cm,  14.0 * cm,
      13.0 * cm,     8.5 * cm,   8.0 * cm,   6.0 * cm,
       1.5 * cm,     1.2 * cm,   1.0 * cm,   .65 * cm,
        .4 * cm,     .37 * cm,   .32 * cm,   .28 * cm,
        .22 * cm,    .215 * cm,  .00005*cm,
      .00005* cm
    };
    mptsilica->AddProperty("ABSLENGTH", abs_energysilica, absLength);
	
	silica->SetMaterialPropertiesTable(mptsilica);
	

	
	
	
	
	//FakeGrid
	FakeGrid = new G4Material("FakeGrid", 2.953*g/cm3,1);
	FakeGrid->AddElement(nist->FindOrBuildElement("Xe"), 1);




	//GasXenon
	GXe = nist->FindOrBuildMaterial("G4_Xe");
	G4MaterialPropertiesTable *mptGXe = new G4MaterialPropertiesTable();
	
    // REFRACTIVE INDEX
    std::vector<G4double> ri_energyGas;
    for (int i=0; i<ri_entries; i++) {
      ri_energyGas.push_back(optPhotMinE_ + i * eWidth);
    }

	G4double pressure=2*bar;
    G4double densityGas = GXeDensity(pressure);
    std::vector<G4double> rIndexGas;
    for (int i=0; i<ri_entries; i++) {
      rIndexGas.push_back(XenonRefractiveIndex(ri_energyGas[i], densityGas));
      // G4cout << "* GXe rIndex:  " << std::setw(7)
      //        << ri_energy[i]/eV << " eV -> " << rIndex[i] << G4endl;
    }
    mptGXe->AddProperty("RINDEX", ri_energyGas, rIndexGas, ri_entries);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energyGas = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLengthGas  = {noAbsLength_, noAbsLength_};
    mptGXe->AddProperty("ABSLENGTH", abs_energyGas, absLengthGas);
	
	GXe->SetMaterialPropertiesTable(mptGXe);

	
}














G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
	constexpr G4double optPhotMinE_ = 6 * eV;
	constexpr G4double optPhotMaxE_ = 9  * eV;
	
	G4NistManager *nist = G4NistManager::Instance();

	
	//WORLD
	solidWorld = new G4Box("solidWorld",0.3*m,0.3*m,0.3*m);
	logicWorld = new G4LogicalVolume(solidWorld,lXe,"logicWorld");
	physWorld = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),logicWorld,"physWorld",0,false,0);	

	//ESFERA
	solidSphere = new G4Sphere("solidSphere",0.2*m,0.22*m,90*degree ,360*degree,90*degree ,180*degree);	
	logicSphere = new G4LogicalVolume(solidSphere,Steel,"logicalSphere");
	physSphere = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),logicSphere,"physSphere",logicWorld,false,0);
	
	
	
	
	G4double bigradius = 38*mm;
	G4double smallradius = 32*mm;
	
	//DETECTORES DIFERENTES CASOS//
	if(PMTR){
		// Dimensions
		G4double front_body_diam_=76. * mm;
		G4double front_body_length_ =38. * mm;
		G4double rear_body_diam_ =53. * mm;
		G4double rear_body_length_ =76. * mm;
		G4double body_thickness_ =.5 * mm;   
		G4double window_thickness_ =2. * mm;
		G4double photocathode_diam_ =64. * mm;
		G4double photocathode_thickness_ =.1 * mm;
		G4bool visibility_=1;
		G4int sd_depth_=-1;
		G4double binning_=100.*nanosecond;
		
		G4RotationMatrix* Rotation = new G4RotationMatrix();
		Rotation->rotateX(0*deg);
		Rotation->rotateY(180*deg);
		Rotation->rotateZ(0*deg);
	
		// PMT BODY //////////////////////////////////////////////////////
		G4Tubs* front_body_solid = new G4Tubs("FRONT_BODY", 0., front_body_diam_/2., front_body_length_/2.,0., 360*degree);
		G4Tubs* rear_body_solid = new G4Tubs("REAR_BODY", 0., rear_body_diam_/2., rear_body_length_/2.,0., 360*degree);

		// Union of the two volumes of the phototube body
		G4double z_transl = -front_body_length_/2. - rear_body_length_/2.;
		G4ThreeVector transl(0., 0., z_transl);
		G4UnionSolid* pmt_solid = new G4UnionSolid("PMT_R11410",front_body_solid,rear_body_solid,0,transl);


		G4LogicalVolume* pmt_logic = new G4LogicalVolume(pmt_solid, Kovar, "PMT_R11410");


		// PMT GAS  //////////////////////////////////////////////////////
		G4double front_body_gas_diam = front_body_diam_ - 2. * body_thickness_;
		G4double front_body_gas_length = front_body_length_ - body_thickness_;
		G4Tubs* front_body_gas_solid = new G4Tubs("FRONT_BODY_GAS", 0., front_body_gas_diam/2., front_body_gas_length/2., 0., 360*degree);

		G4double rear_body_gas_diam = rear_body_diam_ - 2. * body_thickness_;
		G4double rear_body_gas_length = rear_body_length_;
		G4Tubs* rear_body_gas_solid = new G4Tubs("REAR_BODY_GAS", 0., rear_body_gas_diam/2., rear_body_gas_length/2., 0., 360*degree);

		// Union of the two volumes of the phototube body gas
		G4double z_gas_transl = -front_body_gas_length/2. - rear_body_gas_length/2.;
		G4ThreeVector gas_transl(0., 0., z_gas_transl);
		G4UnionSolid* pmt_gas_solid = new G4UnionSolid("PMT_GAS", front_body_gas_solid,rear_body_gas_solid, 0, gas_transl);

		G4Material* pmt_gas_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
		pmt_gas_mat->SetMaterialPropertiesTable(new G4MaterialPropertiesTable());		
		G4LogicalVolume* pmt_gas_logic = new G4LogicalVolume(pmt_gas_solid, pmt_gas_mat, "PMT_GAS");

		G4double pmt_gas_posz = body_thickness_/2.;
		new G4PVPlacement(0, G4ThreeVector(0., 0., pmt_gas_posz), pmt_gas_logic,"PMT_GAS", pmt_logic, false, 0);


		// PMT WINDOW ////////////////////////////////////////////////////
		G4double  window_diam_ = front_body_gas_diam;
		G4Tubs* window_solid = new G4Tubs("PMT_WINDOW", 0, window_diam_/2., window_thickness_/2., 0., 360*degree);

		G4LogicalVolume* window_logic = new G4LogicalVolume(window_solid, silica, "PMT_WINDOW");

		G4double window_posz = front_body_gas_length/2. - window_thickness_/2.;
		new G4PVPlacement(0, G4ThreeVector(0.,0.,window_posz), window_logic,"PMT_WINDOW", pmt_gas_logic, false, 0);


		// PMT PHOTOCATHODE  /////////////////////////////////////////////
		G4Tubs* photocathode_solid = new G4Tubs("PMT_PHOTOCATHODE", 0, photocathode_diam_/2., photocathode_thickness_/2.,0., 360*degree);

		G4Material* aluminum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");	
		aluminum->SetMaterialPropertiesTable(new G4MaterialPropertiesTable());
		
		photocathode_logic = new G4LogicalVolume(photocathode_solid, aluminum, "PMT_PHOTOCATHODE");

		G4double photocathode_posz = window_posz - window_thickness_/2. - photocathode_thickness_/2.;
		new G4PVPlacement(0, G4ThreeVector(0., 0.,photocathode_posz), photocathode_logic,"PMT_PHOTOCATHODE", pmt_gas_logic, false, 0);


		// Optical properties
		G4OpticalSurface* pmt_opt_surf = GetPhotOptSurf();
		new G4LogicalSkinSurface("PMT_PHOTOCATHODE", photocathode_logic, pmt_opt_surf);	
		
		new G4PVPlacement(Rotation, G4ThreeVector(0., 0.,DetPos+19), pmt_logic,"PMT_physical", logicWorld, false, 0);
		}
	else if (Patron1){
		G4double r=12.8*mm;
		
		solidDetector = new G4Tubs("solidDetector",0.,8,0.001*m,0*degree ,360*degree);
		logicDetector = new G4LogicalVolume(solidDetector,lXe,"logicDetector");
		solidDetectorTapa = new G4Tubs("solidDetectorTapa",8,r,0.001*m,0*degree ,360*degree);
		logicDetectorTapa = new G4LogicalVolume(solidDetectorTapa,Steel,"logicDetectorTapa");
		
		physDetector  = new G4PVPlacement(0,G4ThreeVector(18.2,0.,DetPos+1),logicDetector,"physDetector",logicWorld,false,0);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(18.2,0.,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,0);
		
		physDetector  = new G4PVPlacement(0,G4ThreeVector(-18.2,0.,DetPos+1),logicDetector,"physDetector",logicWorld,false,1);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(-18.2,0.,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,1);
		
		physDetector  = new G4PVPlacement(0,G4ThreeVector(0.,18.2,DetPos+1),logicDetector,"physDetector",logicWorld,false,2);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(0.,18.2,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,2);
				
		physDetector  = new G4PVPlacement(0,G4ThreeVector(0.,-18.2,DetPos+1),logicDetector,"physDetector",logicWorld,false,3);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(0.,-18.2,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,3);
				
		
		//Eficiencia
		G4MaterialPropertiesTable* star_mpt = new G4MaterialPropertiesTable();
		G4double starenergy[2]={6*eV,9*eV};
		G4double starefficiency[2]={0.95,0.95};
		star_mpt->AddProperty("REFLECTIVITY", starenergy, starefficiency, 2);

		G4OpticalSurface* star_surf = new G4OpticalSurface("PHOTOCATHODE", unified, polished, dielectric_metal);
		star_surf->SetMaterialPropertiesTable(star_mpt);
		new G4LogicalSkinSurface("PMT_PHOTOCATHODE", logicDetector, star_surf);	
		
		}
	else if (Patron2){
		bigradius = 38.6*mm;
		G4double r=12.8*mm;
		
		solidDetector = new G4Tubs("solidDetector",0.,8,0.001*m,0*degree ,360*degree);
		logicDetector = new G4LogicalVolume(solidDetector,lXe,"logicDetector");
		solidDetectorTapa = new G4Tubs("solidDetectorTapa",8,r,0.001*m,0*degree ,360*degree);
		logicDetectorTapa = new G4LogicalVolume(solidDetectorTapa,Steel,"logicDetectorTapa");
		
		physDetector  = new G4PVPlacement(0,G4ThreeVector(0.,0.,DetPos+1),logicDetector,"physDetector",logicWorld,false,0);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(0.,0.,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,0);
		
		physDetector  = new G4PVPlacement(0,G4ThreeVector(-2*r,0.,DetPos+1),logicDetector,"physDetector",logicWorld,false,1);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(-2*r,0.,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,1);
		
		physDetector  = new G4PVPlacement(0,G4ThreeVector(2*r,0.,DetPos+1),logicDetector,"physDetector",logicWorld,false,2);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(2*r,0.,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,2);
				
		physDetector  = new G4PVPlacement(0,G4ThreeVector(r,r*sqrt(3),DetPos+1),logicDetector,"physDetector",logicWorld,false,3);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(r,r*sqrt(3),DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,3);
				
		physDetector  = new G4PVPlacement(0,G4ThreeVector(r,-r*sqrt(3),DetPos+1),logicDetector,"physDetector",logicWorld,false,4);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(r,-r*sqrt(3),DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,4);
				
		physDetector  = new G4PVPlacement(0,G4ThreeVector(-r,r*sqrt(3),DetPos+1),logicDetector,"physDetector",logicWorld,false,5);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(-r,r*sqrt(3),DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,5);
				
		physDetector  = new G4PVPlacement(0,G4ThreeVector(-r,-r*sqrt(3),DetPos+1),logicDetector,"physDetector",logicWorld,false,6);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(-r,-r*sqrt(3),DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,6);
				
		
		//Eficiencia
		G4MaterialPropertiesTable* star_mpt = new G4MaterialPropertiesTable();
		G4double starenergy[2]={6*eV,9*eV};
		G4double starefficiency[2]={0.95,0.95};
		star_mpt->AddProperty("REFLECTIVITY", starenergy, starefficiency, 2);

		G4OpticalSurface* star_surf = new G4OpticalSurface("PHOTOCATHODE", unified, polished, dielectric_metal);
		star_surf->SetMaterialPropertiesTable(star_mpt);
		new G4LogicalSkinSurface("PMT_PHOTOCATHODE", logicDetector, star_surf);
		
		}
	else if (Patron3){
		bigradius = 64*mm;
		G4double r=12.8*mm;
		
		solidDetector = new G4Tubs("solidDetector",0.,8,0.001*m,0*degree ,360*degree);
		logicDetector = new G4LogicalVolume(solidDetector,lXe,"logicDetector");
		solidDetectorTapa = new G4Tubs("solidDetectorTapa",8,r,0.001*m,0*degree ,360*degree);
		logicDetectorTapa = new G4LogicalVolume(solidDetectorTapa,Steel,"logicDetectorTapa");
		
		//0
		physDetector  = new G4PVPlacement(0,G4ThreeVector(0.,0.,DetPos+1),logicDetector,"physDetector",logicWorld,false,0);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(0.,0.,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,0);
		//1
		physDetector  = new G4PVPlacement(0,G4ThreeVector(-2*r,0.,DetPos+1),logicDetector,"physDetector",logicWorld,false,1);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(-2*r,0.,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,1);
		//2
		physDetector  = new G4PVPlacement(0,G4ThreeVector(2*r,0.,DetPos+1),logicDetector,"physDetector",logicWorld,false,2);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(2*r,0.,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,2);
		//3		
		physDetector  = new G4PVPlacement(0,G4ThreeVector(r,r*sqrt(3),DetPos+1),logicDetector,"physDetector",logicWorld,false,3);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(r,r*sqrt(3),DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,3);
		//4		
		physDetector  = new G4PVPlacement(0,G4ThreeVector(r,-r*sqrt(3),DetPos+1),logicDetector,"physDetector",logicWorld,false,4);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(r,-r*sqrt(3),DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,4);
		//5		
		physDetector  = new G4PVPlacement(0,G4ThreeVector(-r,r*sqrt(3),DetPos+1),logicDetector,"physDetector",logicWorld,false,5);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(-r,r*sqrt(3),DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,5);
		//6		
		physDetector  = new G4PVPlacement(0,G4ThreeVector(-r,-r*sqrt(3),DetPos+1),logicDetector,"physDetector",logicWorld,false,6);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(-r,-r*sqrt(3),DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,6);
		
		//7		
		physDetector  = new G4PVPlacement(0,G4ThreeVector(-4*r,0.,DetPos+1),logicDetector,"physDetector",logicWorld,false,7);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(-4*r,0.,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,7);
		//8
		physDetector  = new G4PVPlacement(0,G4ThreeVector(4*r,0.,DetPos+1),logicDetector,"physDetector",logicWorld,false,8);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(4*r,0.,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,8);
		//9		
		physDetector  = new G4PVPlacement(0,G4ThreeVector(3*r,-sqrt(3)*r,DetPos+1),logicDetector,"physDetector",logicWorld,false,9);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(3*r,-sqrt(3)*r,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,9);
		//10	
		physDetector  = new G4PVPlacement(0,G4ThreeVector(3*r,sqrt(3)*r,DetPos+1),logicDetector,"physDetector",logicWorld,false,10);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(3*r,sqrt(3)*r,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,10);
		//11
		physDetector  = new G4PVPlacement(0,G4ThreeVector(-3*r,-sqrt(3)*r,DetPos+1),logicDetector,"physDetector",logicWorld,false,11);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(-3*r,-sqrt(3)*r,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,11);
		//12
		physDetector  = new G4PVPlacement(0,G4ThreeVector(-3*r,sqrt(3)*r,DetPos+1),logicDetector,"physDetector",logicWorld,false,12);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(-3*r,sqrt(3)*r,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,12);
		//13	
		physDetector  = new G4PVPlacement(0,G4ThreeVector(0.,2*sqrt(3)*r,DetPos+1),logicDetector,"physDetector",logicWorld,false,13);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(0.,2*sqrt(3)*r,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,13);
		//14		
		physDetector  = new G4PVPlacement(0,G4ThreeVector(0.,-2*sqrt(3)*r,DetPos+1),logicDetector,"physDetector",logicWorld,false,14);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(0.,-2*sqrt(3)*r,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,14);
		//15		
		physDetector  = new G4PVPlacement(0,G4ThreeVector(2*r,2*sqrt(3)*r,DetPos+1),logicDetector,"physDetector",logicWorld,false,15);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(2*r,2*sqrt(3)*r,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,15);
		//16		
		physDetector  = new G4PVPlacement(0,G4ThreeVector(-2*r,2*sqrt(3)*r,DetPos+1),logicDetector,"physDetector",logicWorld,false,16);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(-2*r,2*sqrt(3)*r,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,16);
		//17		
		physDetector  = new G4PVPlacement(0,G4ThreeVector(-2*r,-2*sqrt(3)*r,DetPos+1),logicDetector,"physDetector",logicWorld,false,17);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(-2*r,-2*sqrt(3)*r,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,17);
		//18		
		physDetector  = new G4PVPlacement(0,G4ThreeVector(2*r,-2*sqrt(3)*r,DetPos+1),logicDetector,"physDetector",logicWorld,false,18);
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(2*r,-2*sqrt(3)*r,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,18);
				
		
		//Eficiencia
		G4MaterialPropertiesTable* star_mpt = new G4MaterialPropertiesTable();
		G4double starenergy[2]={6*eV,9*eV};
		G4double starefficiency[2]={0.95,0.95};
		star_mpt->AddProperty("REFLECTIVITY", starenergy, starefficiency, 2);

		G4OpticalSurface* star_surf = new G4OpticalSurface("PHOTOCATHODE", unified, polished, dielectric_metal);
		star_surf->SetMaterialPropertiesTable(star_mpt);
		new G4LogicalSkinSurface("PMT_PHOTOCATHODE", logicDetector, star_surf);
		
		}
	else{
		//Para PMT pequeño
		//bigradius = 12.8*mm;
		//smallradius = 8*mm;
	
		//DETECTOR GENERIC SENSOR
		solidDetector = new G4Tubs("solidDetector",0.,smallradius,0.001*m,0*degree ,360*degree);
		logicDetector = new G4LogicalVolume(solidDetector,lXe,"logicDetector");
		physDetector  = new G4PVPlacement(0,G4ThreeVector(0.,0.,DetPos+1),logicDetector,"physDetector",logicWorld,false,0);
		
		solidDetectorTapa = new G4Tubs("solidDetectorTapa",smallradius,bigradius,0.001*m,0*degree ,360*degree);
		logicDetectorTapa = new G4LogicalVolume(solidDetectorTapa,Steel,"logicDetectorTapa");
		physDetectorTapa  = new G4PVPlacement(0,G4ThreeVector(0.,0.,DetPos+1),logicDetectorTapa,"physDetectorTapa",logicWorld,false,0);
		}
	
	

	//TUBO1
	solidTub1 = new G4Tubs("solidTub1",bigradius,bigradius+1*mm,0.1*m,0*degree ,360*degree);	
	logicTub1 = new G4LogicalVolume(solidTub1,Steel,"logicalTub1");
	physTub1 = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.08*m),logicTub1,"physTub1",logicWorld,false,0);
		
	//TUBO2 (TAPA DE ESFERA)
	solidTub2 = new G4Tubs("solidTub2",bigradius+1*mm,0.22*m,0.01*m,0*degree ,360*degree);	
	logicTub2 = new G4LogicalVolume(solidTub2,Steel,"logicalTub2");
	physTub2 = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.01*m),logicTub2,"physTub2",logicWorld,false,0);
				
	//REFLECTIVITY STEEL
	G4OpticalSurface *SteelSurf = new G4OpticalSurface("SteelSurf");
	new G4LogicalBorderSurface("SteelSurface",physWorld,physSphere,SteelSurf);
	new G4LogicalBorderSurface("SteelSurface",physWorld,physTub1,SteelSurf);
	new G4LogicalBorderSurface("SteelSurface",physWorld,physTub2,SteelSurf);

	SteelSurf->SetModel(glisur);
	SteelSurf->SetType(dielectric_metal);
	SteelSurf->SetFinish(polished);
 
	G4MaterialPropertiesTable *Smpt = new G4MaterialPropertiesTable();
	
	G4double energyreflx[2]={6.2*eV,8.2*eV};
	G4double Reflect[2]={0.0,0.0};
	Smpt->AddProperty("REFLECTIVITY", energyreflx,Reflect,2);
	SteelSurf->SetMaterialPropertiesTable(Smpt);

	
	
	
	
	
	
	
	//TEFLON CASOS
	if(CaseA){
		//ESFERA TEFLON
		solidSphereTeflon = new G4Sphere("solidSphereTeflon",0.199*m,0.2*m,90*degree ,360*degree,90*degree ,180*degree);	
		logicSphereTeflon = new G4LogicalVolume(solidSphereTeflon,Teflon,"logicalSphereTeflon");
		physSphereTeflon  = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),logicSphereTeflon,"physSphereTeflon",logicWorld,false,0);
	}
	
	if(CaseB){
		if (!Gas) {
		//TUBO2 (TAPA) TEFLON
		solidTub2Teflon = new G4Tubs("solidTub2Teflon",bigradius+1*mm,0.198*m,0.001*m,0*degree ,360*degree);	
		logicTub2Teflon = new G4LogicalVolume(solidTub2Teflon,Teflon,"logicalTub2Teflon");
		physTub2Teflon  = new G4PVPlacement(0,G4ThreeVector(0.,0.,-0.001*m),logicTub2Teflon,"physTub2Teflon",logicWorld,false,0);
		}
	}	

	if(CaseC){
		//TUBO1 TEFLON INTERIOR
		solidTub1TeflonInt = new G4Tubs("solidTub1TeflonExt",bigradius-0.1*mm,bigradius,DetPos/2.0+1*cm,0*degree ,360*degree);	
		logicTub1TeflonInt = new G4LogicalVolume(solidTub1TeflonInt,Teflon,"logicTub1TeflonInt");
		physTub1TeflonInt  = new G4PVPlacement(0,G4ThreeVector(0.,0.,-1*cm+DetPos/2.0),logicTub1TeflonInt,"physTub1TeflonInt",logicWorld,false,0);
	}	
		
	if(CaseD){
		//TUBO1 TEFLON EXTERIOR
		solidTub1TeflonExt = new G4Tubs("solidTub1TeflonExt",bigradius+1*mm,bigradius+2*mm,0.009*m,0*degree ,360*degree);	
		logicTub1TeflonExt = new G4LogicalVolume(solidTub1TeflonExt,Teflon,"logicTub1TeflonExt");
		physTub1TeflonExt  = new G4PVPlacement(0,G4ThreeVector(0.,0.,-0.011*m),logicTub1TeflonExt,"physTub1TeflonExt",logicWorld,false,0);
	}
	
	//REFLECTIVITY TEFLON
	G4OpticalSurface *TeflonSurf = new G4OpticalSurface("TeflonSurf");
	new G4LogicalBorderSurface("TeflonSurface",physWorld,physSphereTeflon,TeflonSurf);
	new G4LogicalBorderSurface("TeflonSurface",physWorld,physTub2Teflon,TeflonSurf);
	new G4LogicalBorderSurface("TeflonSurface",physWorld,physTub1TeflonInt,TeflonSurf);
	new G4LogicalBorderSurface("TeflonSurface",physWorld,physTub1TeflonExt,TeflonSurf);

	TeflonSurf->SetModel(glisur);
	TeflonSurf->SetType(dielectric_metal);
	TeflonSurf->SetFinish(polished);
 
	G4MaterialPropertiesTable *Teflonmpt = new G4MaterialPropertiesTable();
	
	// REFLECTIVITY IN LXE (from https://link.springer.com/content/pdf/10.1140/epjc/s10052-020-7800-6.pdf)
	std::vector<G4double> ENERGIESTeflon = {optPhotMinE_,optPhotMaxE_};
	std::vector<G4double> REFLECTIVITYTeflon = {0.98, 0.98};
	Teflonmpt->AddProperty("REFLECTIVITY", ENERGIESTeflon, REFLECTIVITYTeflon);

	// REFLEXION BEHAVIOR
	// Specular reflection about the normal to a microfacet.
	// Such a vector is chosen according to a gaussian distribution with
	// sigma = SigmaAlhpa (in rad) and centered in the average normal.
	std::vector<G4double> specularlobe  = {0., 0.};
	// specular reflection about the average normal
	std::vector<G4double> specularspike = {0., 0.};
	// 180 degrees reflection.
	std::vector<G4double> backscatter   = {0., 0.};
	// 1 - the sum of these three last parameters is the percentage of Lambertian reflection
	Teflonmpt->AddProperty("SPECULARLOBECONSTANT", ENERGIESTeflon, specularlobe);
	Teflonmpt->AddProperty("SPECULARSPIKECONSTANT",ENERGIESTeflon, specularspike);
	Teflonmpt->AddProperty("BACKSCATTERCONSTANT",  ENERGIESTeflon, backscatter);

	TeflonSurf->SetMaterialPropertiesTable(Teflonmpt);
	
	
	
	
	
	

	
	
	
	
	
	/////FAKEGRID PROPERTIES//////
	G4MaterialPropertiesTable *mptFakeGrid = new G4MaterialPropertiesTable();
	
	const G4int ri_entries = 200;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;
	
	std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }
    G4double density = 2.953 * g/cm3;
    std::vector<G4double> ri_index;
    for (G4int i=0; i<ri_entries; i++) {
      ri_index.push_back(XenonRefractiveIndex(ri_energy[i], density));
    }
	
	//REFRACTIVE INDEX
	mptFakeGrid->AddProperty("RINDEX",ri_energy,ri_index);
	//ABSORTION LENGTH

	if (transparency != 1.0) {
		G4double thickness=0.1*mm;
		std::vector<G4double> abs_energyGrid = {optPhotMinE_, optPhotMaxE_};
		std::vector<G4double> abs_lengthGrid = {-thickness/log(transparency),-thickness/log(transparency)};
		mptFakeGrid->AddProperty("ABSLENGTH", abs_energyGrid, abs_lengthGrid);
	}
	FakeGrid->SetMaterialPropertiesTable(mptFakeGrid);	
	
	//MESH
	solidMesh = new G4Tubs("solidMesh",0.,bigradius-0.1*mm,0.05*mm,0*degree ,360*degree);
	logicMesh = new G4LogicalVolume(solidMesh,FakeGrid,"logicMesh");
	physMesh  = new G4PVPlacement(0,G4ThreeVector(0.,0.,zGun),logicMesh,"physMesh",logicWorld,false,0);
	physMesh  = new G4PVPlacement(0,G4ThreeVector(0.,0.,zGun+5-0.05),logicMesh,"physMesh",logicWorld,false,1);
	physMesh  = new G4PVPlacement(0,G4ThreeVector(0.,0.,zGun-5+0.05),logicMesh,"physMesh",logicWorld,false,2);

	
	
	
	

	//XENON GAS EN LA TAPA DE ESFERA
	if(Gas){
		solidGasfilm = new G4Tubs("solidGasfilm",bigradius+2*mm,0.195*m,0.125*cm,0*degree ,360*degree);	
		logicGasfilm = new G4LogicalVolume(solidGasfilm,GXe,"logicGasfilm");
		physGasfilm  = new G4PVPlacement(0,G4ThreeVector(0.,0.,-0.125*cm),logicGasfilm,"physGasfilm",logicWorld,false,0);
	}		

	return physWorld;		
}









G4OpticalSurface* MyDetectorConstruction::GetPhotOptSurf()
  {
    const G4int entries = 57;
    G4double ENERGIES[entries] =
      {	2.06640321682 *eV, 2.10142700016 *eV,  2.13765850016 *eV,  2.17516128086 *eV,
	2.21400344659 *eV, 2.25425805471 *eV,  2.29600357425 *eV,  2.3393243964 *eV,
	2.38431140402 *eV, 2.43106260802 *eV,  2.47968386018 *eV,  2.53028965325 *eV,
	2.58300402103 *eV, 2.63796155339 *eV,  2.69530854368 *eV,  2.75520428909 *eV,
	2.81782256839 *eV, 2.8833533258 *eV,   2.95200459546 *eV,  3.02400470754 *eV,
	3.09960482523 *eV, 3.17908187203 *eV,  3.2627419213 *eV,   3.35092413538 *eV,
	3.44400536137 *eV, 3.54240551455 *eV,  3.64659391204 *eV,  3.75709675786 *eV,
	3.87450603154 *eV, 3.99949009707 *eV,  4.13280643364 *eV,  4.27531700032 *eV,
	4.42800689319 *eV, 4.59200714849 *eV,  4.76862280805 *eV,  4.95936772037 *eV,
	5.16600804205 *eV, 5.39061708736 *eV,  5.63564513678 *eV,  5.90400919092 *eV,
	6.19920965046 *eV, 6.358163744063561 *eV, 6.525483842591549 *eV, 6.7018482707697 *eV,
	6.888010722735524 *eV, 7.084811029099397 *eV, 7.293187824072908 *eV, 7.314701652462505 *eV,
	7.336342781611801 *eV, 7.358112344761985 *eV, 7.380011488645204 *eV, 7.402041373685937 *eV,
	7.424203174205955 *eV, 7.446498078633 *eV, 7.4689272897132195 *eV, 7.491492024727459 *eV,
	7.51419351571 *eV};

    G4double EFFICIENCY[entries] =
      { 0.0530,	0.0625, 0.0720, 0.0850,
	0.1050, 0.1190, 0.1335, 0.1550,
	0.1770, 0.1970, 0.2100, 0.2200,
	0.2300,	0.2430, 0.2580, 0.2770,
	0.2920,	0.3050, 0.3150, 0.3270,
	0.3320, 0.3400, 0.3480, 0.3500,
	0.3530,	0.3600, 0.3680, 0.3650,
	0.3640, 0.3640, 0.3560, 0.3420,

	0.3280, 0.3180, 0.3050, 0.2980,
	0.2920,	0.2900, 0.2920, 0.2945,
	0.3100,	0.3280, 0.3560, 0.3880,
	0.3920,	0.3900, 0.4040, 0.3930,
	0.3700,	0.3500, 0.3300, 0.3150,
	0.2950,	0.2750, 0.2550, 0.2450,
	0.2400 };
	
	
    G4double REFLECTIVITY[entries] =
     { 0.0,	0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0,

	0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0,
	0.0 };
   

    G4MaterialPropertiesTable* phcath_mpt = new G4MaterialPropertiesTable();
    //phcath_mpt->AddProperty("EFFICIENCY", ENERGIES, EFFICIENCY, entries);
    phcath_mpt->AddProperty("REFLECTIVITY", ENERGIES, EFFICIENCY, entries);

    G4OpticalSurface* opt_surf = new G4OpticalSurface("PHOTOCATHODE", unified, polished, dielectric_metal);
    opt_surf->SetMaterialPropertiesTable(phcath_mpt);

    return opt_surf;
  }








void MyDetectorConstruction::ConstructSDandField()
{
	MySensitiveDetector *sensDet=new MySensitiveDetector("SensitiveDetector");
	
	if(PMTR){
		photocathode_logic->SetSensitiveDetector(sensDet);
	}
	else{
		logicDetector->SetSensitiveDetector(sensDet);

	}

}




