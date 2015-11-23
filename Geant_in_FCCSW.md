Simulation with Geant4 in FCCSW
====
1. Full simulation
----

To run full simulation:

    ./run gaudirun.py config/geant_fullsim.py

The configuration file (config/geant_fullsim.py) contains:
  * reading an event saved in a HepMC file

            reader = HepMCReader("Reader", Filename="example_MyPythia.dat")
            reader.Outputs.hepmc.Path = "hepmc"

  * translating a HepMC event to the EDM
   
            from Configurables import HepMCConverter
            hepmc_converter = HepMCConverter("Converter")
            hepmc_converter.Inputs.hepmc.Path="hepmc"
            hepmc_converter.Outputs.genparticles.Path="all_genparticles"
            hepmc_converter.Outputs.genvertices.Path="all_genvertices"

  * construction of geometry using DD4hep

            from Configurables import GeoSvc
            geoservice = GeoSvc("GeoSvc", OutputLevel = VERBOSE)

  * simulation with Geant

            from Configurables import Geant4Simulation
            geant4simulation = Geant4Simulation("Geant4Simulation", simtype="full")
            geant4simulation.Inputs.genparticles.Path="all_genparticles"
            geant4simulation.Outputs.particles.Path = "recparticles"
            geant4simulation.Outputs.particleassociation.Path = "particleMCparticle"

  * saving the output

            from Configurables import AlbersWrite, AlbersOutput
            out = AlbersOutput("out", OutputLevel=DEBUG)
            out.outputCommands = ["keep *"]



2. Fast simulation
----

To run fast simulation:

    ./run gaudirun.py config/geant_fastsim.py
   
The differences between configuration file of the fast simulation (config/geant_fastsim.py) and the full simulation:

  * simulation in Geant with different properties

            geant4simulation = Geant4Simulation("Geant4Simulation", simtype="fast",
                                                smearingtoolname = "SimpleSmear")

  * additional tool for smearing
    
            from Configurables import SimpleSmear
            smear = SimpleSmear("SimpleSmear", sigma = 0.015)
            geant4simulation.addTool(smear)

3. Main algorithm: Geant4Simulation
----

(Simulation/Simulation/Geant4Simulation.h)

Main algorithm for simulation with Geant4 Geant4Simulation derives from G4RunManager. Own implementation of G4RunManager was necessary in order to leave the event flow to GAUDI.


### 3.1. Properties

Gean4Simulation has two properties:

    simtype
    smearingtoolname

Propery simtype describes which simulation type the user wants to perform: "fast" or "full". The default is "full".
Propery smearingtoolname is valid only for fast simulation. It is a name of a class derived from ISmearingTool (GeantSim/GeantComponents/GeantComponents/ISmearingTool.h). It describes how the momentum/energy is smeared.

### 3.2. Geometry

Geometry is provided by DD4hep via GAUDI's service GeoSvc.

  * In fast simulation user wants a specific behaviour in certain geometry volumes. That specific behaviour is described in classes derived from G4VFastSimulationModel. Geant will not perform normal transportation inside volumes with fast simulation models attached (providing that particle triggers that model).

  The fast simulation model needs to be attached to a G4Region object. That G4Region can contain one or many logical volumes (parts of the detector).

  Currently a fast simulation model FastSimModelTracker (GeantSim/GeantFast/GeantFast/FastSimModelTracker.h) will be created and attached automatically for every G4LogicalVolume that contains 'tracker' in its name in DD4hep geometry constructor (with the intermidiate step of creating the G4Region). Work in progress to change that to be more customizable (in particular set from xml file).

### 3.3. Physics List

Geant simulation requires the particles and processes definitions - the so-called physics list. The currently used physics list is FTFP_BERT, which is recommended by Geant4 for HEP.

  * additionaly for fast simulation:
     - "Coupled Transportation" is used which allows to invoke G4PathFinder which propagates the particle in the magnetic field.
     It is used within fast simulation model for tracker to compute the exit position of a particle from the volume (taking the momentum from the entrance to the tracker volume).
     - FastSimPhysics is registered as an additional process (GeantSim/GeantFast/GeantFast/FastSimPhysics.h). It attaches the fast simulation manager process to ALL the particles. That means that along with the standard processes such as transportation, multiple scattering etc. the particle can encounter 'fast simulation' process. That happens if a particle enters a volume with fast simulation model attached and if that particle fullfils trigger conditions (is charged in case of FastSimModelTracker).

     Details may be found at
                http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/ch05s02.html#sect.PhysProc.Param


### 3.4. User Actions

Geant allows to specify actions that will be invoked at different stages of event processing.

TrackingAction (GeantSim/GeantGeneral/GeantGeneral/TrackingAction.h) saves to the EDM the information on simulated particle (at the end of processing the particle with Geant).

In order to be able to preserve the association between the MCParticle and Particle in EDM, the association is created while the event is translated from EDM to G4Event. A handle to MCParticle and associated Particle (at that time empty) is stored in ParticleInformation (GeantSim/GeantGeneral/GeantGeneral/ParticleInformation.h). It is an object that Geant allows to be created for each g4particle object. At the end of tracking, the information about the simulated particle (momentum, vertex position) is saved in the (EDM) Particle using the handle from ParticleInformation.

The user actions are the last (and optional) element of G4RunManager initialization.


### 3.5. Event Processing

For each execution of the algorithm an EDM MCParticleCollection is translated into G4Event using the method G4Simulation::EDM2G4().

At the same time the ParticleCollection is created so that the association between MCParticle and Particle in EDM can be preserved. ParticleCollection will be the collection of the so-called 'reconstructed' particles. 'Reconstructed' in this case means the particles have undergone simulation process and the resulting changes to the momentum are taken into account. It does not involve any reconstruction process.

