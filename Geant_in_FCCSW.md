Simulation with Geant4 in FCCSW
====

Full simulation uses the detailed detector description and simulates the particles passage through it taking into account all the physics processes they may encounter. It is very time and CPU consuming.

Therefore for many tasks, especially in the early stage of detector design, the fast simulation is used. It takes a less detailed description of the detector and does not simulate every particle step-by-step. Instead, it simulates the overall response of the (particular) detector in a parametric way.

Generated particles are transported inside the detector and they are smeared taking into account the resolutions and efficiency (their 4-momentum and/or position). Those smeared particles may be analysed and treated as reconstructed particles, even though no hits were produced and no reconstruction was performed. All the detector effects (both physics processes that may encounter and detector resolutions) come from the smearing process.

The resolutions used in the smearing may come arbitrary from our knowledge of the detectors. In that case one applies a Gaussian smearing with a given standard deviation. That approach may be also used by the physicists to test how the detector resolution affect the results. That smearing is currently implemented in FCCSW.

More complex approach involves construction of the tables with detector resolutions (pseudorapidity/momentum/particle dependent). They are calculated from a small (relatively) sample of full simulations of single-particle events. Single-particle events simplify the reconstruction process (they don't involve the pattern recognition etc.). Such resolutions are valid for tested detectors hence they may be used for smearing the particles with a better accuracy. Implementation of this approach is still in progress.

Both full and fast simulation can be performed in FCCSW using Geant4. Since the same tools are used for both of them, each simulation may be an interplay of both with full simulation performed in some volumes and fast simulation in others.


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
    - path to the XML file with geometry is specified as a property "detector"

                  from Configurables import GeoSvc
                  geoservice = GeoSvc("GeoSvc",
                                       detector='file:DetectorDescription/Detectors/compact/ParametricSimTracker.xml',
                                       OutputLevel = VERBOSE)

  * simulation with Geant
    - simulation type is specified in a property "simtype"

                 from Configurables import Geant4Simulation
                 geant4simulation = Geant4Simulation("Geant4Simulation", simtype="full")
                 geant4simulation.Inputs.genparticles.Path="all_genparticles"
                 geant4simulation.Outputs.particles.Path = "recparticles"
                 geant4simulation.Outputs.particleassociation.Path = "particleMCparticle"

  * saving the output

            from Configurables import AlbersWrite, AlbersOutput
            out = AlbersOutput("out", OutputLevel=DEBUG)
            out.outputCommands = ["keep *"]

  * GAUDI's main component
      - stating which algorithms (TopAlg) should be initialised, executed and finalised in the event loop
      - specifying how many events from the input file should be processed (EvtMax)
      - stating which services (ExtSvc) should be created at initialisation

                from Configurables import ApplicationMgr
                ApplicationMgr( TopAlg = [reader, hepmc_converter, geant4simulation, out],
                                EvtSel = 'NONE',
                                EvtMax   = 1,
                                ExtSvc = [albersevent, geoservice],
                                OutputLevel=DEBUG
                                )

2. Fast simulation
----

To run fast simulation:

    ./run gaudirun.py config/geant_fastsim.py

The differences between configuration file of the fast simulation (config/geant_fastsim.py) and the full simulation:

  * simulation in Geant with different properties
    - simulation type is specified in a property "simtype"
    - name of the tool used for smearing is specified in a property "smearingtoolname"

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

Geant4Simulation has two properties:

    simtype
    smearingtoolname

Property simtype describes which simulation type the user wants to perform: "fast" or "full". The default is "full".
Property smearingtoolname is valid only for fast simulation. It is a name of a class derived from ISmearingTool (GeantSim/GeantComponents/GeantComponents/ISmearingTool.h). It describes how the momentum/energy is smeared.

### 3.2. Geometry

Geometry is provided by DD4hep via GAUDI's service GeoSvc.

DD4hep translates the geometry into Geant4 geometry. That translation is done based on xml file as well as on the corresponding constructors:

       DetectorDescription/Detectors/compact/ParametricSimTracker.xml
       DetectorDescription/Detectors/src/ParametricSimTracker_geo.cxx

Path to the XML file with detector description is a property of service GeoSvc.

  * In fast simulation user wants a specific behaviour in certain geometry volumes. That specific behaviour is described in classes derived from G4VFastSimulationModel. Geant will not perform normal transportation inside volumes with fast simulation models attached (providing that particle triggers that model).

  The fast simulation model needs to be attached to a G4Region object. That G4Region can contain one or many logical volumes (parts of the detector). Logical volumes are created by DD4hep translation.

  Currently a fast simulation model FastSimModelTracker (GeantSim/GeantFast/GeantFast/FastSimModelTracker.h) will be created and attached automatically for every G4LogicalVolume that contains 'tracker' in its name in DD4hep geometry constructor (with the intermediate step of creating the G4Region). Work is in progress to change that to be easier to modify (in particular to be set from xml file).

### 3.3. Physics List

Geant simulation requires the particles and processes definitions - the so-called physics list. The currently used physics list is FTFP_BERT, which is recommended by Geant4 for HEP.

  * additionally for fast simulation:
     - "Coupled Transportation" is used which allows to invoke G4PathFinder which propagates the particle in the magnetic field.
     It is used within fast simulation model for tracker to compute the exit position of a particle from the volume (taking the momentum from the entrance to the tracker volume).
     - FastSimPhysics is registered as an additional process (GeantSim/GeantFast/GeantFast/FastSimPhysics.h). It attaches the fast simulation manager process to ALL the particles. That means that along with the standard processes such as transportation, multiple scattering etc. the particle can encounter 'fast simulation' process. That happens if a particle enters a volume with fast simulation model attached and if that particle fulfils trigger conditions (is charged in case of FastSimModelTracker).

     Details may be found at
                http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/ch05s02.html#sect.PhysProc.Param


### 3.4. User Actions

Geant allows to specify actions that will be invoked at different stages of event processing.

TrackingAction (GeantSim/GeantGeneral/GeantGeneral/TrackingAction.h) saves to the EDM the information on simulated particle (at the end of processing the particle with Geant).

To be able to compare generated particles (EDM MCParticle) with reconstructed particles (EDM Particle), one needs to create the association between those objects. That association is created while the event is translated from EDM to G4Event. A handle to MCParticle and Particle (the latter at that time empty) is stored in ParticleInformation (GeantSim/GeantGeneral/GeantGeneral/ParticleInformation.h). It is an object that Geant allows to be created for each G4PrimaryParticle. At the end of tracking, the information about the simulated particle (momentum, vertex position) is saved to the EDM Particle using the handle from ParticleInformation.

The user actions are the last (and optional) element of G4RunManager initialisation.


### 3.5. Event Processing

For each execution of the algorithm an EDM MCParticleCollection is translated into G4Event using the method G4Simulation::EDM2G4().

At the same time the ParticleCollection is created so that the association between MCParticle and Particle in EDM can be preserved. ParticleCollection will be the collection of the so-called 'reconstructed' particles. 'Reconstructed' in this case means the particles have undergone simulation process and the resulting (final) changes to the momentum are taken into account. At this stage it does not involve any reconstruction process.


4. Smearing
----

In case of running the simulation with simtype=='fast', user can define a smearing tool that should be used.

Smearing is performed using the GAUDI tool derived from ISmearingTool. User may specify how to smear a vector (CLHEP::ThreeVector) or a scalar (double). Currently there is only one implementation of a smearing tool called SimpleSmear (GeantSim/GeantComponents/src/SimpleSmear.h). The momentum/energy are smeared by multiplying Gaussian distribution with the mean $\mu=1$ and the standard deviation $\sigma$ as specified in the tool properties (default: 0.01).
