Fast Simulation with Geant4 in FCCSW
====

Instruction on how to use Geant4 fast simulation within FCCSW (GAUDI framework).
It is an addition to the instruction on Geant4 usage [see](Geant4_fullsim.md).

## Contents
* [Ovierview](#1-overview)
  * [Fast simulation](#fast-simulation)
* [Example of fast sim configuration](#2-example)
* [Geant configuration via GAUDI service `G4SimSvc`](#3-geant-configuration-via-gaudi-service-g4simsvc)
  * [Geometry construction](#3-1-geometry-construction)
  * [Physics list](#3-2-physics-list)
  * [User actions](#3-3-user-actions)
* [Simulation in GAUDI algorithm `G4SimAlg`](#4-simulation-in-gaudi-algorithm-g4simalg)
  * [Events](#4-1-event-processing)
  * [Output](#4-2-output)



## How to



## 1. Overview

For details on GAUDI and Geant see [full simulation overview](#Geant4_fullsim.md#overview)

### Fast simulation

Full simulation uses the detailed detector description and simulates the particles passage through it, taking into account all the physics processes they may encounter. It is very time and CPU consuming.

Therefore for many tasks, especially in the early stage of detector design, the fast simulation is used. It takes a less detailed description of the detector and does not simulate every particle step-by-step. Instead, it simulates the overall response of the (particular) detector in a parametric way.

Generated particles are transported inside the detector and they (their 4-momentum and/or position) are smeared taking into account the resolutions and efficiency. Those smeared particles may be analysed and treated as reconstructed particles, even though no hits were produced and no reconstruction was performed. All the detector effects (both physics processes that may encounter and detector resolutions) come from the smearing process (or rather the resolutions that were used for smearing).

The resolutions used in the smearing may come arbitrary from our knowledge of the detectors. In that case one applies a Gaussian smearing with a given standard deviation. That approach may be also used by the physicists to test how the detector resolution affect the results. That smearing is currently implemented in FCCSW.

More complex approach involves construction of the tables with the detector resolutions (pseudorapidity/momentum/particle dependent). They are calculated from a small (relatively) sample of full simulations of single-particle events. Single-particle events simplify the reconstruction process (they don't involve the pattern recognition etc.). Such resolutions are valid for tested detectors hence they may be used for smearing the particles with a better accuracy. Implementation of this approach is still in progress.

TODO: gflash

Both full and fast simulation can be performed in FCCSW using Geant4. Since the same tools are used for both of them, each simulation may be an interplay of both with full simulation performed in some volumes and fast simulation in others.


## 2. Example

To run fast simulation:

~~~{.sh}
./run gaudirun.py config/geant_fastsim.py
~~~

The differences between configuration file of the fast simulation (options/geant_fastsim.py) and the full simulation:

  * Geant configuration
    - physics list for the fast simulation is an overlay of the full sim physics list.
      It registers the parametrisation process to "fullphysics" list.

    ~~~{.py}
    from Configurables import G4FastSimPhysicsList
    physicslisttool = G4FastSimPhysicsList("Physics", fullphysics="G4FtfpBert")
    ~~~

    - action initialization contains the list of user actions to be set in Geant
      (in particular InitializeModelsRunAction essential for the fast sim models registration to the geometry).
      Also a simple smearing tool is created with a constant resolution.

    ~~~{.py}
    from Configurables import G4FastSimActions, G4ParticleSmearSimple
    smeartool = G4ParticleSmearSimple("Smear", sigma = 0.15)
    actionstool = G4FastSimActions("Actions", smearing=smeartool)
    ~~~

    - specification of tools that will construct the detector, physics list and user actions

    ~~~{.py}
    from Configurables import G4SimSvc
    geantservice = G4SimSvc("G4SimSvc",
                             detector='G4DD4hepDetector',
                             physicslist=physicslisttool,
                             actions=actionstool)
    ~~~


  * simulation
    - as for full sim, but using different tools to store the output of the simulation

    ~~~{.py}
    from Configurables import G4SimAlg,G4SaveSmearedParticles
    saveparticlestool = G4SaveSmearedParticles("G4SaveSmearedParticles")
    saveparticlestool.DataOutputs.particles.Path = "smearedParticles"
    saveparticlestool.DataOutputs.particlesMCparticles.Path = "particleMCparticleAssociation"
    geantsim = G4SimAlg("G4SimAlg",
                         outputs = ["G4SaveSmearedParticles/G4SaveSmearedParticles"])
    geantsim.DataInputs.genParticles.Path="allGenParticles"
    ~~~



3. Geant configuration: via GAUDI service G4SimSvc
----

(Sim/SimG4Components/src/G4SimSvc.h)

Main service for simulation with Geant4, G4SimSvc, contains the RunManager which derives from G4RunManager. Own implementation of G4RunManager was necessary in order to leave the event flow to GAUDI.
Any communication to the Geant Run Manager is handled by tthe service:
    - configuration:
      - geometry construction
      - physics list initialization
      - user action initialization
    - simulation:
      - passing an event to G4EventManager
      - retrieving a simulated event (with hits collections and other information)


### 3.1. Geometry construction

Geometry is provided by DD4hep via GAUDI's service GeoSvc. There is also a tool handling the detector construction using the GDML file though it is only for test purposes and it is not meant to be used in FCC simulation.

DD4hep translates the geometry into Geant4 geometry. That translation is done based on xml file as well as on the corresponding constructors, eg.

       DetectorDescription/Detectors/compact/ParametricSimEnvelopes.xml
       DetectorDescription/Detectors/src/SimpleTube_geo.cxx

Path to the XML file with detector description is a property of the service GeoSvc.
-----
  * In fast simulation user wants a specific behaviour in certain geometry volumes. That specific behaviour is described in classes derived from G4VFastSimulationModel. Geant will not perform normal transportation inside volumes with fast simulation models attached (providing that particle triggers that model).

  The fast simulation model needs to be attached to a G4Region object. That G4Region can contain one or many logical volumes (parts of the detector). Logical volumes are created by DD4hep translation.

  Currently a fast simulation model FastSimModelTracker (GeantSim/GeantFast/GeantFast/FastSimModelTracker.h) will be created and attached automatically for every G4LogicalVolume that contains 'tracker' in its name in DD4hep geometry constructor (with the intermediate step of creating the G4Region). Work is in progress to change that to be easier to modify (in particular to be set from xml file).

-----
  * In fast simulation user wants a specific behaviour in certain geometry volumes. That specific behaviour is described in classes derived from G4VFastSimulationModel. Geant will not perform normal transportation inside volumes with fast simulation models attached (providing that particle triggers that model).

  The fast simulation model needs to be attached to a G4Region object. That G4Region can contain one or many logical volumes (parts of the detector). Logical volumes are created by DD4hep translation. G4Region object is created automatically in InitializeModelsRunAction for any detector that has in its name 'Tracker', 'ECal', 'EMCal' or 'HCal'. Name of the detector is specified in DD4hep xml file under the tag:

    ~~~{.xml}
    <detector name ="CentralTracker">
    ~~~

  Currently a fast simulation model FastSimModelTracker will be created and attached automatically for every G4LogicalVolume that contains 'tracker'.
  There is also an ongoing work on implementation of the GFlash parametrisation that could be attached to the calorimeters.
  Attaching the models is done in InitializeModelsRunAction and can be controlled (switched on/off for tracker/ECal/HCal detectors) via flags set in a constructor (set to true by default).

TODO: What if I want another smearing


###
Geant4Simulation has two properties:

    simtype
    smearingtoolname

Property simtype describes which simulation type the user wants to perform: "fast" or "full". The default is "full".
Property smearingtoolname is valid only for fast simulation. It is a name of a class derived from ISmearingTool (GeantSim/GeantComponents/GeantComponents/ISmearingTool.h). It describes how the momentum/energy is smeared.



### 3.3. Physics List

Geant simulation requires the particles and processes definitions - the so-called physics list. The currently used physics list is FTFP_BERT, which is recommended by Geant4 for HEP.
-----
  * additionally for fast simulation:
     - "Coupled Transportation" is used which allows to invoke G4PathFinder which propagates the particle in the magnetic field.
     It is used within fast simulation model for tracker to compute the exit position of a particle from the volume (taking the momentum from the entrance to the tracker volume).
     - FastSimPhysics is registered as an additional process (GeantSim/GeantFast/GeantFast/FastSimPhysics.h). It attaches the fast simulation manager process to ALL the particles. That means that along with the standard processes such as transportation, multiple scattering etc. the particle can encounter 'fast simulation' process. That happens if a particle enters a volume with fast simulation model attached and if that particle fulfils trigger conditions (is charged in case of FastSimModelTracker).

     Details may be found at
                http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/ch05s02.html#sect.PhysProc.Param

-----


  * additionally for fast simulation:
    Fast simulation requires registration of another process, called in Geant FastSimManagerProcess. It can be attached to any existing physics list (like FTFP_BERT). The tool for the fast sim is called G4FastSimPhysicsList. It has a property 'fullphysics' that accepts a name of a tool with full sim physics definition.
Fast simulation manager process is attached to ALL the particles defined by the full sim physics list. That means that along with the standard processes such as transportation, multiple scattering etc. the particle can encounter 'fast simulation' process. That happens if a particle enters a volume with fast simulation model attached (which happens in InitializeModelsRunAction) and if that particle fulfils trigger conditions (is charged in case of FastSimModelTracker).

Furthermore, the fast simulation tool uses "Coupled Transportation"  which allows to invoke G4PathFinder which propagates the particle in the magnetic field. It is used within fast simulation model for tracker to compute the exit position of a particle from the volume (taking the momentum from the entrance to the tracker volume).

     Details may be found at
                http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/ch05s02.html#sect.PhysProc.Param



### 3.4. User Actions

Geant allows to specify actions that will be invoked at different stages of event processing.


-----
TrackingAction (GeantSim/GeantGeneral/GeantGeneral/TrackingAction.h) saves to the EDM the information on simulated particle (at the end of processing the particle with Geant).

To be able to compare generated particles (EDM MCParticle) with reconstructed particles (EDM Particle), one needs to create the association between those objects. That association is created while the event is translated from EDM to G4Event. A handle to MCParticle and Particle (the latter at that time empty) is stored in ParticleInformation (GeantSim/GeantGeneral/GeantGeneral/ParticleInformation.h). It is an object that Geant allows to be created for each G4PrimaryParticle. At the end of tracking, the information about the simulated particle (momentum, vertex position) is saved to the EDM Particle using the handle from ParticleInformation.

The user actions are the last (and optional) element of G4RunManager initialisation.
-----


For the fast simulation purposes, a InitializeModelsRunAction is created. It loops over all logical volumes in the world volume and creates G4Region for any volume with 'Tracker' in the detector name. All those regions become envelopes of the FastSimModelTracker.
Each model will use smearing tool for the momentum/energy smearing, so the name of the smearing tool needs to be passed to the InitializeModelsRunAction through the property of G4FastSimActions.

Smearing is performed using the GAUDI tool derived from ISmearingTool. User may specify how to smear a vector (CLHEP::ThreeVector) or a scalar (double). Currently there is only one implementation of a smearing tool called SimpleSmear (GeantSim/GeantComponents/src/SimpleSmear.h). The momentum/energy are smeared by multiplying Gaussian distribution with the mean $\mu=1$ and the standard deviation $\sigma$ as specified in the tool properties (default: 0.01).

Also, there is a tracking user action SaveParticlesUserAction that sets the information about the momentum, position and status of each particle at the end of tracking. The information ParticleInformation is associated with any G4PrimaryPaarticle and can be retrieved after the simulation from an event (see ... TODO).


## 4. Simulation in GAUDI algorithm G4SimAlg

Simulation algorithm handles all the communication between other algorithms and G4SimSvc.
It takes as input the EDM MCParticleCollection.
Also, a list of names to the output-saving tools can be specified as "outputs".


### 4.1. Event Processing

For each execution of the algorithm an EDM MCParticleCollection is translated into G4Event using the method G4Simulation::EDM2G4().

At the same time the ParticleCollection is created so that the association between MCParticle and Particle in EDM can be preserved. ParticleCollection will be the collection of the so-called 'reconstructed' particles. 'Reconstructed' in this case means the particles have undergone simulation process and the resulting (final) changes to the momentum are taken into account. At this stage it does not involve any reconstruction process.

### 4.2. Output

Saving the output from a simulated G4Event is performed by tools deriving from an interface IG4SaveTool.
Tools may have the outputs specified.
A method ::SaveOutput(const G4Event &aEvent) is meant to retrieve any useful information and save it to EDM.
Useful information means eg. hits collections (G4HCofThisEvent) or anything stored in an implementation of G4VUserEventInformation, G4VUserEventInformation, G4VUserTrackInformation, G4VUserPrimaryParticleInformation etc.

Existing tools store hits collections from tracker detectors (G4SaveTrackerHits) or calorimeters (G4SaveCalHits).

G4SaveTrackerHits stores trackClusters (EDM TrackClusterCollection), trackHits (EDM TrackHitCollection) and trackHitsClusters (EDM TrackClusterHitsAssociationCollection) for any hit collection with "Tracker" in its name.

G4SaveCalHits takes a property "caloType" that can be either "ECal", "EMCal" or "HCal. It stores caloClusters (EDM CaloClusterCollection) and caloHits (EDM CaloHitCollection) for any hit collection with name that contains the string defined in "caloType" property.

Hits collections' names are defined in DD4hep as readouts' names.
For additional information about DD4hep see TODO.


FAST SIM

G4SaveSmearedParticles tool stores particles (EDM ParticleCollection) and particlesMCparticles (EDM ParticleMCParticleAssociationCollection). At the end of tracking information about the particles momenta, status and vertex position is stored in ParticleInformation (SaveParticlesUserAction::PostUserTrackingAction()). Particle Information contains as well the reference to the EDM MCParticle. Therefore, for each particle in the event, based on ParticleInformation, a EDM output is created (including the association between MCParticle and Particle).
Those particles may be treated as 'reconstructed' particles. 'Reconstructed' in this case means the particles have undergone simulation process and the resulting (final) changes to the momentum are taken into account. At this stage it does not involve any reconstruction process.


5. Smearing
----

In case of running the simulation with simtype=='fast', user can define a smearing tool that should be used.

Smearing is performed using the GAUDI tool derived from ISmearingTool. User may specify how to smear a vector (CLHEP::ThreeVector) or a scalar (double). Currently there is only one implementation of a smearing tool called SimpleSmear (GeantSim/GeantComponents/src/SimpleSmear.h). The momentum/energy are smeared by multiplying Gaussian distribution with the mean $\mu=1$ and the standard deviation $\sigma$ as specified in the tool properties (default: 0.01).
