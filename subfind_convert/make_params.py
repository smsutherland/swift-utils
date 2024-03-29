def make_param_file(params):
    return f"""%% parameters for cosmo_box_star_formation_3d

%----  Relevant files
InitCondFile                            /cosma7/data/dp004/dc-love2/data/CAMELS/ICs/ics_camels_cv_0.hdf5
OutputDir                               ./
SnapshotFileBase                        subsnap
OutputListFilename                      /cosma7/data/dp004/dc-love2/codes/CAMELS-SWIFT/arepo_subfind/output_list_arepo_subfind.txt

%---- File formats
ICFormat                                3
SnapFormat                              3

%---- CPU-time limits
TimeLimitCPU                            90000
CpuTimeBetRestartFile                   12000
% FlushCpuTimeDiff                        120

ResubmitOn                              0
ResubmitCommand                         my-scriptfile

%----- Memory alloction
MaxMemSize                              16000

%---- Characteristics of run
TimeBegin                               0.0078125  % Begin of the simulation; z=127
TimeMax                                 1.0  % End of the simulation z=0

%---- Basic code options that set the type of simulation
ComovingIntegrationOn                   1
PeriodicBoundariesOn                    1
CoolingOn                               0
StarformationOn                         0

%---- Cosmological parameters (Planck cosmology)
Omega0                                  {params["Omega_m"]}
OmegaLambda                             {1 - params["Omega_m"]}
OmegaBaryon                             {params["Omega_b"]}
HubbleParam                             {params["H0"]}
BoxSize                                 {params["BoxSize"]}

%---- Output frequency and output parameters
OutputListOn                            1
TimeBetSnapshot                         0.0
TimeOfFirstSnapshot                     0.0
TimeBetStatistics                       0.01
NumFilesPerSnapshot                     1
NumFilesWrittenInParallel               1

%---- Accuracy of time integration
TypeOfTimestepCriterion                 0
ErrTolIntAccuracy                       0.012
CourantFac                              0.3
MaxSizeTimestep                         0.005
MinSizeTimestep                         1.0e-9

%---- Treatment of empty space and temperature limits
InitGasTemp                             244.8095
MinGasTemp                              5.0
MinimumDensityOnStartUp                 1.0e-20
LimitUBelowThisDensity                  0.0
LimitUBelowCertainDensityToThisValue    0.0
MinEgySpec                              0.0

%---- Tree algorithm, force accuracy, domain update frequency
TypeOfOpeningCriterion                  1
ErrTolTheta                             0.7
ErrTolForceAcc                          0.0025
MultipleDomains                         8
TopNodeFactor                           2.5
ActivePartFracForNewDomainDecomp        0.01
 
%---- Initial density estimate
DesNumNgb                               64
MaxNumNgbDeviation                      4

%---- System of units
UnitLength_in_cm                        3.085678e21    %  1.0 kpc
UnitMass_in_g                           1.989e43       %  1.0e10 solar masses
UnitVelocity_in_cm_per_s                1e5           %  1 km/sec
GravityConstantInternal                 0

%---- Gravitational softening lengths

SofteningComovingType0    {params["SofteningComoving"]}
SofteningComovingType1    {params["SofteningComovingType1"]}
SofteningComovingType2    {params["SofteningComoving"]}
SofteningComovingType3    {params["SofteningComoving"]}
SofteningComovingType4    {params["SofteningComoving"]}
SofteningComovingType5    {params["SofteningComoving"]}

SofteningMaxPhysType0     {params["SofteningMaxPhys"]}
SofteningMaxPhysType1     {params["SofteningMaxPhysType1"]}
SofteningMaxPhysType2     {params["SofteningMaxPhys"]}
SofteningMaxPhysType3     {params["SofteningMaxPhys"]}
SofteningMaxPhysType4     {params["SofteningMaxPhys"]}
SofteningMaxPhysType5     {params["SofteningMaxPhys"]}

SofteningTypeOfPartType0        1
SofteningTypeOfPartType1        1
SofteningTypeOfPartType2        1
SofteningTypeOfPartType3        1
SofteningTypeOfPartType4        1
SofteningTypeOfPartType5        2

GasSoftFactor             2.5
MinimumComovingHydroSoftening   0.25
AdaptiveHydroSofteningSpacing   1.2


%----- Mesh regularization options
CellShapingSpeed                        0.5
CellMaxAngleFactor                      2.25
ReferenceGasPartMass                    0
TargetGasMassFactor                     1
RefinementCriterion                     1
DerefinementCriterion                   1

%----- Subfind
ErrTolThetaSubfind                      0.7
DesLinkNgb                              20

%---- Parameters for star formation model
% CritPhysDensity                         0       % critical physical density for star formation (in cm^(-3))
% MaxSfrTimescale                         2.27    % in internal time units
% CritOverDensity                         57.7    % overdensity threshold value
% TempSupernova                           5.73e7  % in Kelvin
% TempClouds                              1000.0  % in Kelvin
% FactorEVP                               573.0
% TemperatureThresh                       1e+06
% FactorSN                                0.1

% TreecoolFile                            ./TREECOOL_ep
"""
