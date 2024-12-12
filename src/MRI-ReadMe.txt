#================================================================================

                This document describes how the MRI is implemented
                into GENEC and how one can use it to run models
#================================================================================

2 MRI implementations: With circulation by meridional currents as Advective and the other with Diffusive.

Advective Implementation: (Must set imagn==0)
    -Set parameter mri==1.
    -The chemical gradient term in minimum condition can be weighted by fmu, default set to 1.
    - In diffusion.f90 one must turn on Advective section at line 580.
    -MRI diffusion coefficient added in diffusion.f90 subroutine with D_shear coefficient.

Diffusive Implementation: (Must set imagn==1)
    -For just MRI set parameter mri==1.
    -For just TS set parameter mri==2.
    -For both set parameter mri==3.
    -Warning fmu parameter only impacts MRI minimum condition.
    -Warning TS chemical mixing not perfected, best to switch it off for first tests.
    -MRI and TS diffusion coefficents computed in magmod.f90 alongside the relevant magnetic parameters

There are changes in inputparam with the adding of variables fmu and mri. These changes also are made in the inputparam for the
inifile generation.

#For any questions contact Adam Griffiths/ adam.griffiths@uv.es
