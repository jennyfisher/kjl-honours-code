!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: methanol_mod
!
! !DESCRIPTION: Module METHANOL\_MOD contains subroutines for methanol sources and
!  sinks in the full chemistry simulation.
!\\
!\\
! !INTERFACE:
!
      MODULE METHANOL_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: OCEAN_SOURCE_METHANOL
      PUBLIC :: OCEAN_SINK_METHANOL
      PUBLIC :: CLOUD_SINK_METHANOL
!
! !REMARKS:
!  10 Jun 2013 - J.A. Fisher - copied from offline methanol simulation
!EOP
!------------------------------------------------------------------------------
!BOC
! !DEFINED PARAMETERS:
!
      ! Avogadro's number
      REAL*8,  PARAMETER   :: AVO        = 6.022d23

      ! MOH moelcular weight
      REAL*8,  PARAMETER   :: MW_MOH     = 32d-3

      ! Molecules MOH / kg MOH 
      REAL*8,  PARAMETER   :: XNUMOL_MOH = AVO / MW_MOH

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ocean_source_methanol
!
! !DESCRIPTION: Subroutine OCEAN\_SOURCE\_METHANOL specifies the ocean source of 
!  methanol.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OCEAN_SOURCE_METHANOL( I, J, METHANOL, State_Met )
!
! !USES:
!
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE ERROR_MOD, ONLY : CHECK_VALUE
      USE DIAG_MOD,  ONLY : AD64
      USE GRID_MOD,  ONLY : GET_AREA_CM2
      USE TIME_MOD,  ONLY : GET_TS_EMIS

      USE CMN_SIZE_MOD                  ! Size parameters
      USE CMN_DIAG_MOD                  ! ND64
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: I         ! GEOS-Chem longitude index
      INTEGER, INTENT(IN)    :: J         ! GEOS-Chem latitude index
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(INOUT) :: METHANOL   ! Methanol emission [atoms C/s]
!
! !REMARKS:
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8                 :: HENCONST, KKL, Cl
      REAL*8                 :: DTSRCE, OCEAN_MOH, AREA_CM2, FOCEAN

      !=================================================================
      ! OCEAN_SOURCE_METHANOL begins here!
      !=================================================================

      ! Emission timestep in seconds
      DTSRCE   = GET_TS_EMIS() * 60d0

      ! Fraction of (I,J) that is ocean
      FOCEAN   = 1d0 - State_Met%FRCLND(I,J)

      ! Area of grid box (I,J) in cm^2
      AREA_CM2 = GET_AREA_CM2( I, J, 1 )

      !===========================================================
      ! Only compute ocean source if there is more than 50% ocean
      ! in the grid box, and if it is not ice (albedo > 0.4)
      ! (mje, rvm, bmy, 11/26/01)
      !===========================================================
      IF ( FOCEAN > 0.5d0 .AND. State_Met%ALBD(I,J) <= 0.4d0 ) THEN 

         ! Calculate Henry constant and transfer velocity
         CALL GET_MOH_OCN_PARAMS( I, J, HENCONST, KKL, State_Met )

!***************************************************************************
! dbm
!
! MOH air-sea exchange
! -------------------------------
! The flux of gas (MOH) through the air-sea interface can be expressed as:
!
!  Liss & Slater, Nature 247, 181-184, 1974:
!  F = KKl * ( Cg / H - Cl)
!
! where
! (a)KKl    :  total transfer velocity expressed on a liquid phase (seawater) 
!               basis, 
!
!                 1      1        1
!               ----- = ---  + --------
!                KKl     kl     H * kg
!
!               where 
!                (a.1) kl : transfer velocity in liquid phase (seawater) (cm/hr)
!
!                  Nightingale et al., Global Biogeochem. Cyc., 14, 373-387 (2000):
!                   
!                  kl/kref = ( Sc / Scref )^(-1/2)
!                     Parameterized relationship for available observations:
!                     k600 = 0.222*U^2 + 0.333*U 
!                           where k600 = kl @ Sc = 600
!                                 U = wind speed @ 10m (m/s)
!
!                  So that k1 = k600 * (Sc / 600 )^(-1/2)
!
!                  Sc is the Schmidt number for MOH:
!
!                    Sc = A1 * T^4 + A2 * T^3 + A3 * T^2 + A4 * T + A5
!
!                (a.2) kg : transfer velocity in gas phase
!
!                  Asher (1997):                    
!                  kg = (5.2E-5 + 3.2E-3 * u) * SQRT(WTH2O / WTMOH)
!                  
!                  where WTH2O =18.0 and WTMOH =32.0 are the molecular weight 
!                  of H2O and MOH in (g/mole), respectively. kg is in m/s.
! 
! (b) Cg, Cl :  gas (MOH) concentration in gas and liquid phase, respectively 
! (c) H      :  dimensionless Henry's law constant 
!                 [mass MOH/volume air]/[mass MOH/volume H2O]
!
!                        1
!               H = ------------- 
!                    KH * Rg * T
!
!               where
!                (c.1) KH  :  Henry's law constant 
!                              [X]aq / P(X) (moles / litr / atm) 
!                (c.2) Rg  :  Universal gas constant
!                (c.3) T   :  Temperature
!
! Temperature dependence of KH (Seinfeld & Pandis):
!
!               KH(T2) = KH(T1) * EXP ( (DELTAH / R) * (1/T1 - 1/T2))
!
! At equilibrium, Cl = Cg / H
!               Introduce S = saturation ratio (0-1)
!               Cl = S*Cg/H
!               At 0, Cl = 0.  At S = 1, Cl = Cg/H (equilibrium).
!
! And the flux becomes
!
!  F = KKl * ( Cg / H - Cl) = KKl * ( Cg / H - S * Cg / H )
!                           = KKl * ( 1 - S ) * Cg/H
!
! (d)Reference:
!
!             Liss, P.S., and P.G. Slater, Flux of gases across the air-
!                sea interface, Nature, 247, 181-184, 1974
!             Liss, P.S., and L. Merlivat, Air-sea gas exchange rates: 
!                Introduction and synthesis, in The Role of Air-Sea Exchange 
!                in Geochemical Cycling, pp113-127, edited by P. Buat-Menard, 
!                D.Reidel, Norwell, MA, 1986
!             Edwards, T.J., G. Maurer, J. Newman, and J.M. Prausnitz, Vapor-
!                Liquid Equilibria in Multicomponent Aqueous Solutions of
!                Volatile Weak Electrolytes, AIChE J., 24, 966-976, 1978    
!             Asher, W., The sea-surface microlayer and its effect on global
!                air-sea gas transfer, in The sea surface and global change,
!                edited by P.S. Liss and R.A. Duce, Cambridge Univ. Press, 
!                Cambridge, UK, 1997  
!             Wilke, C.R., and P. Chang, Correlation of diffusion coefficients
!                in dilute solutions, AIChE J., 1, 264-270, 1955
!             Jahne, B., G. Heinz, and W. Dietrich, Measurement of the 
!                diffusion coefficients of sparingly soluable gases in 
!                water with a modified Barrer method, J. Geophys. Res., 
!                92, 10,767-10,776, 1987                      
!             Wanninkhof,R., Relationship Between Wind Speed and Gas Exchange 
!                Over the Ocean, J. Geophys. Res., 97, 7373-7382, 1992
!
! -- qli, 1/24/1999, 4/16/99, dbm 05/20/2005
!*****************************************************************************

         !-------------------------------------------------------
         ! Compute gross ocean source
         ! sea to air FLUX(kg/cm2/s) = KKl(cm/s) * Cl(kg/cm3)
         !-------------------------------------------------------

            ! Assume a constant seawater concentration of 118 nM
            !   (from Williams et al. (2004), GRL, VOL. 31, L23S06, 
            !   doi:10.1029/2004GL020012) and convert to kg/cm3
            Cl = 118.0d-9 * MW_MOH / 1000.0d0

            ! Correct for the fraction of the grid cell that is ocean
            ! and compute flux (kg/cm2/s)
            OCEAN_MOH = KKL * FOCEAN * Cl

            ! Convert to kg / box / step
            OCEAN_MOH = OCEAN_MOH * DTSRCE * AREA_CM2

      ELSE

         ! If there is less than 50% water in the grid box, zero 
         ! the ocean source from methanol (mje, rvm, bmy, 11/26/01)
         OCEAN_MOH = 0d0

      ENDIF

      ! Add ocean source to total biogenic source in [molec/box/s]
      METHANOL = METHANOL + ( OCEAN_MOH * XNUMOL_MOH / DTSRCE )

      !=================================================================
      ! ND64 diag -- save METHANOL from the ocean in [molec/cm2/s]
      !=================================================================
      ! kg/box/timestep * (molec/kg) * (1/s) * (1/cm2)
      IF ( ND64 > 0 ) 
     &     AD64(I,J,1,4) = AD64(I,J,1,4) + ( OCEAN_MOH * XNUMOL_MOH ) /
     &                                     ( AREA_CM2  * DTSRCE   )

      END SUBROUTINE OCEAN_SOURCE_METHANOL
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ocean_sink_methanol
!
! !DESCRIPTION: Subroutine OCEAN\_SINK\_METHANOL applies the ocean sink to global
!  methanol concentrations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OCEAN_SINK_METHANOL( METHANOL, State_Met )
!
! !USES:
!
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE DIAG_MOD, ONLY : AD64
      USE GRID_MOD, ONLY : GET_AREA_CM2
      USE TIME_MOD, ONLY : GET_TS_CHEM

      USE CMN_SIZE_MOD                              ! Size parameters 
      USE CMN_DIAG_MOD                              ! ND64
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL*8, INTENT(INOUT) :: METHANOL(IIPAR,JJPAR)   ! methanol mass [kgC]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER                :: I,        J
      REAL*8                 :: HENCONST, KKL,       RKT,      DV
      REAL*8                 :: DTCHEM,   OCEAN_MOH, AREA_CM2, FOCEAN

      !=================================================================
      ! OCEAN_SINK_METHANOL begins here! 
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Loop over grid boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box area in cm2
         AREA_CM2 = GET_AREA_CM2( I, J, 1 )

         ! Fraction of grid box that is ocean
         FOCEAN = 1d0 - State_Met%FRCLND(I,J)

         !===========================================================
         ! Only compute ocean sink if there is more than 50% ocean
         ! in the grid box, and if it is not ice (albedo > 0.4)
         ! (mje, rvm, bmy, 11/26/01)
         !===========================================================
         IF ( FOCEAN > 0.5d0 .and. State_met%ALBD(I,J) <= 0.4d0 ) THEN

            ! Calculate Henry constant and transfer velocity
            CALL GET_MOH_OCN_PARAMS( I, J, HENCONST, KKL, State_Met )

            ! Net air-to-sea deposition velocity
            ! cm/s * unitless/ unitless = cm/s
            DV = KKL * FOCEAN / HENCONST

            ! Pre-compute the exponential term [1/timestep]
            ! Just use surface box for simplicity
            ! Convert BXHEIGHT from m to cm
            ! (mixing done elsewhere?)
            RKT = DV * DTCHEM / (100d0 * State_Met%BXHEIGHT(I,J,1) )

            ! Ocean sink in [kg]
            OCEAN_MOH = METHANOL(I,J) * ( 1 - EXP(-RKT) )

            ! Error trap (loss can't be greater than total)
            IF ( OCEAN_MOH > METHANOL(I,J) ) OCEAN_MOH = METHANOL(I,J)

         ELSE

            ! If there is less than 50% water in the grid box, or  
            ! if there is ice on the ocean, then zero the ocean sink
            ! for methanol (mje, rvm, bmy, 11/26/01)
            OCEAN_MOH = 0d0

         ENDIF

         ! Apply exponential loss
         METHANOL(I,J) = METHANOL(I,J) - OCEAN_MOH

         !===========================================================
         ! Diagnostics: save METHANOL lost to ocean in [molec/cm2/s]
         !===========================================================
            IF ( ND64 > 0 ) AD64(I,J,1,5) = 
     &           AD64(I,J,1, 5) + OCEAN_MOH * XNUMOL_MOH /
     &                           ( AREA_CM2 * DTCHEM )
      ENDDO
      ENDDO

      END SUBROUTINE OCEAN_SINK_METHANOL
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cloud_sink_methanol
!
! !DESCRIPTION: Subroutine CLOUD\_SINK\_METHANOL applies the cloud uptake sink
!   to methanol concentrations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLOUD_SINK_METHANOL( METHANOL, State_Met )
!
! !USES:
!
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE COMODE_MOD,    ONLY : CSPEC, JLOP
      USE DIAG_MOD,      ONLY : AD64
      USE GRID_MOD,      ONLY : GET_AREA_CM2
      USE TIME_MOD,      ONLY : GET_TS_CHEM
      USE TRACERID_MOD,  ONLY : IDOH

      USE CMN_SIZE_MOD                              ! Size parameters 
      USE CMN_DIAG_MOD                              ! ND64, LD64
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL*8, INTENT(INOUT) :: METHANOL(IIPAR,JJPAR,LLPAR)   ! methanol mass [kgC]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER                :: I,      J,         L,       JLOOP
      REAL*8                 :: LIQUID_FACTOR,     L_WATER
      REAL*8                 :: K_METH, K_AQ,  K,  CLDF1
      REAL*8                 :: DTCHEM, CLOUD_MOH, AREA_CM2
      REAL*8                 :: OH_AQ,  OH_MOLEC_CM3

      !=================================================================
      ! CLOUD_SINK_METHANOL begins here! 
      !=================================================================
!------------------------------------------------------------------------------
! dbm
!
! First order loss rate (s^-1) for atmospheric MOH:
!
! k = H * L * R * T * k * [OH(aq)]
!
!              Where
!                   H = Henry's Law constant (M / atm)
!                   R = gas constant
!                   T = Temperature in K
!                   k = Aqueous phase reaction rate
!                   [OH(aq)] = OH concentration in cloud drops (moles/L)
!                   L = dimensionless liquid water content
!                           Where                             
!                                L = 4*OPTD*r/(3*Q*BXHEIGHT)
!                                OPTD = cloud optical depth
!                                r = cloud effective drop radius (10
!                                microns)
!                                Q = cloud droplet extinction efficiency
!                                (2)
!                                BXHEIGHT = height of grid cell (m)
!                                --> multiply by 1E-6 to get length
!                                units to cancel
!------------------------------------------------------------------------------

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Loop over grid boxes
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box area in cm2
         AREA_CM2 = GET_AREA_CM2( I, J, 1 )

         ! JLOOP = SMVGEAR 1-D grid box index
         JLOOP = JLOP(I,J,L)

         ! Take OH from the SMVGEAR array CSPEC
         ! OH is defined only in the troposphere
         IF ( JLOOP > 0 ) THEN 
            OH_MOLEC_CM3 = CSPEC(JLOOP,IDOH)
         ELSE 
            OH_MOLEC_CM3 = 0d0
         ENDIF

         ! Assume all liquid above 268K
         IF ( State_Met%T(I,J,L) > 268d0 ) THEN 
            LIQUID_FACTOR = 1D0
         ! Assume linear transition from 0-1 between 248 & 268K
         ELSE IF ( State_Met%T(I,J,L) >= 248d0 .AND.
     &             State_Met%T(I,J,L) <= 268d0 ) THEN 
            LIQUID_FACTOR = (State_Met%T(I,J,L)-248D0)/20D0
         ! Assume no liquid below 248K
         ELSE
            LIQUID_FACTOR = 0D0
         ENDIF

         ! Henry's Law constant (M/atm)
!           KH(T2) = KH(T1) * EXP ( (DELTAH / R) * (1/T1 - 1/T2))
!                      KH(298K) = 220 M atm^-1 (Seinfeld and Pandis)
!                      DELTAH = 9.7 kcal/mole; DELTAH/R = -4881 K
         K_meth = 220d0 * EXP( 4900d0 * ( 298d0 - State_Met%T(I,J,L) ) /
     &          ( 298d0 * State_Met%T(I,J,L) ) )

!------------------------------------------------------------------------------
! dbm:
! Aqueous reaction rate MOH + OH --> products
! Rate constant is from Elliot and McCracken, Radiat. Phys. Chem. V.33,
! 69-74, 1989
!
! Ea = 4.8 kJ/mole
!      R = 8.31 J/mole/K --> E/R = 578 K
!
! K = A*EXP(-E/(RT))
!      K @ 20C = 10.1 E 8 dm^3/mol/s (i.e. M^-1 s^-1)
!      So that
!      A = 10.1 E 8 / EXP(-E/R/293)
!
! K @ T(I,J,L):
!      K = { 10.1 E 8 / EXP(-E/R/293) } * EXP(-E/(RT))
!        = 10.1 E 8 * EXP{ -578 * (1/T - 1/293) }
!------------------------------------------------------------------------------
         K_AQ  = 1.01d9 * EXP( -578d0 * ( 293d0 - State_Met%T(I,J,L) ) /
     &         ( 293d0 * State_Met%T(I,J,L) ) )

         !kcw: need to multiply OPTD by cloud frac^1.5 to 
         !account for the fact that OPTD is in-cloud OPTD
         !1.5 factor accounts for random overlap of clouds
         CLDF1 = State_Met%CLDF(L,I,J)
         IF ( CLDF1 < 0d0 ) CLDF1 = 0d0 !remove negative values
         CLDF1 = CLDF1**1.5

         ! Calculate Liquid water content
         ! kcw: multiply by cldf here
         L_WATER = 4d0 * State_Met%OPTD(L,I,J) * CLDF1 * 10d0 * 1D-6 /
     &           ( 3D0 * 2D0 * State_Met%BXHEIGHT(I,J,L) )

         ! Assume aqueous phase OH (moles/L) = 1E-19*atmosheric OH
         ! (in molec/cm3)
         OH_AQ = 1D-19 * OH_MOLEC_CM3

         ! Compute first order loss rate
         ! H (M/atm) * L(unitless) * R(L*atm/K/mole) * T (K)
         !  * k(M^-1 s^-1) * OH_AQ (M) = s^-1
         K = K_METH * L_WATER * 8.21D-2 * State_Met%T(I,J,L) *
     &       K_AQ   * OH_AQ   * LIQUID_FACTOR

         ! Error trap - force positive
         IF ( K < 0d0 ) K = 0d0

         ! Calculate sink [kg]
         CLOUD_MOH = METHANOL(I,J,L) * ( 1 - EXP( -K * DTCHEM ) )

         ! Error trap (loss can't be greater than total)
         IF ( CLOUD_MOH > METHANOL(I,J,L) ) CLOUD_MOH = METHANOL(I,J,L)

         ! Apply exponential loss
         METHANOL(I,J,L) = METHANOL(I,J,L) - CLOUD_MOH

         !===========================================================
         ! Diagnostics: save METHANOL cloud sink in [molec/cm3/s]
         !===========================================================
            IF ( ND64 > 0 .and. L <= LD64 ) AD64(I,J,L,6) = 
     &           AD64(I,J,L,6) + CLOUD_MOH * XNUMOL_MOH / (DTCHEM *
     &           AREA_CM2 * State_Met%BXHEIGHT(I,J,L) * 100d0 )
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE CLOUD_SINK_METHANOL
!EOC
!
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_moh_ocn_params
!
! !DESCRIPTION: Subroutine GET\_MOH\_OCN\_PARAMS calculates the Henry constant
!  and air-to-sea transfer velocity (dbm, 06/13/06)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_MOH_OCN_PARAMS( I, J, HENCONST, KKL, State_Met )
! !USES:
!
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: I         ! GEOS-Chem longitude index
      INTEGER, INTENT(IN)    :: J         ! GEOS-Chem latitude index
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT)   :: HENCONST ! Henry's Law constant
      REAL*8,  INTENT(OUT)   :: KKL      ! Air-to-sea transfer velocity
!
! !REMARKS:
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      ! Coefficients for fitting the Schmidt number for methanol
      ! Derived using polynomial fit (code provided by qli)
      REAL*8,  PARAMETER     :: A0 =  2292.677d0
      REAL*8,  PARAMETER     :: A1 = -94.9917d0
      REAL*8,  PARAMETER     :: A2 =  1.538656d0
      REAL*8,  PARAMETER     :: A3 = -0.009837149d0

      REAL*8                 :: KL, KL600, KG
      REAL*8                 :: SC, U, TC
      !kcw (6/2012) Add these parameters for Johnson [2010] update
      REAL*8                 :: USTAR, CD

      ! External functions
      REAL*8, EXTERNAL       :: SFCWINDSQR

      !=================================================================
      !  GET_MOH_OCN_PARAMS begins here!
      !=================================================================

      !-------------------------------------------------------
      ! Compute H
      !-------------------------------------------------------

      ! For MOH use H(298K) = 220(M atm-1) from Seinfeld + Pandis
      ! 220 M atm-1 * Rg * 298 K = 5382
      !
      ! And use -9.7(kcal mol-1) for heat of dissolution coef
      ! R = 1.9872 E-3 kcal / K / mole
      ! -9.7 / R = -4881
      ! Compute 1/H at local temperature 
      HENCONST = (5382d0)*
     &     EXP(-4881d0*(1./298d0 - 1./State_Met%TS(I,J)))

      ! Want Henry constant exressed as
      ! concentration in air/concentration in water.
      ! Take reciprocal.  
      HENCONST = 1d0/HENCONST
      ! Now HENCONST = dimensionless H 
      ! [mass MOH/volume air]/[mass MOH/volume H2O]

      !-------------------------------------------------------
      ! Compute schmidt number
      !-------------------------------------------------------
      ! Magnitude of surface wind [m/s]
      ! SFCWINDSQR(I,J) is needed since this will compute the 
      ! square of the surface wind correctly for all GEOS models 
      U     = SQRT( SFCWINDSQR( I, J ) )

      ! TC is temperature in Celsius
      TC    = State_Met%TS(I,J) - 273.15d0
      ! Bug fix: make sure -5 <= TC <= 30 (mje, rvm, bmy, 11/26/01) 
      TC    = MIN( MAX( TC, -5d0 ), 30d0 )

      ! SC is Schmidt # for MOH [unitless]
      SC    = A0 + A1*TC + A2*TC**2 + A3*TC**3

      ! KL is conductance for mass transfer in liquid phase 
      ! (Wanninkhof 1992), which has units of [cm/hr]
      !KL600 = ( 0.222d0 * U * U ) + ( 0.333D0 * U )
      !KL    = KL600 * ( SC / 600d0 )**( -.5d0 )
      !Update (kcw, 6/2012) Nightingale [2010]
      KL = ( 0.24d0*U*U + 0.061d0*U ) * SQRT( 600d0/Sc )

      ! Updated KG to the Johnson [2010] parameterization (kcw, 6/2012)
      ! USTAR = friction velocity (U* in Johnson [2010]
      USTAR  = SQRT(6.1d-4 + U*6.3d-5)*U
      ! CD = drag coefficient
      CD = (USTAR/U)**2

      ! KG is airside transfer velocity (Johnson 2010)
      ! Multiply KG by 360000 to convert from [m/s] to [cm/hr]
      KG = 1d-3 + (USTAR/(13.3*SC**(1/2) + CD**(-1/2)-5+LOG(SC)
     &         /0.8))
      KG = KG * 360000d0

      !-------------------------------------------------------
      ! Compute KKl (total air to sea transfer velocity expressed on a 
      ! liquid phase basis)
      !-------------------------------------------------------

      ! KKL is the air-to-sea transfer velocity (Liss and Slater 1974)
      ! Divide KKL by 3600 to convert from [cm/hr] to [cm/s]
      KKL   = 1d0 / ( 1d0/KL + 1d0/( HENCONST * KG ) )
      KKL   = KKL / 3600d0

      ! Return to calling program
      END SUBROUTINE GET_MOH_OCN_PARAMS
!EOC
      END MODULE METHANOL_MOD
