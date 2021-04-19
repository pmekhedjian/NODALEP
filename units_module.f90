!
!.....This is an overly simplified, light version of units_module.f!
!.....Consult original Agile-IDSA's units_module.f for a full
!.....description of all these values....

      module units_module

      type units_in_cgs
        real :: alpha,pi,sinsqthetw,MeV          !numbers
        real :: c,mcsq_e,mb,NCSCS,G        !cgs constants
        real :: h,Q,n_mass,p_mass,G_F_IDSA,G_F_Burrows,G_F_Bruenn  !non-cgs constants
      end type units_in_cgs

      type(units_in_cgs), parameter :: units = units_in_cgs(-1.23,&
      3.1415926535898d+00,0.2325,1.602192d-06,2.997924562d+10,0.511,1.674d-24,&
      1.705d-44,6.67384e-08,4.1356943d-21,1.29332,939.565,&
      938.272,8.957d-44,8.8179d-44,5.18d-44)

!     -1.23                  ! alpha (axial-vector coupling constant)
!     3.1415926535898d+00,   ! pi (non-dim, obvs.)
!     0.2325,                ! sinsqthetw
!     1.602192d-06           ! convert erg to MeV (1MeV=1.602192e-06erg)
!     2.997924562d+10,       ! c (in cm/s)
!     0.511,                 ! mcsq_e (in MeV)
!     1.674d-24,             ! mb (in grams)
!     1.705d-44,             ! NCSCS(sigma_o) [cm^2] (neutrino cross sec.)
!     6.67384e-08            ! G [cm^3/s^2/g]
!     4.1356943d-21,         ! h [MeV*s] (Planck's, obvs. not cgs)
!     1.29332                ! Q=neutron-proton mass [MeV]
!     939.565                ! n_mass=neutron mass [MeV]
!     938.272                ! p_mass=proton mass [MeV]
!     8.957d-44              ! G_F_IDSA [MeV*cm^3]
!     8.8179d-44             ! G_F_Burrows (2006) [MeV*cm^3]
!     5.18d-44               ! G_F_Bruenn (1985) [cm^2*MeV^-2]

      end module units_module
