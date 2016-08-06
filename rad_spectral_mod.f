module rad_spectral
use rad_param
use rad_config
use rad_fields
use rad_gas_species, only : gasInfoType, rad_gas_species_init, getGasInfo
use rad_solid_species, only : solidParInfoType, rad_solid_species_init, getSolidInfo, smax
implicit none
private

! ---- Interfaces ------
public :: rad_spectral_init
public :: rad_spectral_calc
public :: rad_spectral_srad
public :: rad_spectral_final

interface rad_spectral_init
    module procedure init1
end interface

interface rad_spectral_final
    module procedure finalize1
end interface

interface rad_spectral_calc
    module procedure calculateSpectralFields
end interface

interface rad_spectral_srad
    module procedure calcRadSources
end interface
! ---- Data members ----
character(len=strLen) :: modelName_
integer :: modelId_

! ---- Parameters ------
! model ids
integer, parameter :: unknown=0, gray=1

contains

subroutine init1(config)
    use rad_spectral_gray, only : rad_spectral_gray_init
    implicit none
    type(configuration), intent(in) :: config
    modelName_ = config%SpectralModelName
    ! parse options
    select case (modelName_)
        case ("GRAY")
            modelId_ = gray
            call rad_spectral_gray_init(config)
        case default
            modelId_ = unknown
    end select
    call rad_gas_species_init(config)
    call rad_solid_species_init(config)
end subroutine

subroutine finalize1()
    use rad_spectral_gray, only : rad_spectral_gray_final
    implicit none
    select case (modelId_)
        case(gray)
            call rad_spectral_gray_final()
        case default
    end select
end subroutine
subroutine calculateSpectralFields()
    use functions, only : fluid_at
    implicit none
    type(gasInfoType) :: gasinfo
    type(solidParInfoType), dimension(smax) :: solidinfos
    real(dp),dimension(nq) :: kg ! gas absorption coefficient cm^-1
    real(dp),dimension(SMAX,nq) :: ks ! solid absorption coefficient cm^-1
    real(dp),dimension(SMAX,nq) :: scats ! solid scattering coefficient cm^-1
    real(dp),dimension(0:SMAX,nq) :: emiss !spectral emission 
    real(dp),dimension(1) :: T,kco2,kh2o
    
    integer :: ijk, m
    do ijk = 1, dimension_3
        if (fluid_at(ijk)) then
            call getGasInfo(ijk, gasinfo)
            do m=1, smax
                call getSolidInfo(ijk,m, solidinfos(m))
            end do
            call spectralCalc(gasinfo, solidinfos, kg, ks, scats, emiss)
            
            ! put cell calculations into fields
            k_g(ijk,:) = kg
            E(ijk,0,:) = emiss(0,:)
            do m=1, smax
                E(ijk,m,:) = emiss(m,:)
                k_s(ijk,m,:) = ks(m,:)
                scat(ijk,m,:) = scats(m,:)
            end do
        endif    
    enddo
    
    call extrapolation()

    ! unit conversions
    if (UNITS=='CGS') then
        ! gas and solid absorption coeffs are in cm^-1 already so no change
        scat = scat/1.d2
    else if (UNITS=='SI') then
        k_g = k_g*1.d2
        k_s = k_s*1.d2
    else 
        print*, "Radiative property calculation: unknown units: ", UNITS
        print*, "Known units are CGS and SI"
    endif
end subroutine

subroutine calcRadSources()
    implicit none
    integer iq, m
    Srad = 0.d0
    do iq = 1, nq
        Srad(:,0) = Srad(:,0)+k_g(:,iq)*(G(:,iq)-4.d0*pi*E(:,0,iq)) 
        do m=1, smax
            Srad(:,m) = Srad(:,m)+k_s(:,m,iq)*(G(:,iq)-4.d0*pi*E(:,m,iq))
        end do
    end do
    ! by far Srad units are W/(m^3) in SI or erg/(cm^3.s) in CGS
    ! C_pg uses J/kg.K in SI or cal/g.K in CGS
    ! 1 erg = 2.39005736d-8 cal
    if (UNITS=='CGS') Srad = Srad * 2.39005736d-8
end subroutine

subroutine spectralCalc(gasinfo, solidinfos, kg,ks, scats, emiss)
    use rad_spectral_gray, only : rad_spectral_gray_calc
    implicit none
    type(gasInfoType), intent(in) :: gasinfo
    type(solidParInfoType),dimension(smax), intent(in) :: solidinfos
    real(dp),dimension(nq), intent(out) :: kg 
        ! gas absorption coefficient cm^-1
    real(dp),dimension(SMAX,nq), intent(out) :: ks 
        ! solid absorption coefficient cm^-1
    real(dp),dimension(SMAX,nq), intent(out) :: scats 
        ! solid scattering coefficient cm^-1
    real(dp),dimension(0:SMAX,nq), intent(out) :: emiss 
        !spectral emission 
    
    select case(modelId_)
        case (gray)
            call rad_spectral_gray_calc(gasinfo, solidinfos, kg, ks, scats, emiss)
        case default
    end select 
end subroutine

subroutine extrapolation()
    use rad_util, only : extrapolateGhostCells
    implicit none
    integer :: iq, is
    do iq=1, nq
        call extrapolateGhostCells(k_g(:,iq))
        do is = 1, size(k_s, 3)
            call extrapolateGhostCells(k_s(:,is,iq))
            call extrapolateGhostCells(scat(:,is,iq))
        end do
    end do
            
end subroutine
end module rad_spectral
