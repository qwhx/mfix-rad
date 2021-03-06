module rad_gas_species
use rad_param
use rad_config
use rad_fields
use rxns, only : species_g, dim_m, dim_n_all, dim_n_g
use fldvar, only : ep_g, T_g, P_g, X_g
implicit none
private

! ---- Data types -----
type gasInfoType
    real(dp) :: p ! pressure unit bar
    real(dp) :: T ! temperature unit Kelvin
    real(dp), dimension(nSp) :: C ! mole fractions
!    real(dp) :: sootfv ! soot volume fractions unit ppm
end type gasInfoType

! ---- Interfaces ------
public :: gasInfoType
public :: rad_gas_species_init
public :: getGasInfo

interface rad_gas_species_init
    module procedure init1
end interface

! ---- Data members ----
! -- species data --
integer, dimension(DIM_N_G) :: radSpeciesId

logical, dimension(DIM_N_G) :: isRadSpecies


contains

subroutine init1(config)
    implicit none
    type(configuration), intent(in) :: config
    integer :: i
    isRadSpecies =.false.
    radSpeciesId = 0
    do i =1, dim_n_g
        select case (trim(species_g(i)))
        case ("CO2")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = CO2
        case ("H2O")
            isRadSpecies(i) = .true.
            radSpeciesId(i) = H2O
        case default
            isRadSpecies(i) = .false.
            radSpeciesId(i) = 0
        end select
    end do
end subroutine

subroutine getGasInfo (cell, gasinfo)
    implicit none
    integer, intent(in) :: cell
    type(gasInfoType), intent(out) :: gasinfo
    real(dp) :: xg
    integer :: i
    gasinfo%T = T_g(cell)
    gasinfo%P = P_g(cell)
    do i=1, NMAX(0) ! mass to mole fraction
        if (isRadSpecies(i)) then
	    ! for other gas species, gas mole fraction is uesed
		    xg= X_g(cell,i)*MW_MIX_g(cell)/MW_g(i)
            gasinfo%C(radSpeciesId(i)) = xg
        endif
    end do
end subroutine

end module rad_gas_species
