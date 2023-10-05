
!***********************************************************************
!*      SUMMARY OF MASCLET UNITS                             *
!***********************************************************************
!!!! MAIN UNITS
! 1 u.l. = 10.98 Mpc
! 1 u.m. = 9.1717e18 Msun
! 1 u.t. = 1.129e15 sec
!!!! Therefore
! c = 1
! G = 1/(8pi)
! H = 3.66e-3 h (u.t.^-1  or   u.v./u.l.)
! 1 u.dens. = 4.72e-21 kg/m^3

!***********************************************************************
!*      SUBROUTINES IN EXTERNAL FILES                                  *
!***********************************************************************
include 'reader.f90'
include 'patches.f90'
include 'SZ.f90'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program cluster_SZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COMPUTES THE SUNYAEV-ZELDOVICH EFFECT FOR A CLUSTER OF GALAXIES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use omp_lib
implicit none
include 'cluster_SZ_parameters.dat'

!GAS
real, dimension(nmax, nmay, nmaz) :: U1, U2, U3, U4, temp
real, dimension(namrx, namry, namrz, npalev) :: U11, U21, U31, U41, temp1
integer, dimension(nmax, nmay, nmaz) :: cr0amr
integer, dimension(namrx, namry, namrz, npalev) :: cr0amr1, solap

!GRID
integer :: npatch(0:nlevels) !from 0 to nlevels, not from 1 to nlevels
integer, dimension(npalev) :: patch_level
integer, dimension(npalev) :: pare, patchnx, patchny, patchnz, patchx, patchy, patchz
real, dimension(npalev) :: patchrx, patchry, patchrz

!ASOHF HALOES
integer, dimension(nmaxclus) :: halID
real, dimension(nmaxclus) :: halRx, halRy, halRz, halVirMass, halVirRad, halVx, halVy, halVz
integer :: nhaloes
real(8), dimension(6) :: box
real(8), dimension(2) :: xlims, ylims, zlims
real(8) :: res
integer ::  nx, ny, nz
real(8), allocatable, dimension(:) :: grid_faces_x, grid_faces_y, grid_faces_z, grid_centers_x, grid_centers_y, grid_centers_z
integer, dimension(npalev) :: patch_inside
integer, dimension(npalev) :: patchid
integer :: npatches_inside
integer, dimension(:), allocatable :: which_patches
real(8) :: bulkVx, bulkVy, bulkVz

!PARAMETERS
character(len=5) iter_string
character(len=5) ih_string
integer :: first, last, every, iter, is_magnetic_field
real :: zeta
real(8) :: ache, omega0, fdm, zeta8, rho_background
real :: lado0
real :: cio_mass, cio_lenght, cio_speed
real :: thres_mass, factor_Rvir

!GAS DENSITY
real(8), parameter :: mu = 0.6 !mean molecular weight, for ionized gas
real(8), parameter :: mp = 1.6726e-27 !kg
real(8), parameter :: density_from_masclet2isu = 4.72e-21 !kg/m^3

!OUTPUT
real(8), allocatable, dimension(:,:) :: local_kSZ_x, global_kSZ_x, local_kSZ_y, global_kSZ_y, local_kSZ_z, global_kSZ_z
real(8), allocatable, dimension(:,:) :: tSZ_x, tSZ_y, tSZ_z

!LOOP INDICES
integer :: i, j, k, ix, jy, kz, ih

!**************************************************************
!*      OPENING FILES
!**************************************************************
OPEN(1, file = './cluster_SZ.dat', STATUS='UNKNOWN', ACTION='READ')

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! READING INITIAL DATA
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
READ(1,*) 
READ(1,*)
READ(1,*) 
READ(1,*) !Files: first, last, every -------------------------------------------->
READ(1,*) first,last,every
READ(1,*) !Cells per direction (NX,NY,NZ) --------------------------------------->
READ(1,*) !NX,NY,NZ
READ(1,*) !Hubble constant (h), omega matter, fraction of DM to total mass ------>
READ(1,*) ache,omega0,fdm
READ(1,*) !Max box sizelength (in length units specified below) ----------------->
READ(1,*) lado0
READ(1,*) !Reading: is_magnetic_field (=0, no; =1, yes) ------------------------->
READ(1,*) is_magnetic_field
READ(1,*) !Input units: MASS (Msun), LENGTH (cMpc), SPEED (km/s) ------------ --->
READ(1,*) cio_mass,cio_lenght,cio_speed
READ(1,*) !Cluster mass threshold (Msun) ---------------------------------------->
READ(1,*) thres_mass
READ(1,*) !Cluster box half sidelength factor (cMpc) (to multiply Rvir) --------->
READ(1,*) factor_Rvir
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
CLOSE(1)
!///////////////////////////////////
!///////////////////////////////////
!///////////////////////////////////
! MAIN DO LOOP (over iterations)
!///////////////////////////////////
do iter = first, last, every

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! READING THE cluster CATALOGUE
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! FIRST, INITIALIZE VARIABLES
halID = 0
halRx = 0.0
halRy = 0.0
halRz = 0.0
halVirMass = 0.0
halVirRad = 0.0
nhaloes = 0
! CALL READER
call read_ASOHF_catalogue(iter, nhaloes, halID, halRx, halRy, halRz, halVirMass, halVirRad, &
                          halVx, halVy, halVz)
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! READING grid AND clus FILES
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

write(*,*) 'reading iter', iter

! FIRST, INITIALIZE VARIABLES
patchnx = 0
patchny = 0
patchnz = 0
patchrx = 0.0
patchry = 0.0
patchrz = 0.0
patchx = 0
patchy = 0
patchz = 0
pare = 0
npatch = 0

call create_patch_ids(npalev, patchid)

do k = 1, nmaz
do j = 1, nmay
do i = 1, nmax
    U1(i,j,k) = -1.0
    U2(i,j,k) = 0.0
    U3(i,j,k) = 0.0
    U4(i,j,k) = 0.0
    temp(i,j,k) = 0.0
    cr0amr(i,j,k) = 0
enddo
enddo
enddo

do i=1,npalev
 do kz=1,namrz
 do jy=1,namry
 do ix=1,namrx
    U11(ix,jy,kz,i) = -1.0
    U21(ix,jy,kz,i) = 0.0
    U31(ix,jy,kz,i) = 0.0
    U41(ix,jy,kz,i) = 0.0
    temp1(ix,jy,kz,i) = 0.0
    cr0amr1(ix,jy,kz,i) = 0
    solap(ix,jy,kz,i) = 0
 enddo
 enddo
 enddo
enddo

! NOW READ THE clus FILE
call read_clus(iter, is_magnetic_field, U1, U2, U3, U4, temp, cr0amr, &
                U11, U21, U31, U41, temp1, solap, cr0amr1, &
                patchnx, patchny, patchnz, patchx, patchy, patchz, &
                patchrx, patchry, patchrz, pare, npatch, zeta)
            
! CREATE patch_level ARRAY
call create_patch_levels(npalev, nlevels, npatch, patch_level)



! BACKGROUND DENSITY of the UNIVERSE at this redshift
zeta8 = dble(zeta)
rho_background = 3.D0 * omega0 * (ache*3.66D-3)**2 * (1.D0+zeta8)**3 ! in u.mass / u.lenght^3

! U1 and U11 (density_contrast) to density
U1 = (1. + U1) * rho_background
U11 = (1. + U11) * rho_background
! Then conversion from density in units of background density to density in units of kg/m^3
U1 = U1 * density_from_masclet2isu
U11 = U11 * density_from_masclet2isu
! From gas density to number density of electrons
U1 = U1 / (mu*mp)
U11 = U11 / (mu*mp)

! write(*,*) 'max density', maxval(abs(U11))
! write(*,*) 'max vx', maxval(abs(U21))
! write(*,*) 'max vy', maxval(abs(U31))
! write(*,*) 'max vz', maxval(abs(U41))
! write(*,*) 'max patch level', maxval(patch_level)


! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! FINISHED READING  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  
 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  
! ! ! ! ! LOOP OVER CLUSTERS ! ! ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! !  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  

write(*,*) 'Calculating the S-Z effect for desired clusters...'

do ih=1,nhaloes

    if (halVirMass(ih) < thres_mass) cycle

    ! create the box for the cluster
    box(1) = halRx(ih) - factor_Rvir*halVirRad(ih)
    box(2) = halRx(ih) + factor_Rvir*halVirRad(ih)
    box(3) = halRy(ih) - factor_Rvir*halVirRad(ih)
    box(4) = halRy(ih) + factor_Rvir*halVirRad(ih)
    box(5) = halRz(ih) - factor_Rvir*halVirRad(ih)
    box(6) = halRz(ih) + factor_Rvir*halVirRad(ih)

    !Define the uniform grid
    res = (lado0/nmax)/2**nlevels
    xlims(1) = box(1)
    xlims(2) = box(2) + (res - mod(box(2) - box(1), res) ) !Number of cells between the box boundaries must be an integer, thus we modify the right boundary
    ylims(1) = box(3)
    ylims(2) = box(4) + (res - mod(box(4) - box(3), res) )
    zlims(1) = box(5)
    zlims(2) = box(6) + (res - mod(box(6) - box(5), res) )

    nx = int((xlims(2) - xlims(1))/res)
    ny = int((ylims(2) - ylims(1))/res)
    nz = int((zlims(2) - zlims(1))/res)
    
    ! Allocate and initialize the uniform grid
    allocate(grid_faces_x(nx+1), grid_faces_y(ny+1), grid_faces_z(nz+1), grid_centers_x(nx), grid_centers_y(ny), grid_centers_z(nz))

    do i=1,nx+1
        grid_faces_x(i) = xlims(1) + (i-1)*res
    enddo
    do i=1,ny+1
        grid_faces_y(i) = ylims(1) + (i-1)*res
    enddo
    do i=1,nz+1
        grid_faces_z(i) = zlims(1) + (i-1)*res
    enddo

    grid_centers_x = (grid_faces_x(1:nx) + grid_faces_x(2:nx+1))/2.0
    grid_centers_y = (grid_faces_y(1:ny) + grid_faces_y(2:ny+1))/2.0
    grid_centers_z = (grid_faces_z(1:nz) + grid_faces_z(2:nz+1))/2.0

    ! Search for the patches that are inside the box
    patch_inside = 0
    call which_patches_inside_box(box, lado0, npatch, patchnx, patchny, patchnz, &
                                  patchrx, patchry, patchrz, patch_level, patch_inside)
                                  
    npatches_inside = count(patch_inside /= 0)

    ! Now, parallel loop over the patches inside the box                
    allocate(which_patches(npatches_inside))
    which_patches = pack(patchid, patch_inside /= 0)

    write(*,*) '    Now, calculating the SZ effect...'

    ! SUNYAEV-ZELDOVICH EFFECT calculation
    allocate(local_kSZ_x(ny, nz), local_kSZ_y(nx, nz), local_kSZ_z(nx, ny))
    allocate(global_kSZ_x(ny, nz), global_kSZ_y(nx, nz), global_kSZ_z(nx, ny))
    allocate(tSZ_x(ny, nz), tSZ_y(nx, nz), tSZ_z(nx, ny))

    local_kSZ_x = 0.0
    local_kSZ_y = 0.0
    local_kSZ_z = 0.0
    global_kSZ_x = 0.0
    global_kSZ_y = 0.0
    global_kSZ_z = 0.0
    tSZ_x = 0.0
    tSZ_y = 0.0
    tSZ_z = 0.0

    bulkVx = dble(halVx(ih))
    bulkVy = dble(halVy(ih))
    bulkVz = dble(halVz(ih))

    call SZ_effect(nx, ny, nz, res, grid_centers_x, grid_centers_y, grid_centers_z, &
                    npalev, namrx, namry, namrz, patch_level, patchrx, patchry, patchrz, &
                    patchnx, patchny, patchnz, npatches_inside, which_patches, &
                    lado0, nmax, solap, U1, U11, U2, U21, U3, U31, U4, U41, temp, temp1, &
                    bulkVx, bulkVy, bulkVz, &
                    local_kSZ_x, local_kSZ_y, local_kSZ_z, tSZ_x, tSZ_y, tSZ_z, &
                    global_kSZ_x, global_kSZ_y, global_kSZ_z)

    write(*,*) '    Cluster', ih, 'finished , saving and deallocating arrays...'

    ! Save data in single precision (sngl)
    write(iter_string, '(i5.5)') iter
    write(ih_string, '(i5.5)') ih
    open(unit=33, file='output_files/SZ_effect'//'_it_'//iter_string//'_ih_'//ih_string, form='unformatted')

    !header
    write(33) nx, ny, nz 
    !sz effect
    write(33) ((sngl(local_kSZ_x(jy,kz)), jy=1,ny), kz=1,nz)
    write(33) ((sngl(local_kSZ_y(ix,kz)), ix=1,nx), kz=1,nz)
    write(33) ((sngl(local_kSZ_z(ix,jy)), ix=1,nx), jy=1,ny)
    write(33) ((sngl(global_kSZ_x(jy,kz)), jy=1,ny), kz=1,nz)
    write(33) ((sngl(global_kSZ_y(ix,kz)), ix=1,nx), kz=1,nz)
    write(33) ((sngl(global_kSZ_z(ix,jy)), ix=1,nx), jy=1,ny)
    write(33) ((sngl(tSZ_x(jy,kz)), jy=1,ny), kz=1,nz)
    write(33) ((sngl(tSZ_y(ix,kz)), ix=1,nx), kz=1,nz)
    write(33) ((sngl(tSZ_z(ix,jy)), ix=1,nx), jy=1,ny)
    close(33)

    ! Deallocation
    deallocate(grid_faces_x, grid_faces_y, grid_faces_z, grid_centers_x, grid_centers_y, grid_centers_z)
    deallocate(which_patches)
    deallocate(local_kSZ_x, local_kSZ_y, local_kSZ_z, global_kSZ_x, global_kSZ_y, global_kSZ_z)
    deallocate(tSZ_x, tSZ_y, tSZ_z)

enddo !loop over clusters

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 


!//// END MAIN DO LOOP /////////////
enddo
!///////////////////////////////////
!///////////////////////////////////
!///////////////////////////////////

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!  END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program cluster_SZ



