
subroutine find_corresponding_AMR_cell(xpoint, ypoint, zpoint, npalev, namrx, namry, namrz, &
                                        patch_level, patchrx, patchry, patchrz, &
                                        patchnx, patchny, patchnz, npatches_inside, &
                                        which_patches, L, ncoarse, solap, cr0amr1, &
                                        patch_id, amr_ix, amr_iy, amr_iz)
implicit none

!input
integer :: npalev, namrx, namry, namrz, npatches_inside, ncoarse
real :: L
integer, dimension(npalev) :: patchnx, patchny, patchnz, patch_level
real, dimension(npalev) :: patchrx, patchry, patchrz
integer, dimension(npatches_inside) :: which_patches
integer, dimension(namrx, namry, namrz, npalev) :: solap, cr0amr1
real(8) :: xpoint, ypoint, zpoint
!output
integer :: patch_id, amr_ix, amr_iy, amr_iz
!iterators
integer :: ipatch,iipatch
!local variables
real :: patch_rx, patch_ry, patch_rz, patch_res
integer :: patch_l
integer :: patch_nx, patch_ny, patch_nz


!Find the finest patch that contains the point
patch_id = 0 !coarse grid
do ipatch=npatches_inside,1,-1 !start by the finest patch
    iipatch = which_patches(ipatch)
    patch_l = patch_level(iipatch)
    patch_nx = patchnx(iipatch)
    patch_ny = patchny(iipatch)
    patch_nz = patchnz(iipatch)
    patch_rx = patchrx(iipatch)
    patch_ry = patchry(iipatch)
    patch_rz = patchrz(iipatch)
    patch_res = L / ncoarse / 2.0**(patch_l)
    !check if point is inside patch limits
    if ( (xpoint > patch_rx - patch_res) .and. (xpoint < patch_rx + patch_res*(patch_nx-1)) .and. &
         (ypoint > patch_ry - patch_res) .and. (ypoint < patch_ry + patch_res*(patch_ny-1)) .and. &
         (zpoint > patch_rz - patch_res) .and. (zpoint < patch_rz + patch_res*(patch_nz-1))   ) then

        !Find the AMR cell of the patch that contains the point
        amr_ix = int((xpoint - patch_rx) / patch_res + 1) + 1
        amr_iy = int((ypoint - patch_ry) / patch_res + 1) + 1
        amr_iz = int((zpoint - patch_rz) / patch_res + 1) + 1

        !Check if it is the master cell
        if (solap(amr_ix, amr_iy, amr_iz, iipatch) == 1  .and. cr0amr1(amr_ix, amr_iy, amr_iz, iipatch) == 1) then
            patch_id = iipatch
            exit
        endif
    end if
end do

!Calculate amr_ix, amr_iy, amr_iz if the point is in the coarse grid
if (patch_id == 0) then
    patch_res = L / ncoarse
    patch_rx = - L/2 
    patch_ry = - L/2 
    patch_rz = - L/2 
    amr_ix = int((xpoint - patch_rx) / patch_res) + 1
    amr_iy = int((ypoint - patch_ry) / patch_res) + 1
    amr_iz = int((zpoint - patch_rz) / patch_res) + 1
end if
end subroutine




subroutine SZ_effect(nx, ny, nz, res, grid_centers_x, grid_centers_y, grid_centers_z, &
                    npalev, namrx, namry, namrz, patch_level, patchrx, patchry, patchrz, &
                    patchnx, patchny, patchnz, npatches_inside, which_patches, &
                    L, ncoarse, solap, cr0amr1, U1, U11, U2, U21, U3, U31, U4, U41, temp, temp1, &
                    cx, cy, cz, radius, &
                    bulkVx, bulkVy, bulkVz, &
                    zeta8, re0, &
                    local_kSZ_x, local_kSZ_y, local_kSZ_z, tSZ_x, tSZ_y, tSZ_z, &
                    global_kSZ_x, global_kSZ_y, global_kSZ_z, &
                    dens_x, dens_y, dens_z, los_vx, los_vy, los_vz, &
                    temp_cutoff, dens_cutoff, &
                    T_SZ_cells, dens_SZ_cells, thrash_mass_fraction)

use omp_lib
implicit none

!input, output
integer :: nx, ny, nz
integer :: npalev, namrx, namry, namrz, npatches_inside, ncoarse
real :: L, temp_cutoff, dens_cutoff
integer, dimension(npalev) :: patchnx, patchny, patchnz, patch_level
real, dimension(npalev) :: patchrx, patchry, patchrz
integer, dimension(npatches_inside) :: which_patches
integer, dimension(namrx, namry, namrz, npalev) :: solap, cr0amr1
real, dimension(ncoarse, ncoarse, ncoarse) :: U1, U2, U3, U4, temp
real, dimension(namrx, namry, namrz, npalev) :: U11, U21, U31, U41, temp1
real :: cx, cy, cz, radius
real(8) :: bulkVx, bulkVy, bulkVz
real(8), dimension(nx) :: grid_centers_x
real(8), dimension(ny) :: grid_centers_y
real(8), dimension(nz) :: grid_centers_z
real(8) :: res
real(8) :: zeta8
real(8) :: re0
real(8), dimension(ny, nz) :: local_kSZ_x, global_kSZ_x, tSZ_x, dens_x, los_vx
real(8), dimension(nx, nz) :: local_kSZ_y, global_kSZ_y, tSZ_y, dens_y, los_vy
real(8), dimension(nx, ny) :: local_kSZ_z, global_kSZ_z, tSZ_z, dens_z, los_vz
real(8), dimension(nx*ny*nz) :: T_SZ_cells, dens_SZ_cells
real(8) :: thrash_mass_fraction
!constants
real(8), parameter :: mpc2m = 3.0856e22 !Mpc to m
real(8), parameter :: c = 2.9979e8 !m/s
real(8), parameter :: kb = 1.3806e-23 !m^2 kg s^-2 K^-1
real(8), parameter :: me = 9.1094e-31 !kg
real(8), parameter :: sigma_T = 6.6524e-29 !m^2
real(8), parameter :: factor_kSZ = sigma_T * mpc2m
real(8), parameter :: factor_tSZ = kb*sigma_T/(me*c**2) * mpc2m
real(8), parameter :: mu = 0.6 !mean molecular weight, for ionized gas
real(8), parameter :: mp = 1.6726e-27 !kg
real(8), parameter :: density_to_isu = 4.6693e-22 !kg/m^3

!local
integer :: ix, jy, kz
real(8) :: thrashVol_x, thrashMass_x, volTot_x, massTot_x

!ray variables
real(8) :: xpoint, ypoint, zpoint
integer :: patch_id, amr_ix, amr_iy, amr_iz
real(8) :: dens, T, LOSvel

thrashVol_x = 0.0
thrashMass_x = 0.0
volTot_x = 0.0
massTot_x = 0.0

bulkVx = bulkVx/(c/1.e3) !to c units
bulkVy = bulkVy/(c/1.e3)
bulkVz = bulkVz/(c/1.e3)

!X direction
!$OMP PARALLEL DEFAULT(SHARED), &
!$OMP PRIVATE(jy, kz, zpoint, ypoint, xpoint, patch_id, amr_ix, amr_iy, amr_iz, dens, T, LOSvel)
!$OMP DO REDUCTION(+:global_kSZ_x, local_kSZ_x, tSZ_x, los_vx, dens_x, volTot_x, massTot_x, thrashVol_x, thrashMass_x, &
!$OMP T_SZ_cells, dens_SZ_cells)
do kz=1,nz
    zpoint = grid_centers_z(kz)
    do jy=1,ny
        ypoint = grid_centers_y(jy)
        global_kSZ_x(jy,kz) = 0.0
        local_kSZ_x(jy,kz) = 0.0
        tSZ_x(jy,kz) = 0.0
        do ix = 1,nx ! ray in x direction
            xpoint = grid_centers_x(ix)
            call find_corresponding_AMR_cell(xpoint, ypoint, zpoint, npalev, namrx, namry, namrz, &
                                            patch_level, patchrx, patchry, patchrz, &
                                            patchnx, patchny, patchnz, npatches_inside, &
                                            which_patches, L, ncoarse, solap, cr0amr1, &
                                            patch_id, amr_ix, amr_iy, amr_iz)
            
            if (patch_id == 0) then !coarse grid
                dens = U1(amr_ix, amr_iy, amr_iz)
                T = temp(amr_ix, amr_iy, amr_iz)
                LOSvel = U2(amr_ix, amr_iy, amr_iz)
            else
                dens = U11(amr_ix, amr_iy, amr_iz, patch_id)
                T = temp1(amr_ix, amr_iy, amr_iz, patch_id)
                LOSvel = U21(amr_ix, amr_iy, amr_iz, patch_id)
            endif

            !PHASE DIAGRAM
            T_SZ_cells(ix + (jy-1)*nx + (kz-1)*nx*ny) = T
            dens_SZ_cells(ix + (jy-1)*nx + (kz-1)*nx*ny) = dens

            !check density and temperature (to prove that it is indeed ICM)
            !control for avoided cells
            if (dens > dens_cutoff*(1+zeta8)**3 .or. T < temp_cutoff) then
                thrashVol_x = thrashVol_x + res**3
                thrashMass_x = thrashMass_x + dens * res**3
            endif
            volTot_x = volTot_x + res**3
            massTot_x = massTot_x + dens * res**3

            ! !avoid cluster center for unusually high values 
            ! if ((xpoint-cx)**2 + (ypoint-cy)**2 + (zpoint-cz)**2 < (0.1*radius)**2) then
            !     cycle
            ! endif

            !density threshold for avoiding unphyisical clumps
            if (dens > dens_cutoff*(1+zeta8)**3) then
                dens = dens_cutoff*(1+zeta8)**3
            endif
            
            !tempreature threshold for selecting just ICM (ionized)
            if (T < temp_cutoff) then
                cycle
            endif

            global_kSZ_x(jy,kz) = global_kSZ_x(jy,kz) + dens * LOSvel * res
            local_kSZ_x(jy,kz) = local_kSZ_x(jy,kz) + dens * (LOSvel - bulkVx) * res
            tSZ_x(jy,kz) = tSZ_x(jy,kz) + dens * T * res
            los_vx(jy,kz) = los_vx(jy,kz) + LOSvel * res
            dens_x(jy,kz) = dens_x(jy,kz) + dens * res

        enddo
    enddo
enddo
!$OMP END DO
!$OMP END PARALLEL

!Y direction
!$OMP PARALLEL DEFAULT(SHARED), &
!$OMP PRIVATE(ix, kz, zpoint, ypoint, xpoint, patch_id, amr_ix, amr_iy, amr_iz, dens, T, LOSvel)
!$OMP DO REDUCTION(+:global_kSZ_y, local_kSZ_y, tSZ_y, los_vy, dens_y)
do kz=1,nz
    zpoint = grid_centers_z(kz)
    do ix=1,nx
        xpoint = grid_centers_x(ix)
        global_kSZ_y(ix,kz) = 0.0
        local_kSZ_y(ix,kz) = 0.0
        tSZ_y(ix,kz) = 0.0
        do jy = 1,ny ! ray in y direction
            ypoint = grid_centers_y(jy)
            call find_corresponding_AMR_cell(xpoint, ypoint, zpoint, npalev, namrx, namry, namrz, &
                                            patch_level, patchrx, patchry, patchrz, &
                                            patchnx, patchny, patchnz, npatches_inside, &
                                            which_patches, L, ncoarse, solap, cr0amr1, &
                                            patch_id, amr_ix, amr_iy, amr_iz)

            if (patch_id == 0) then !coarse grid
                dens = U1(amr_ix, amr_iy, amr_iz)
                T = temp(amr_ix, amr_iy, amr_iz)
                LOSvel = U3(amr_ix, amr_iy, amr_iz)
            else
                dens = U11(amr_ix, amr_iy, amr_iz, patch_id)
                T = temp1(amr_ix, amr_iy, amr_iz, patch_id)
                LOSvel = U31(amr_ix, amr_iy, amr_iz, patch_id)
            endif

            ! !avoid cluster center for unusually high values 
            ! if ((xpoint-cx)**2 + (ypoint-cy)**2 + (zpoint-cz)**2 < (0.1*radius)**2) then
            !     cycle
            ! endif

            !check density and temperature (to prove that it is indeed ICM)
            if (dens > dens_cutoff*(1+zeta8)**3) then
                dens = dens_cutoff*(1+zeta8)**3
            endif

            if (T < temp_cutoff) then
                cycle
            endif

            global_kSZ_y(ix,kz) = global_kSZ_y(ix,kz) + dens * LOSvel * res
            local_kSZ_y(ix,kz) = local_kSZ_y(ix,kz) + dens * (LOSvel - bulkVy) * res
            tSZ_y(ix,kz) = tSZ_y(ix,kz) + dens * T * res
            los_vy(ix,kz) = los_vy(ix,kz) + LOSvel * res
            dens_y(ix,kz) = dens_y(ix,kz) + dens * res

        enddo
    enddo
enddo
!$OMP END DO
!$OMP END PARALLEL

!Z direction
!$OMP PARALLEL DEFAULT(SHARED), &
!$OMP PRIVATE(ix, jy, zpoint, ypoint, xpoint, patch_id, amr_ix, amr_iy, amr_iz, dens, T, LOSvel)
!$OMP DO REDUCTION(+:global_kSZ_z, local_kSZ_z, tSZ_z, los_vz, dens_z)
do jy=1,ny
    ypoint = grid_centers_y(jy)
    do ix=1,nx
        xpoint = grid_centers_x(ix)
        global_kSZ_z(ix,jy) = 0.0
        local_kSZ_z(ix,jy) = 0.0
        tSZ_z(ix,jy) = 0.0
        do kz=1,nz
            zpoint = grid_centers_z(kz)
            call find_corresponding_AMR_cell(xpoint, ypoint, zpoint, npalev, namrx, namry, namrz, &
                                            patch_level, patchrx, patchry, patchrz, &
                                            patchnx, patchny, patchnz, npatches_inside, &
                                            which_patches, L, ncoarse, solap, cr0amr1, &
                                            patch_id, amr_ix, amr_iy, amr_iz)
                                            
            if (patch_id == 0) then !coarse grid
                dens = U1(amr_ix, amr_iy, amr_iz)
                T = temp(amr_ix, amr_iy, amr_iz)
                LOSvel = U4(amr_ix, amr_iy, amr_iz)
            else
                dens = U11(amr_ix, amr_iy, amr_iz, patch_id)
                T = temp1(amr_ix, amr_iy, amr_iz, patch_id)
                LOSvel = U41(amr_ix, amr_iy, amr_iz, patch_id)
            endif

            ! !avoid cluster center for unusually high values 
            ! if ((xpoint-cx)**2 + (ypoint-cy)**2 + (zpoint-cz)**2 < (0.1*radius)**2) then
            !     cycle
            ! endif

            !check density and temperature (to prove that it is indeed ICM)
            if (dens > dens_cutoff*(1+zeta8)**3) then
                dens = dens_cutoff*(1+zeta8)**3
            endif

            if (T < temp_cutoff) then
                cycle
            endif


            global_kSZ_z(ix,jy) = global_kSZ_z(ix,jy) + dens * LOSvel * res
            local_kSZ_z(ix,jy) = local_kSZ_z(ix,jy) + dens * (LOSvel - bulkVz) * res
            tSZ_z(ix,jy) = tSZ_z(ix,jy) + dens * T * res
            los_vz(ix,jy) = los_vz(ix,jy) + LOSvel * res
            dens_z(ix,jy) = dens_z(ix,jy) + dens * res

        enddo
    enddo
enddo
!$OMP END DO
!$OMP END PARALLEL

!First, density to ISU
global_kSZ_x = global_kSZ_x * density_to_isu
global_kSZ_y = global_kSZ_y * density_to_isu
global_kSZ_z = global_kSZ_z * density_to_isu
local_kSZ_x = local_kSZ_x * density_to_isu
local_kSZ_y = local_kSZ_y * density_to_isu
local_kSZ_z = local_kSZ_z * density_to_isu
tSZ_x = tSZ_x * density_to_isu
tSZ_y = tSZ_y * density_to_isu
tSZ_z = tSZ_z * density_to_isu
dens_x = dens_x * density_to_isu
dens_y = dens_y * density_to_isu
dens_z = dens_z * density_to_isu

!Second, mass density to number density of electrons
global_kSZ_x = global_kSZ_x / (mu*mp)
global_kSZ_y = global_kSZ_y / (mu*mp)
global_kSZ_z = global_kSZ_z / (mu*mp)
local_kSZ_x = local_kSZ_x / (mu*mp)
local_kSZ_y = local_kSZ_y / (mu*mp)
local_kSZ_z = local_kSZ_z / (mu*mp)
tSZ_x = tSZ_x / (mu*mp)
tSZ_y = tSZ_y / (mu*mp)
tSZ_z = tSZ_z / (mu*mp)
dens_x = dens_x / (mu*mp)
dens_y = dens_y / (mu*mp)
dens_z = dens_z / (mu*mp)

!Third, multiply by the factors outside the integral
global_kSZ_x = global_kSZ_x * factor_kSZ
global_kSZ_y = global_kSZ_y * factor_kSZ 
global_kSZ_z = global_kSZ_z * factor_kSZ 
local_kSZ_x = local_kSZ_x * factor_kSZ 
local_kSZ_y = local_kSZ_y * factor_kSZ 
local_kSZ_z = local_kSZ_z * factor_kSZ 
tSZ_x = tSZ_x * factor_tSZ
tSZ_y = tSZ_y * factor_tSZ
tSZ_z = tSZ_z * factor_tSZ

!Last, from internal units to Mpc
global_kSZ_x = global_kSZ_x / (1.0+zeta8)
global_kSZ_y = global_kSZ_y / (1.0+zeta8)
global_kSZ_z = global_kSZ_z / (1.0+zeta8)
local_kSZ_x = local_kSZ_x / (1.0+zeta8)
local_kSZ_y = local_kSZ_y / (1.0+zeta8)
local_kSZ_z = local_kSZ_z / (1.0+zeta8)
tSZ_x = tSZ_x / (1.0+zeta8)
tSZ_y = tSZ_y / (1.0+zeta8)
tSZ_z = tSZ_z / (1.0+zeta8)

!Thrash mass fraction to analyze the effect of the density and temperature thresholds
thrash_mass_fraction = thrashMass_x/massTot_x


! write(*,*) 'Thrash mass fraction', thrashMass_x/massTot_x, thrashMass_x, massTot_x
! write(*,*) 'Thrash volume fraction', thrashVol_x/volTot_x, thrashVol_x, volTot_x

end subroutine