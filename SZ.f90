
subroutine find_corresponding_AMR_cell(xpoint, ypoint, zpoint, npalev, namrx, namry, namrz, &
                                        patch_level, patchrx, patchry, patchrz, &
                                        patchnx, patchny, patchnz, npatches_inside, &
                                        which_patches, L, ncoarse, solap, &
                                        patch_id, amr_ix, amr_iy, amr_iz)
implicit none

!input
integer :: npalev, namrx, namry, namrz, npatches_inside, ncoarse
real :: L
integer, dimension(npalev) :: patchnx, patchny, patchnz, patch_level
real, dimension(npalev) :: patchrx, patchry, patchrz
integer, dimension(npatches_inside) :: which_patches
integer, dimension(namrx, namry, namrz, npalev) :: solap
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
        amr_ix = int((xpoint - patch_rx) / patch_res + 0.5) + 1
        amr_iy = int((ypoint - patch_ry) / patch_res + 0.5) + 1
        amr_iz = int((zpoint - patch_rz) / patch_res + 0.5) + 1

        !Check if it is the master cell
        if (solap(amr_ix, amr_iy, amr_iz, iipatch) == 1) then
            patch_id = iipatch
            exit
        endif
    end if
end do

!Calculate amr_ix, amr_iy, amr_iz if the point is in the coarse grid
if (patch_id == 0) then
    patch_res = L / ncoarse
    patch_rx = -L/2 + patch_res
    patch_ry = -L/2 + patch_res
    patch_rz = -L/2 + patch_res
    amr_ix = int((xpoint - patch_rx) / patch_res + 0.5) + 1
    amr_iy = int((ypoint - patch_ry) / patch_res + 0.5) + 1
    amr_iz = int((zpoint - patch_rz) / patch_res + 0.5) + 1
end if
end subroutine




subroutine SZ_effect(nx, ny, nz, res, grid_centers_x, grid_centers_y, grid_centers_z, &
                    npalev, namrx, namry, namrz, patch_level, patchrx, patchry, patchrz, &
                    patchnx, patchny, patchnz, npatches_inside, which_patches, &
                    L, ncoarse, solap, U1, U11, U2, U21, U3, U31, U4, U41, temp, temp1, &
                    bulkVx, bulkVy, bulkVz, &
                    local_kSZ_x, local_kSZ_y, local_kSZ_z, tSZ_x, tSZ_y, tSZ_z, &
                    global_kSZ_x, global_kSZ_y, global_kSZ_z)

use omp_lib
implicit none

!input, output
integer :: nx, ny, nz
integer :: npalev, namrx, namry, namrz, npatches_inside, ncoarse
real :: L
integer, dimension(npalev) :: patchnx, patchny, patchnz, patch_level
real, dimension(npalev) :: patchrx, patchry, patchrz
integer, dimension(npatches_inside) :: which_patches
integer, dimension(namrx, namry, namrz, npalev) :: solap
real, dimension(ncoarse, ncoarse, ncoarse) :: U1, U2, U3, U4, temp
real, dimension(namrx, namry, namrz, npalev) :: U11, U21, U31, U41, temp1
real(8) :: bulkVx, bulkVy, bulkVz
real(8), dimension(nx) :: grid_centers_x
real(8), dimension(ny) :: grid_centers_y
real(8), dimension(nz) :: grid_centers_z
real(8) :: res
real(8), dimension(ny, nz) :: local_kSZ_x, global_kSZ_x, tSZ_x
real(8), dimension(nx, nz) :: local_kSZ_y, global_kSZ_y, tSZ_y
real(8), dimension(nx, ny) :: local_kSZ_z, global_kSZ_z, tSZ_z

!constants
real(8), parameter :: mpc2m = 3.0856e22 !Mpc to m
real(8), parameter :: c = 2.9979e8 !m/s
real(8), parameter :: kb = 1.3806e-23 !m^2 kg s^-2 K^-1
real(8), parameter :: me = 9.1094e-31 !kg
real(8), parameter :: sigma_T = 6.6524e-29 !m^2
real(8), parameter :: factor_kSZ = sigma_T * mpc2m
real(8), parameter :: factor_tSZ = kb*sigma_T/(me*c**2) * mpc2m
real(8), parameter :: max_ksz = 1e-4 !maximum value of a cell contribution to kSZ

!iterators
integer :: ix, jy, kz

!ray variables
real(8) :: xpoint, ypoint, zpoint
integer :: patch_id, amr_ix, amr_iy, amr_iz
real(8) :: dens, T, LOSvel

!X direction
!$OMP PARALLEL DEFAULT(SHARED), &
!$OMP PRIVATE(jy, kz, zpoint, ypoint, xpoint, patch_id, amr_ix, amr_iy, amr_iz, dens, T, LOSvel)
!$OMP DO REDUCTION(+:global_kSZ_x, local_kSZ_x, tSZ_x)
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
                                            which_patches, L, ncoarse, solap, &
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
            
            !control for large values (discontinuities) in kSZ
            if (factor_kSZ * dens * abs(LOSvel) * res > max_ksz) then
                global_kSZ_x(jy,kz) = global_kSZ_x(jy,kz) + sign(max_ksz, LOSvel)
            else
                global_kSZ_x(jy,kz) = global_kSZ_x(jy,kz) + factor_kSZ * dens * LOSvel * res
            endif

            if (factor_kSZ * dens * abs(LOSvel - bulkVx/(c/1.e3)) * res > max_ksz) then
                local_kSZ_x(jy,kz) = local_kSZ_x(jy,kz) + sign(max_ksz, LOSvel - bulkVx/(c/1.e3))
            else
                local_kSZ_x(jy,kz) = local_kSZ_x(jy,kz) + factor_kSZ * dens * (LOSvel - bulkVx/(c/1.e3)) * res
            endif

            tSZ_x(jy,kz) = tSZ_x(jy,kz) + factor_tSZ * dens * T * res

        enddo
    enddo
enddo
!$OMP END DO
!$OMP END PARALLEL

!Y direction
!$OMP PARALLEL DEFAULT(SHARED), &
!$OMP PRIVATE(ix, kz, zpoint, ypoint, xpoint, patch_id, amr_ix, amr_iy, amr_iz, dens, T, LOSvel)
!$OMP DO REDUCTION(+:global_kSZ_y, local_kSZ_y, tSZ_y)
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
                                            which_patches, L, ncoarse, solap, &
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

            !control for large values (discontinuities) in kSZ
            if (factor_kSZ * dens * abs(LOSvel) * res > max_ksz) then
                global_kSZ_y(ix,kz) = global_kSZ_y(ix,kz) + sign(max_ksz, LOSvel)
            else
                global_kSZ_y(ix,kz) = global_kSZ_y(ix,kz) + factor_kSZ * dens * LOSvel * res
            endif

            if (factor_kSZ * dens * abs(LOSvel - bulkVy/(c/1.e3)) * res > max_ksz) then
                local_kSZ_y(ix,kz) = local_kSZ_y(ix,kz) + sign(max_ksz, LOSvel - bulkVy/(c/1.e3))
            else
                local_kSZ_y(ix,kz) = local_kSZ_y(ix,kz) + factor_kSZ * dens * (LOSvel - bulkVy/(c/1.e3)) * res
            endif

            tSZ_y(ix,kz) = tSZ_y(ix,kz) + factor_tSZ * dens * T * res

        enddo
    enddo
enddo
!$OMP END DO
!$OMP END PARALLEL

!Z direction
!$OMP PARALLEL DEFAULT(SHARED), &
!$OMP PRIVATE(ix, jy, zpoint, ypoint, xpoint, patch_id, amr_ix, amr_iy, amr_iz, dens, T, LOSvel)
!$OMP DO REDUCTION(+:global_kSZ_z, local_kSZ_z, tSZ_z)
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
                                            which_patches, L, ncoarse, solap, &
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

            !control for large values (discontinuities) in kSZ
            if (factor_kSZ * dens * abs(LOSvel) * res > max_ksz) then
                global_kSZ_z(ix,jy) = global_kSZ_z(ix,jy) + sign(max_ksz, LOSvel)
            else
                global_kSZ_z(ix,jy) = global_kSZ_z(ix,jy) + factor_kSZ * dens * LOSvel * res
            endif

            if (factor_kSZ * dens * abs(LOSvel - bulkVz/(c/1.e3)) * res > max_ksz) then
                local_kSZ_z(ix,jy) = local_kSZ_z(ix,jy) + sign(max_ksz, LOSvel - bulkVz/(c/1.e3))
            else
                local_kSZ_z(ix,jy) = local_kSZ_z(ix,jy) + factor_kSZ * dens * (LOSvel - bulkVz/(c/1.e3)) * res
            endif

            tSZ_z(ix,jy) = tSZ_z(ix,jy) + factor_tSZ * dens * T * res

        enddo
    enddo
enddo
!$OMP END DO
!$OMP END PARALLEL


end subroutine