subroutine gas_core_bulk_velocity(npalev, namrx, namry, namrz, &
                                patch_level, patchrx, patchry, patchrz, &
                                patchnx, patchny, patchnz, npatches_inside, &
                                which_patches, L, ncoarse, solap, cr0amr1, &
                                U11, U21, U31, U41, &
                                x, y, z, R_core, vx, vy, vz, mass)
implicit none
!input
integer :: npalev, namrx, namry, namrz, npatches_inside, ncoarse
real :: L, R_core, x, y, z
integer, dimension(npalev) :: patchnx, patchny, patchnz, patch_level
real, dimension(npalev) :: patchrx, patchry, patchrz
integer, dimension(npatches_inside) :: which_patches
integer, dimension(namrx, namry, namrz, npalev) :: solap, cr0amr1
real, dimension(namrx, namry, namrz, npalev) :: U11, U21, U31, U41
!output
real(8) :: vx, vy, vz, mass
!iterators
integer :: ipatch,iipatch, ix, iy, iz
real :: xcell, ycell, zcell, cell_volume, cell_mass
!local variables
real :: patch_rx, patch_ry, patch_rz, patch_res
integer :: patch_l
integer :: patch_nx, patch_ny, patch_nz

!U11 must be density, not delta
do ipatch=1, npatches_inside
    iipatch = which_patches(ipatch)
    patch_l = patch_level(iipatch)
    patch_nx = patchnx(iipatch)
    patch_ny = patchny(iipatch)
    patch_nz = patchnz(iipatch)
    patch_rx = patchrx(iipatch)
    patch_ry = patchry(iipatch)
    patch_rz = patchrz(iipatch)
    patch_res = L / ncoarse / 2.0**(patch_l)
    do ix=1,patch_nx
        xcell = patch_rx - 0.5*patch_res + (ix-1)*patch_res
        if ((x - xcell)**2 < R_core**2) then
            do iy=1,patch_ny
                ycell = patch_ry - 0.5*patch_res + (iy-1)*patch_res
                if ((y - ycell)**2 < R_core**2) then
                    do iz=1,patch_nz
                        zcell = patch_rz - 0.5*patch_res + (iz-1)*patch_res
                        if ((x - xcell)**2 + (y - ycell)**2 + (z - zcell)**2 < R_core**2) then
                            cell_volume = patch_res**3
                            !if it is refined or solaped, then cell_mass is zero and it doesn't matter
                            cell_mass = U11(ix,iy,iz,iipatch) * cell_volume * cr0amr1(ix,iy,iz,iipatch) * solap(ix,iy,iz,iipatch)
                            vx = vx + U21(ix,iy,iz,iipatch) * cell_mass
                            vy = vy + U31(ix,iy,iz,iipatch) * cell_mass
                            vz = vz + U41(ix,iy,iz,iipatch) * cell_mass
                            mass = mass + cell_mass
                        end if
                    end do
                end if
            end do
        end if
    end do
end do

vx = vx / mass
vy = vy / mass
vz = vz / mass

end subroutine