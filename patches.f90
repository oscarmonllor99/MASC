subroutine which_patches_inside_box(box, lado0, npatch, patchnx, patchny, patchnz, patchrx, patchry, patchrz, & 
                                    patch_level, patch_inside)
implicit none
include 'cluster_SZ_parameters.dat'

integer ::  npatch(0:nlevels)
integer, dimension(npalev) :: patch_level
integer, dimension(npalev) :: patchnx, patchny, patchnz
real, dimension(npalev) :: patchrx, patchry, patchrz
integer, dimension(npalev) :: patch_inside
real(8), dimension(6) :: box
real :: lado0
real(8) :: cellsize
real(8), dimension(6) :: patch_boundaries
integer :: i, nx, ny, nz
real :: rx,ry,rz

do i=1,sum(npatch)
    !Get patch coordinates and size
    rx = patchrx(i)
    ry = patchry(i)
    rz = patchrz(i)
    nx = patchnx(i)
    ny = patchny(i)
    nz = patchnz(i)

    !Calculate patch boundaries
    cellsize = lado0 / nmax / 2.0**patch_level(i)
    patch_boundaries(1) = rx - cellsize
    patch_boundaries(2) = rx + cellsize*nx
    patch_boundaries(3) = ry - cellsize
    patch_boundaries(4) = ry + cellsize*ny
    patch_boundaries(5) = rz - cellsize
    patch_boundaries(6) = rz + cellsize*nz

    !Check if patch is inside box
    if ( (patch_boundaries(1) .le. box(2)) .and. (patch_boundaries(2) .ge. box(1)) .and. &
         (patch_boundaries(3) .le. box(4)) .and. (patch_boundaries(4) .ge. box(3)) .and. &
         (patch_boundaries(5) .le. box(6)) .and. (patch_boundaries(6) .ge. box(5)) ) then
      
        patch_inside(i) = 1

    else
        
        patch_inside(i) = 0

    endif

enddo
end subroutine



subroutine create_patch_levels(npalev, nlevels, npatch, patch_level)
implicit none 
integer :: npalev, nlevels
integer :: npatch(0:nlevels) !from 0 to nlevels, not from 1 to nlevels
integer :: lvl, i
integer, dimension(npalev) :: patch_level
patch_level = 0
lvl = 1
do i=1,sum(npatch)
    if (i .gt. sum(npatch(0:lvl))) then
        lvl = lvl + 1
    endif
    patch_level(i) = lvl
enddo
end subroutine


subroutine create_patch_ids(npalev, patchid)
implicit none
integer :: npalev, i
integer, dimension(npalev) :: patchid
patchid = 0
do i=1,npalev
    patchid(i) = i
enddo
end subroutine