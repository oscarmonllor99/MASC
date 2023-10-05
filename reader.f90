subroutine read_clus(iter, is_magnetic_field, U1, U2, U3, U4, temp, cr0amr, &
                      U11, U21, U31, U41, temp1, solap, cr0amr1, &
                      patchnx, patchny, patchnz, patchx, patchy, patchz, &
                      patchrx, patchry, patchrz, pare, npatch, zeta)
!!!!! THIS IS A SUBROUTINE TO READ THE clus(BARYONIC) AND grid FILE
implicit none
include 'cluster_SZ_parameters.dat'

real, dimension(nmax, nmay, nmaz) :: U1, U2, U3, U4, temp
real, dimension(namrx, namry, namrz, npalev) :: U11, U21, U31, U41, temp1
integer, dimension(nmax, nmay, nmaz) :: cr0amr
integer, dimension(namrx, namry, namrz, npalev) :: cr0amr1, solap

integer :: npatch(0:nlevels) !from 0 to nlevels, not from 1 to nlevels
integer, dimension(npalev) :: pare, patchnx, patchny, patchnz, patchx, patchy, patchz
real, dimension(npalev) :: patchrx, patchry, patchrz

integer :: i, j, k, nl, irr, ir, ndm, low1, low2, ix, n1, n2, n3
real :: t, AAA, BBB, CCC, zeta, map

integer :: iter, is_magnetic_field
character(len=5) :: iter_string

write(*,*) 'reading grid file...'

!reading data
write(iter_string, '(i5.5)') iter
open(33, file='./simu_masclet/grids'//iter_string, status='unknown', action='read')
open(31, file='./simu_masclet/clus'//iter_string, status='unknown', action='read', form='unformatted')

! grid data
read(33,*) irr, t, nl, map
read(33,*) zeta
read(33,*) ir, ndm

do ir = 1,nl 
    read(33,*) irr, npatch(ir)
    read(33,*)
    if (ir .ne. irr) write(*,*) 'Warning: fail in restart'
    low1 = sum(npatch(0:ir-1))+1
    low2 = sum(npatch(0:ir))
    do i = low1,low2
        read(33,*) patchnx(i), patchny(i), patchnz(i)
        read(33,*) patchx(i), patchy(i), patchz(i)
        read(33,*) AAA, BBB, CCC
        patchrx(i) = AAA
        patchry(i) = BBB
        patchrz(i) = CCC
        read(33,*) pare(i) 
    enddo
enddo

close(33)
write(*,*) 'total number of patches = ', sum(npatch(0:nl))
write(*,*) '..done'

write(*,*) 'reading clus file...'

! clus data (baryonic)
read(31)
ir = 0
read(31) (((U1(i,j,k), i=1,nmax), j=1,nmay), k=1,nmaz)
read(31) (((U2(i,j,k), i=1,nmax), j=1,nmay), k=1,nmaz)
read(31) (((U3(i,j,k), i=1,nmax), j=1,nmay), k=1,nmaz)
read(31) (((U4(i,j,k), i=1,nmax), j=1,nmay), k=1,nmaz)
read(31) !pres
read(31) !pot
read(31) !opot
read(31) (((temp(i,j,k), i=1,nmax), j=1,nmay), k=1,nmaz)
read(31) !metalicity
read(31) (((cr0amr(i,j,k), i=1,nmax), j=1,nmay), k=1,nmaz)
if (is_magnetic_field == 1) then
    read(31) !Bx
    read(31) !By
    read(31) !Bz
end if

do ir=1,nl
    write(*,*) '    reading level', ir
    low1 = sum(npatch(0:ir-1))+1
    low2 = sum(npatch(0:ir))
    do i = low1,low2
        n1 = patchnx(i)
        n2 = patchny(i)
        n3 = patchnz(i)
        read(31) (((U11(ix,j,k,i), ix=1,n1), j=1,n2), k=1,n3)
        read(31) (((U21(ix,j,k,i), ix=1,n1), j=1,n2), k=1,n3)
        read(31) (((U31(ix,j,k,i), ix=1,n1), j=1,n2), k=1,n3)
        read(31) (((U41(ix,j,k,i), ix=1,n1), j=1,n2), k=1,n3)
        read(31) !pres
        read(31) !pot
        read(31) !opot
        read(31) (((temp1(ix,j,k,i), ix=1,n1), j=1,n2), k=1,n3)
        read(31) !metalicity
        read(31) (((cr0amr1(ix,j,k,i), ix=1,n1), j=1,n2), k=1,n3)
        read(31) (((solap(ix,j,k,i), ix=1,n1), j=1,n2), k=1,n3)
        if (is_magnetic_field == 1) then
            read(31) !Bx
            read(31) !By
            read(31) !Bz
        end if
    enddo
enddo
close(31)
write(*,*) '..done'

end subroutine 


subroutine read_ASOHF_catalogue(iter, nhaloes, halID, halRx, halRy, halRz, halVirMass, halVirRad, &
                                halVx, halVy, halVz)
implicit none
!!!!! THIS IS A SUBROUTINE TO READ THE ASOHF catalogue
include 'cluster_SZ_parameters.dat'
integer :: iter

integer, dimension(nmaxclus) :: halID
real, dimension(nmaxclus) :: halRx, halRy, halRz, halVirMass, halVirRad, halVx, halVy, halVz
integer :: nhaloes

integer :: i, BASi
real :: BASf

character(len=5) :: iter_string

write(iter_string, '(i5.5)') iter
open(20, file = './ASOHF_results/families'//iter_string, status = 'unknown', action = 'read')

! HEADER
read(20,*)
read(20,*) BASi, BASi, nhaloes
read(20,*)
read(20,*)
read(20,*)
read(20,*)
read(20,*)
do i=1,nhaloes
    read(20,*) halID(i), BASi, halRx(i), halRy(i), halRz(i), halVirMass(i), halVirRad(i), &
               BASf, BASf, BASi, BASi, BASf, BASf, BASf, BASf, BASf, BASf, BASf, &
               BASf, BASf, BASf, BASf, BASf, BASf, BASf, BASf, BASf, &
               halVx(i), halVy(i), halVz(i)
enddo

close(20)
end subroutine
