!
! Version 13.060
!

!
! Controls the size of the character variables
!
module charsize
  integer, parameter :: charsize1 = 6
end module charsize

! Function that prints the version

subroutine version()

write(*,*) 
write(*,*) ' Version 20.0.2'
write(*,*) 

return 
end

!
! This file contains common functions and subroutines for
! all programs. Must be compiled with each of them.
!
! L. Martinez, I. Pasteur, Jun 18, 2008.
!

!
! Gets keyword from input file
!

function keyword(string)

implicit none
integer :: if, il
character(len=200) :: keyword, string

if = 1
do while(string(if:if) <= ' '.and. if < 200) 
  if = if + 1
end do
il = if
do while(string(il:il) > ' '.and.il < 200)
  il = il + 1
end do
il = il - 1
keyword = string(if:il)

return
end    

!
! Gets keyword value from input file
!

function value(string)

implicit none
integer :: if, il, length
character(len=200) :: value, string

if = 0
do while(string(if:if) <= ' '.and.if < 200) 
  if = if + 1
end do
il = if
do while(string(il:il) > ' '.and.il < 200)
  il = il + 1
end do
il = il - 1
if = il + 1 
do while(string(if:if) <= ' '.and.if < 200) 
  if = if + 1
end do
il = if
do while(string(il:il) > ' '.and.il < 200)
  il = il + 1
end do
value = string(if:il)
if(length(value) == 0) then
  write(*,*) ' ERROR: Some keyword without value: '
  write(*,*) string(1:length(string))
  stop
end if

return
end     

!
! subroutine getdim: Simply return the number of atoms of the system
!                    as specified in the psf and the number of
!                    classes in parameter files to allocate arrays.
!
! On input: psffile: name of the psf file
!           inputfile: name of the input file
!
! On return: ndim: minimum required dimension
!

subroutine getdim(psffile,inputfile,ndim)   

implicit none
integer :: ndim, nclass, length, status
character(len=200) :: psffile, inputfile, record, keyword,&
                      value, file
character :: firstchar

open(10,file=psffile,action='read')
read(10,"( a200 )") record
do while(record(length(record)-5:length(record)).ne.'!NATOM')
  read(10,"( a200 )") record
end do
read(record,*) ndim
close(10)

nclass = 0
open(99,file=inputfile,action='read')
do while(.true.)
  read(99,"( a200 )",iostat=status) record
  if(status /= 0) exit
  if(keyword(record) == 'par') then
    file = value(record)
    record(1:9) = '#########'
    open(10,file=file,action='read',status='old',iostat=status)
    if(status /= 0) then
      write(*,*) ' ERROR: Parameter file not found: '
      write(*,*) trim(file)
      stop
    end if
    do while (record(1:9) /= 'NONBONDED')
      read(10,"( a200 )",iostat=status) record   
      if(status /= 0) then
        write(*,*) ' ERROR: Error reading parameter file (after NONBONDED mark): '
        write(*,*) trim(file)
        stop
      end if
    end do
    do while(.true.)
      read(10,"( a200 )",iostat=status) record
      if(status /= 0) exit
      if(firstchar(record) /= '!'.and.firstchar(record) > ' ') then
        if(status /= 0) cycle
        nclass = nclass + 1
      end if
    end do
    close(10)
  end if
end do
close(99)  

write(*,*) ' Number of atoms of the system: ', ndim
write(*,*) ' Number of classes in parameter files: ', nclass
if(nclass > ndim) ndim = nclass

return
end           
                           
!
! subroutine readpsf: Reads a PSF file to get charges and assign
!                     Lennard-Jones parameters that are provided 
!                     as input (which were read previously using
!                     subroutine readpar). In this program the
!                     only used parameter are the masses, actually.
!
! On input: psffile: name of the psf file
!           nclass: Number of atom classes read before by readpar
!           class: array containing the classes of each atom read
!                  by readpar      
!           eps: array containing the epsilon parameter of each type
!                of atom in each residue
!           sig: array containing the sigma parameter of each type
!                of atom in each residue
!           assignpars: .true. if the parameter file was read before
!                      and the parameters are needed, .false. 
!                      otherwise.
!
! On return: natom: total number of atoms of the PSF file
!            segat: name of the segment of each atom
!            resat: name of the residue of each atom 
!            resid: index of the residue of the atom
!            classat: class of each atom
!            typeat: type of each atom
!            q: array containing the charge of each atom
!            e: array containing the epsilon parameter for each atom
!            s: array containing the sigma parameter for each atom
!            mass: array containing the mass of each atom
!

subroutine readpsf(psffile,nclass,class,eps,sig,natom,segat,resat,&
                   resid,classat,typeat,q,e,s,mass,assignpars)   

use charsize
implicit none
integer :: i, j, natom, length, nclass, resid(*), status, iresid,&
           irlast
double precision :: eps(*), sig(*), q(*), e(*), s(*), mass(*)
character(len=charsize1) :: segat(*), resat(*), classat(*),& 
                            typeat(*), class(*)
character(len=200) :: psffile, record
logical :: found, assignpars

! Reading PSF file

open(10,file=psffile,action='read')
read(10,"( a200 )") record
do while(record(length(record)-5:length(record)).ne.'!NATOM')
  read(10,"( a200 )") record
end do
read(record,*) natom

resid(1) = 1
do i = 1, natom
  read(10,"( a200 )") record
  read(record,*,iostat=status) j, segat(i), iresid, resat(i),&
                               typeat(i), classat(i), q(i), mass(i)
  if ( status /= 0 ) then
    write(*,*) ' ERROR: Reading atom line in psf file: '
    write(*,*) record(1:length(record))
    write(*,*) ' Expected: number, segment, residue number, type,',&
               ' class, charge, mass '
    stop
  end if
  if(i.gt.1) then
    if(segat(i) == segat(i-1) .and. &
       resat(i) == resat(i-1) .and. &
       iresid == irlast ) then
      resid(i) = resid(i-1)
    else
      resid(i) = resid(i-1) + 1
    end if   
  end if
  irlast = iresid
end do
close(10)

! Assigning paramaters for each atom

if(.not.assignpars) return

do i = 1, natom
  j = 0                             
  found = .false.
  do while(.not.found.and.j < nclass)
    j = j + 1
    if(class(j) == classat(i)) then
      e(i) = eps(j)
      s(i) = sig(j)
      found = .true.
    end if
  end do   
  if(.not.found) then
    write(*,*) ' ERROR: Could not find Lennard-Jones parameters for atom:',&
                 segat(i),' ',resat(i),' ',typeat(i),' ',classat(i)
    stop
  end if
end do

return
end      


!
! Subroutine that checks if dcd file contains periodic cell information
!

subroutine chkperiod(dcdfile,dcdaxis,readfromdcd) 

implicit none
real :: dummyr, x
double precision :: side(6)
integer :: dummyi, ntotat, i, status
character(len=200) :: dcdfile
character(len=4) :: dummyc
logical :: dcdaxis, readfromdcd

dcdaxis = .false.
open(10,file=dcdfile,action='read',form='unformatted',status='old',iostat=status)
if( status /= 0 ) then
  write(*,*) ' ERROR: Error opening dcd file: ' 
  write(*,*) trim(dcdfile)
  stop
end if
read(10,iostat=status) dummyc, dummyi, (dummyi,i=1,8), dummyr,&
                       (dummyi,i=1,9)
if(status /= 0) then
  write(*,*) ' ERROR: Error reading dcd file: '
  write(*,*) trim(dcdfile)
  stop
end if
read(10,iostat=status) dummyi, dummyr
if(status /= 0) then
  write(*,*) ' ERROR: Error reading dcd file: '
  write(*,*) trim(dcdfile)
  stop
end if
read(10,iostat=status) ntotat
if(status /= 0) then
  write(*,*) ' ERROR: Error reading dcd file: '
  write(*,*) trim(dcdfile)
  stop
end if
read(10,iostat=status) (x,i=1,ntotat)
close(10)

! If periodic cell information was found:

if(status /= 0) then
  open(10,file=dcdfile,action='read',form='unformatted')
  read(10) dummyc, dummyi, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
  read(10) dummyi, dummyr
  read(10) ntotat
  read(10,iostat=status) (side(i),i=1,6)
  close(10)
  if( status /= 0 ) then
    write(*,"(' ERROR: Could not read either coordinates nor ',/,&
              '        periodic cell sizes in first line of ',/,&
              '        first frame of DCD file.' )")
    stop
  end if
  dcdaxis = .true.
  write(*,*) ' DCD file appears to contain periodic cell information. '
  write(*,"( '  Sides in first frame: ',3(f8.2) )") side(1), side(3), side(6)
  if(.not.readfromdcd) then
    write(*,*) ' Warning: Calculation will not use this information! '
    write(*,*) ' To use it set: periodic readfromdcd'
  end if
  return 

! If periodic cell information was not found:

else
  dcdaxis = .false.
  write(*,*) ' DCD file does not contain periodic cell information. '
  if(readfromdcd) then
    write(*,*) ' ERROR: Periodic cell information not found in dcd&
                 & file and periodic=readfromdcd '
    stop
  end if
end if

return
end

!
! Subroutine image: computes the minimum image given the coordinates
!                   and the size of the box. Modifies the coordinates
!                   given
!

subroutine image(x,y,z,sidex,sidey,sidez)

implicit none
real :: x, y, z
double precision :: sidex, sidey, sidez

if(sidex < 1.d-5 .or. sidey < 1.d-5 .or. sidez < 1.d-5) then
  write(*,*) ' ERROR: Found too short box side. '
  stop
end if
x = mod(x,sngl(sidex))
y = mod(y,sngl(sidey))
z = mod(z,sngl(sidez))
if(x > sidex/2) x = x - sidex
if(y > sidey/2) y = y - sidey
if(z > sidez/2) z = z - sidez
if(x < -sidex/2) x = x + sidex
if(y < -sidey/2) y = y + sidey
if(z < -sidez/2) z = z + sidez

return
end


!
! Subroutine readpar: Reads the parameter files in charmm format
!                     and assign for each class of atom its Lennard-Jones
!                     parameters
!
! On input: parfile: name of the parameter file to be read
!           nclass: number of classes of atom read including previous
!                   parameter files
!           class: array containing the class of each atom of each
!                  residue up to nclass
!
! On return: eps: array containing the epsilon parameter class of
!                 atom
!            sig: array containing the sigma parameter class of atom
!            nclass: number of classes of atoms including the current
!                    file
!            class: array containing the atom classes including
!                   the current classes 
!

subroutine readpar(parfile,nclass,class,eps,sig)

use charsize
implicit none
integer :: nclass, status
double precision :: eps(*), sig(*)
double precision :: dummyd
character :: firstchar
character(len=charsize1) :: class(*)
character(len=200) :: record, parfile

! Reading the parameter file

record(1:9) = '#########'
open(10,file=parfile,action='read')
do while (record(1:9) /= 'NONBONDED')
  read(10,"( a200 )") record   
end do
do while(.true.)
  read(10,"( a200 )",iostat=status) record
  if(status /= 0) exit
  if(firstchar(record) /= '!'.and.firstchar(record) > ' ') then
    read(record,*,iostat=status) class(nclass+1), dummyd,&
                                 eps(nclass+1), sig(nclass+1)
    if(status /= 0) cycle
    nclass = nclass + 1
  end if
end do
close(10)

return
end      

! Get the first non-blanck character from a line

function firstchar(record)

implicit none
character(len=200) :: record
character :: firstchar
integer :: i

i = 1
do while(record(i:i) <= ' ')
  i = i + 1
end do
firstchar = record(i:i)

return
end   

!
! Function that sets the length of a string
!

function length(string)

implicit none
integer :: length
character(len=200) :: string

length = 200
do while(string(length:length) <= ' ')
  length = length - 1
end do

return
end      
 
!
! Compute the norm of a vector
!

function xnorm(x1,x2,x3)

implicit none
real :: xnorm, x1, x2, x3

xnorm = sqrt(x1*x1 + x2*x2 + x3*x3)

return 
end

!
! Given two vectors returns the cosine of the
! angle between them 
!

function cosang(x1,y1,z1,x2,y2,z2)

implicit none
real :: cosang, x1, y1, z1, x2, y2, z2, xnorm

cosang = x1*x2 + y1*y2 + z1*z2
cosang = cosang / ( xnorm(x1,y1,z1) * xnorm(x2,y2,z2) )

return
end

! 
! Subroutine that computes the center of mass of a group
! of atoms from a sequential vector of coordinates.
! 
! On input:
!
! n_atoms: Number of atoms of the group.
! i_atoms: Vector with the atom indices in the structure.
! mass: Vector with atom masses of all atoms in the structure.
! mass_group: Total mass of the group.
! x_at, y_at, z_at: Current dcd coordinates.
! 
! On return:
!
! cmx, cmy, cmz: The center of mass of the group.
!

subroutine compute_cm(n_atoms,i_atoms,mass,mass_group,&
                      x_at,y_at,z_at,cmx,cmy,cmz)

implicit none
integer :: i, i_atoms(*), n_atoms
real :: x_at(*), y_at(*), z_at(*)
double precision :: mass(*), mass_group, cmx, cmy, cmz

cmx = 0.d0
cmy = 0.d0
cmz = 0.d0
do i = 1, n_atoms
  cmx = cmx + x_at(i)*mass(i_atoms(i))
  cmy = cmy + y_at(i)*mass(i_atoms(i))
  cmz = cmz + z_at(i)*mass(i_atoms(i))
end do
cmx = cmx / mass_group
cmy = cmy / mass_group
cmz = cmz / mass_group 

return
end 

! 
! Subroutine that computes the center of mass of a group
! of atoms from a frame of a dcd file.
! 
! On input:
!
! n_atoms: Number of atoms of the group.
! i_atoms: Vector with the atom indices in the structure.
! mass: Vector with atom masses of all atoms in the structure.
! mass_group: Total mass of the group.
! xdcd, ydcd, zdcd: Current dcd coordinates.
! iatom: Index of the first atom on this frame - 1
! 
! On return:
!
! cmx, cmy, cmz: The center of mass of the group.
!

subroutine compute_cm_dcd(n_atoms,i_atoms,mass,mass_group,&
                          xdcd,ydcd,zdcd,iatom,&
                          cmx,cmy,cmz)

implicit none
integer :: i, ii, n_atoms, i_atoms(*), iatom
double precision :: mass(*), mass_group, cmx, cmy, cmz
real :: xdcd(*), ydcd(*), zdcd(*)

cmx = 0.d0
cmy = 0.d0
cmz = 0.d0
do i = 1, n_atoms
  ii = iatom + i_atoms(i)
  cmx = cmx + xdcd(ii)*mass(i_atoms(i))
  cmy = cmy + ydcd(ii)*mass(i_atoms(i))
  cmz = cmz + zdcd(ii)*mass(i_atoms(i))
end do
cmx = cmx / mass_group
cmy = cmy / mass_group
cmz = cmz / mass_group 

return
end 



















 
                     
