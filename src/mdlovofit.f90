!
! mdlovofit: Computes the mobility of a PDB trajectory file
!            using a low-order-value-optimization strategy
!            that selects the best aligned subset as a reference.
!
! L. Martinez, Sep 16, 2013
! Institute of Chemistry, State University of Campinas - UNICAMP
!
! Home-Page: http://leandro.iqm.unicamp.br/mdlovofit
!
! Run with: mdlovofit -f [fraction] -t align.pdb file1.pdb file2.pdb file3.pdb ...
!
! where:
!
! [fraction] is (real number) the fraction of the atoms that will be considered
!            explicitly on the fit (that is, that will automatically be chosen
!            by the method as the best aligned atoms). Use any real number
!            between 0 and 1. For example: [fraction] = 0.7 for 70% of the atoms.
!
! align.pdb : This is the PDB file which will containt the aligned trajectory.
!
! file1.pdb etc. : These are the PDB files that contain the trajectory 
!                  to the aligned. Each PDB file may contain more than one
!                  frame of the trajectory, and the files will be considered
!                  as a sequential trajectory in the input order.
!
! Optional parameters are available. Check the Home-Page.
!

program mdlovofit
 
  implicit none
  
  integer :: narg, i, npdb, ioerror, ipdb, nall, iframe, iall, ic, mflash,&
             n_consider, n_not_consider, itrial, ntrial, seed, n_consider_safe,&
             itype, ntype, j, nc, iref
  integer, allocatable :: ialign(:), bije_new(:), bije_old(:), indflash(:),&
                          lflash(:), best_bije(:), bije(:), resnum(:)
  
  double precision :: frac, xread, yread, zread, xref_cm, yref_cm, zref_cm,&
                      xcm, ycm, zcm, u(3,3), random, rmsd_low, rmsd_high,&
                      rmsd_all, occup, xtemp, ytemp, ztemp, beta, &
                      av_rmsd_high, av_rmsd_low, av_rmsd_all, mapstep, &
                      best_rmsd_low, best_rmsd_high, best_rmsd_all, &
                      scores_curr, scores_prev, delta_scores, av_rmsd_low_last
  double precision, allocatable :: x(:), y(:), z(:), xref(:), yref(:), zref(:),&
                                   xall(:), yall(:), zall(:), xm(:), ym(:), zm(:),&
                                   xp(:), yp(:), zp(:), scores(:),&
                                   x_consider(:), y_consider(:), z_consider(:),&
                                   xref_consider(:), yref_consider(:), zref_consider(:),&
                                   x_average(:), y_average(:), z_average(:),&
                                   av_best_bije(:), rmsf(:)
  
  character(len=4) :: attype
  character(len=200) :: record, trajfile, rmsffile, atomsfile
  character(len=200), allocatable :: pdbfile(:) 
  character(len=30), allocatable :: pdbstring_left(:), pdbstring_right(:)
  character(len=4), allocatable :: atomtypes(:)
  
  logical :: printperframe, mapfrac, first_bije, print_rmsf, fileexists, endline
  logical, allocatable :: consider(:)
  
  write(*,"(a)") "#"
  write(*,"(a)") "# MDLovoFit - Version 20.0.0"
  write(*,"(a)") "#"
  
  ! Reading input parameters
  
  narg = iargc()
  fileexists = .false.
  frac = 0.d0
  trajfile = "NONE"
  printperframe = .false.
  print_rmsf = .false.
  mapfrac = .false.
  mapstep = 0.01d0
  npdb = 0
  ntrial = 100
  iref = 1
  itype = 0
  ntype = 0
  i = 1
  do while( i <= narg )
    call getarg(i,record)
    select case ( record )
      case ( "-f" ) 
        call getarg(i+1,record)
        read(record,*,iostat=ioerror) frac
        if ( ioerror /= 0 ) then 
          write(*,*) ' ERROR reading fraction. '
          call arg_error
        end if
        i = i + 2
      case ( "-iref" ) 
        call getarg(i+1,record)
        read(record,*,iostat=ioerror) iref
        if ( ioerror /= 0 ) then 
          write(*,*) ' ERROR reading number of reference frame (-iref [Int]). '
          call arg_error
        end if
        i = i + 2
      case ( "-t" ) 
        call getarg(i+1,record)
        read(record,*,iostat=ioerror) trajfile
        if ( ioerror /= 0 ) then 
          write(*,*) ' ERROR reading output trajectory file name. '
          call arg_error
        end if
        i = i + 2
      case ( "-perframe" , "-pf" ) 
        printperframe = .true.
        i = i + 1 
      case ( "-ntrial" ) 
        call getarg(i+1,record)
        read(record,*,iostat=ioerror) ntrial
        if ( ioerror /= 0 ) then 
          write(*,*) ' ERROR reading ntrial. '
          call arg_error
        end if
        i = i + 2
      case ( "-mapfrac" )
        mapfrac = .true.
        i = i + 1
      case ( "-mapstep" ) 
        call getarg(i+1,record)
        read(record,*,iostat=ioerror) mapstep
        if ( ioerror /= 0 ) then 
          write(*,*) ' ERROR reading mapstep. '
          call arg_error
        end if
        i = i + 2
      case ( "-rmsf" ) 
        call getarg(i+1,record)
        read(record,*,iostat=ioerror) rmsffile
        if ( ioerror /= 0 ) then 
          write(*,*) ' ERROR reading output RMSF file name. '
          call arg_error
        end if
        print_rmsf = .true.
        i = i + 2
      case ( "-atoms" )
        itype = 0
        i = i + 1
        call getarg(i,record)
        do while( i <= narg .and. record(1:1) /= '-' )
          call getarg(i,record)
          itype = itype + 1
          i = i + 1
        end do
        ntype = itype - 1
        i = i - 1
      case ( "-atomsfile" )
        call getarg(i+1,record)
        inquire(file=record,exist=fileexists)
        if ( fileexists ) then
          atomsfile = record
        else
           write(*,*) " ERROR: Could not find atom type file: ", trim(adjustl(record))
           stop
        end if
        i = i + 2
      case default
        npdb = npdb + 1
        i = i + 1
    end select
  end do

  if ( .not. mapfrac ) then
    if ( frac == 0.d0 .or. &
         trajfile == "NONE" .or. &
         npdb == 0 ) then
      write(*,*) ' ERROR: Could not find fraction (-f), output trajectory (-t), or input PDB file information. '
      call arg_error
      stop
    end if
  end if
  if ( .not. mapfrac ) then
    if ( npdb == 0 ) then 
      write(*,*) ' ERROR: Could not find input pdb file. '
      call arg_error
    end if
  end if
  
  if ( .not. mapfrac ) then 
    write(*,"(a,f5.2)") "# Fraction of best aligned atoms to consider: ", frac
    write(*,"( a, a )") "# Output trajectory file: ", trim(trajfile)
    write(*,"( a, L )") "# Different assignment for each frame: ", printperframe
    if ( print_rmsf ) write(*,"( a, a )") "# RMSF output file: ", trim(rmsffile)
  else
    write(*,"( a )") "# Scanning fractions."
  end if
  write(*,"( a, i5 )") "# Number of trials for each alignment: ", ntrial
  
  ! Types of atoms (by default, use CA only)

  if ( .not. fileexists .and. ntype == 0 ) then
    ntype = 1
    allocate( atomtypes(1) )
    atomtypes(1) = "CA"
  else
    if ( fileexists ) then
      write(*,"( a, a )") "# Atom selection file: ", trim(adjustl(atomsfile))
    end if
    if ( ntype > 0 ) then
      allocate( atomtypes(ntype) )
      i = 1
      call getarg(i,record)
      do while( record /= '-atoms' )
        i = i + 1
        call getarg(i,record)
      end do
      do itype = 1, ntype
        i = i + 1
        call getarg(i,record)
        atomtypes(itype) = record
      end do
      write(*,"( a )") "# User defined selection of atom types: "
      do i = 1, ntype
        write(*,"( a, tr1, i4, tr2, a )") "#", i, atomtypes(i)
      end do
    end if
  end if

  ! Read input PDB files
  
  write(*,"(a,/,a)") "#","# Trajectory files: "
  i = 1
  allocate(pdbfile(npdb))
  npdb = 0
  do while( i <= narg )
    call getarg(i,record)
    select case ( record )
      case ( "-f" , &
             "-t" , & 
             "-mapstep", &
             "-ntrial", &
             "-rmsf", &
             "-iref", &
             "-atomsfile" ) 
        i = i + 2
      case ( "-perframe" ,&
             "-pf" ,&
             "-mapfrac" ) 
        i = i + 1 
      case ( "-atoms" )
        i = i + ntype + 1
      case default
        npdb = npdb + 1
        call getarg(i,record)
        pdbfile(npdb) = trim(record)
        write(*,"(a,a)") "#   ",trim(adjustl(pdbfile(npdb)))
        i = i + 1
    end select
  end do
  write(*,"( a )") "#"
  
  ! Reading first frame of the first trajectory to determine the number of atoms
  open(10,file=pdbfile(1),status="old",action="read",iostat=ioerror)
  if ( ioerror /= 0 ) then
    write(*,"( a,a )") " ERROR: Could not open PDB file: ", trim(adjustl(pdbfile(1)))
    write(*,"( a )") " Note: The PDB input files must be the last parameters of the command line."
    stop
  end if
  nall = 0
  do 
    read(10,"( a200 )",iostat=ioerror) record
    if ( ioerror /= 0 ) exit
    if ( record(1:3) == "END" .or. &
         record(1:6) == "ENDMDL" ) exit
    if ( record(1:4) == "ATOM" .or. &
        record(1:6) == "HETATM" ) then
      read(record(31:54),"( 3(f8.3) )",iostat=ioerror) xread, yread, zread
      if ( ioerror /= 0 ) then
        write(*,"( a,a )") " ERROR: Problem reading coordinates in file: ", trim(adjustl(pdbfile(1)))
        write(*,"(a,i8)") "        Line: ", nall
        stop
      end if
      nall = nall + 1
    end if
  end do
  write(*,"( a, i8 )") "# Total number of atoms: ", nall
  allocate( consider(nall) )
  rewind(10)

  !
  ! Now determining which atoms will be considered for the alignment
  !
  do i = 1, nall
    consider(i) = .false.
  end do
 
  ! From the user-defined command-line-argument selection 
  nc = 0
  i = 0
  do
    read(10,"( a200 )",iostat=ioerror) record
    if ( ioerror /= 0 ) exit
    if ( record(1:3) == "END" .or. &
         record(1:6) == "ENDMDL" ) exit
    if ( record(1:4) == "ATOM" .or. &
        record(1:6) == "HETATM" ) then
      read(record(31:54),"( 3(f8.3) )",iostat=ioerror) xread, yread, zread
      i = i + 1
      read(record(13:16),*) attype
      do j = 1, ntype
        if ( attype == atomtypes(j) ) then
          consider(i) = .true.
          nc = nc + 1
          exit
        end if
      end do
    end if
  end do
  close(10)

  ! From the atom selection file
  if ( fileexists ) then
    open(10,file=atomsfile,status="old",iostat=ioerror)
    if ( ioerror /= 0 ) then
      write(*,*) " ERROR: Could not open atom selection file. "
      stop    
    end if
    i = 0
    do
      read(10,"( a200 )",iostat=ioerror) record
      if ( ioerror /= 0 ) exit
      if ( record(1:3) == "END" .or. &
           record(1:6) == "ENDMDL" ) exit
      if ( record(1:4) == "ATOM" .or. &
          record(1:6) == "HETATM" ) then
        read(record(31:54),"( 3(f8.3) )",iostat=ioerror) xread, yread, zread
        i = i + 1
        read(record(55:60),*) xread
        if ( xread == 1.d0 ) then
          nc = nc + 1
          consider(i) = .true.
        end if
      end if
    end do
    close(10)
  end if
  write(*,"( a, i8 )") "# Number of atoms of selection: ", nc
  if ( nc < 4 ) then
    write(*,"( a )") " ERROR: Atom selection has less than four atoms. "
    stop
  end if

  ! Allocate arrays

  allocate( x(nc), y(nc), z(nc), xref(nc), yref(nc), zref(nc),&
            xall(nall), yall(nall), zall(nall), pdbstring_left(nall), pdbstring_right(nall),&
            ialign(nc), xm(nc), ym(nc), zm(nc), xp(nc), yp(nc), zp(nc),&
            bije_new(nc), bije_old(nc), scores(nc),&
            indflash(nc), lflash(nc), rmsf(nc), resnum(nc), &
            x_consider(nc), y_consider(nc), z_consider(nc),&
            xref_consider(nc), yref_consider(nc), zref_consider(nc),&
            x_average(nc), y_average(nc), z_average(nc),&
            best_bije(nc), bije(nc), av_best_bije(nc) )
  mflash = 1 + nc / 10
  
  ! Define auxiliary vector ialign
  
  do i = 1, nc
    ialign(i) = i
    best_bije(i) = 0.d0
    av_best_bije(i) = 0.d0
  end do
  write(*,"(a,i8,/,a)") "# Total number of atoms found: ", nall, "#"
  
  ! Reading the reference frame, now to save the reference coordinates,
  ! compute the baricenter of the reference CA coordinates. Also using
  ! this loop to read the strings of the PDB file for later printing
  ipdb = 0
  ioerror = 1
  i = 1
  endline = .false.
  do while(i <= iref)
    ipdb = ipdb + 1
    if(ipdb > size(pdbfile) .or. i > iref) then
      write(*,*) "ERROR: -iref pointed to a frame that is not available."
      stop
    end if
    open(10,file=pdbfile(ipdb),status="old",action="read",iostat=ioerror)
    if ( ioerror /= 0 ) then
      write(*,"( a,a )") " ERROR: Could not open PDB file: ", trim(adjustl(pdbfile(i)))
      write(*,"( a )") " Note: The PDB input files must be the last parameters of the command line."
      stop
    end if
    ic = 0
    iall = 0
    do while(i <= iref) 
      read(10,"( a200 )",iostat=ioerror) record
      if ( ioerror /= 0 ) then
        if ( .not. endline ) i = i + 1
        close(10)
        exit
      end if
      if ( record(1:3) == "END" ) then
        ic = 0
        iall = 0
        i = i + 1
        endline = .true.
        cycle
      else
        endline = .false.
      end if
      if ( record(1:4) == "ATOM" .or. &
          record(1:6) == "HETATM" ) then
        iall = iall + 1
        pdbstring_left(iall) = record(1:30)
        pdbstring_right(iall) = record(67:80)
      else
        cycle
      end if
      if ( consider(iall) ) then 
        read(record(31:54),"( 3(f8.3) )",iostat=ioerror) xread, yread, zread
        if ( ioerror /= 0 ) then
          write(*,*) ' ERROR reading coordinates of atom: ', iall
          stop
        end if
        ic = ic + 1
        xref(ic) = xread
        yref(ic) = yread
        zref(ic) = zread
        if ( print_rmsf ) then
          read(record(23:26),*,iostat=ioerror) resnum(ic)
          if ( ioerror /= 0 ) then
            write(*,*) ' ERROR reading residue number of atom: ', iall
            stop
          end if
        end if
      end if
    end do
  end do
  
  ! Open output files
  
  if ( .not. mapfrac ) open(20,file=trajfile)
  
  if ( mapfrac ) then
    write(*,"('#  FRACTION   RMSD: BEST ALIGNED          OTHER ATOMS            ALL ATOMS          dBEST/dSTEP')")  
  else
    write(*,"('#     FRAME   RMSD: BEST ALIGNED          OTHER ATOMS            ALL ATOMS')")
  end if
  
  ! Reading the trajectory and aligning
  
  if ( mapfrac ) then
    frac = 0.d0
    av_rmsd_low_last = 0.d0
  end if
  map_fractions : do while( frac <= 2.d0 )
  
    if ( mapfrac ) frac = frac + mapstep
    n_consider = min(nc,int(frac*nc))
    if ( mapfrac ) then
      if ( n_consider < 4 ) cycle map_fractions
    end if
    n_not_consider = nc - n_consider 
  
    iframe = 0
    av_rmsd_low = 0.d0
    av_rmsd_high = 0.d0
    av_rmsd_all = 0.d0
    if ( print_rmsf ) then
      do i = 1, nc
        rmsf(i) = 0.d0
      end do
    end if
    do ipdb = 1, npdb
      open(10,file=pdbfile(ipdb),status="old",action="read",iostat=ioerror)
      if ( ioerror /= 0 ) then
        write(*,"( a,a )") " ERROR: Could not open PDB file: ", trim(adjustl(pdbfile(ipdb)))
        write(*,"( a )") " Note: The PDB input files must be the last parameters of the command line."
        stop
      end if
    
      ic = 0
      iall = 0
      do 
        read(10,"( a200 )",iostat=ioerror) record
        if ( ioerror /= 0 ) exit
    
        ! If this frame has ended, cycle to next frame
    
        if ( record(1:3) == "END" ) then
          ic = 0
          iall = 0
          cycle
        end if
  
        ! Read all coordinates of this frame (for later rotation/translation)
    
        if ( record(1:4) == "ATOM" .or. &
             record(1:6) == "HETATM" ) then
           if ( iall == 0 ) iframe = iframe + 1
           read(record(31:54),"( 3(f8.3) )",iostat=ioerror) xread, yread, zread
           if ( ioerror /= 0 ) then
             write(*,"( a,a )") " ERROR: Problem reading coordinates in file: ", trim(adjustl(pdbfile(ipdb)))
             stop
           end if
           iall = iall + 1
           if ( iall > nall ) cycle 
           xall(iall) = xread
           yall(iall) = yread
           zall(iall) = zread

        ! If this is not a atom line, just cycle

        else
          cycle
        end if
    
        ! Read CA coordinates of this frame
    
        if ( consider(iall) ) then
          read(record(31:54),"( 3(f8.3) )",iostat=ioerror) xread, yread, zread
          ic = ic + 1
          if ( ic == nc + 1 ) then
            write(*,*) ' Warning: frames have different number of atoms to consider. '
          end if
          if ( ic > nc ) cycle 
          x(ic) = xread
          y(ic) = yread
          z(ic) = zread
        end if
    
        ! If the frame is completely read, align 
    
        if ( iall == nall ) then
    
          seed = 12345
          best_rmsd_low = 1.d30
          trial: do itrial = 1, ntrial 
            
            ! Starting random bijections as first guess
    
            do i = 1, nc
              bije_old(i) = i
              scores(i) = random(seed)
            end do
            call flashsort(scores, nc, lflash, mflash, indflash)
            do i = 1, nc
              bije_new(i) = indflash(i)
            end do
    
            first_bije = .true.
            delta_scores = 1.d30
            scores_curr = 1.d30
            do while( delta_scores > 1.d-8 ) 
  
              ! If this is the first bijection (random one) to try in this trial, consider only
              ! the first 4 atoms in the alignment. This will increase the chances that
              ! for a large protein a substructure is used. If we used a large bijection, it would
              ! almost certainly have the same center of mass as the whole protein
  
              if ( first_bije ) then
                n_consider_safe = n_consider
                n_consider = 4
              end if
    
              ! Compute center of mass of considered atoms in this cycle
    
              xcm = 0.d0
              ycm = 0.d0
              zcm = 0.d0
              xref_cm = 0.d0
              yref_cm = 0.d0
              zref_cm = 0.d0
              do i = 1, n_consider
                xcm = xcm + x(bije_new(i))
                ycm = ycm + y(bije_new(i))
                zcm = zcm + z(bije_new(i))
                xref_cm = xref_cm + xref(bije_new(i))
                yref_cm = yref_cm + yref(bije_new(i))
                zref_cm = zref_cm + zref(bije_new(i))
              end do
              xcm = xcm / n_consider
              ycm = ycm / n_consider
              zcm = zcm / n_consider
              xref_cm = xref_cm / n_consider
              yref_cm = yref_cm / n_consider
              zref_cm = zref_cm / n_consider
              do i = 1, n_consider
                x_consider(i) = x(bije_new(i)) - xcm    
                y_consider(i) = y(bije_new(i)) - ycm    
                z_consider(i) = z(bije_new(i)) - zcm    
                xref_consider(i) = xref(bije_new(i)) - xref_cm
                yref_consider(i) = yref(bije_new(i)) - yref_cm
                zref_consider(i) = zref(bije_new(i)) - zref_cm
              end do
    
              ! Align using standard Procrustes method
    
              call align(n_consider,ialign,0,x_consider,y_consider,z_consider,&
                         xref_consider,yref_consider,zref_consider,&
                         u,&
                         xm,ym,zm,xp,yp,zp)
    
              ! Apply transformation to selected atoms to compute RMSD
    
              do i = 1, nc
                xtemp = x(i) - xcm
                ytemp = y(i) - ycm
                ztemp = z(i) - zcm
                call rotate(xtemp,ytemp,ztemp,u,xread,yread,zread)
                scores(i) = ( xread - ( xref(i) - xref_cm ) )**2 +&
                            ( yread - ( yref(i) - yref_cm ) )**2 +&
                            ( zread - ( zref(i) - zref_cm ) )**2
              end do
              call flashsort(scores, nc, lflash, mflash, indflash)
              do i = 1, nc
                bije_old(i) = bije_new(i)
                bije_new(i) = indflash(i)
              end do
  
              ! Restoring the correct number of atoms to be considered after first
              ! random alignment
  
              if ( first_bije ) then
                first_bije = .false.
                n_consider = n_consider_safe
              end if
  
              scores_prev = scores_curr
              scores_curr = 0.d0
              do i = 1, n_consider    
                scores_curr = scores_curr + scores(i)
              end do
              delta_scores = scores_curr - scores_prev
  
            end do

            ! Get inverse of the best assignment
    
            do i = 1, nc
              ic = bije_new(i)
              bije_old(ic) = i
            end do
    
            ! Compute RMSDs of final alignment of this frame
    
            rmsd_low = 0.d0
            do i = 1, n_consider
              rmsd_low = rmsd_low + scores(i)
            end do
            rmsd_low = dsqrt( rmsd_low / n_consider )

            ! Saving the best alignment obtained so far
  
            if ( rmsd_low < best_rmsd_low ) then
  
              best_rmsd_low = rmsd_low
              if ( n_not_consider > 0 ) then
                best_rmsd_high = 0.d0
                do i = n_consider + 1, nc
                  best_rmsd_high = best_rmsd_high + scores(i)
                end do
                best_rmsd_high = dsqrt( best_rmsd_high / n_not_consider )
                best_rmsd_all = (best_rmsd_low*n_consider + best_rmsd_high*n_not_consider)/nc
              else
                best_rmsd_high = 0.d0
                best_rmsd_all = best_rmsd_low
              end if
  
              do i = 1, n_consider
                av_best_bije(bije_new(i)) = av_best_bije(bije_new(i)) - 1.d0
              end do
              do i = 1, nc
                best_bije(i) = bije_new(i)
              end do
  
            end if
  
          end do trial
  
          ! Summing up the best rmsd found to compute average over frames
  
          av_rmsd_low = av_rmsd_low + best_rmsd_low
          av_rmsd_high = av_rmsd_high + best_rmsd_high
          av_rmsd_all = av_rmsd_all + best_rmsd_all
  
          if ( print_rmsf ) then
            do i = 1, nc
              rmsf(i) = rmsf(i) + scores(i)
            end do
          end if
  
          ! If this is a standard trajectory alignment (.not. mapfrac) and the user wants
          ! the RMSD per atom and the considered atoms to be printed at every frame, do it
      
          if ( .not. mapfrac .and. printperframe ) then
            write(*,"( i11, 3(tr2,f19.9) )") iframe, best_rmsd_low, best_rmsd_high, best_rmsd_all
    
            ! Apply transformation to all atoms to print file
  
            write(20,"( a,i8 )") "REMARK FRAME: ", iframe
            ic = 0
            do i = 1, nall
              xtemp = xall(i) - xcm 
              ytemp = yall(i) - ycm 
              ztemp = zall(i) - zcm 
              call rotate(xtemp,ytemp,ztemp,u,xread,yread,zread)
              xread = xread + xref_cm
              yread = yread + yref_cm
              zread = zread + zref_cm
              record = pdbstring_left(i)
              beta = 0.d0
              occup = 0.d0
              if ( consider(i) ) then
                ic = ic + 1
                if ( bije_old(ic) <= n_consider ) occup = 1.00
                beta = dmin1(dsqrt(scores(ic)),99.99d0)
              end if
              write(20,"( a30,3(f8.3),2(f6.2), a30 )") pdbstring_left(i),&
                                              xread, yread, zread,&
                                              occup, beta, pdbstring_right(i)
            end do
            write(20,"( a )") "END"
          end if
    
        end if
      end do
      close(10)
    end do
  
    if ( iframe > 1 ) then
      av_rmsd_low = av_rmsd_low / (iframe-1)
      av_rmsd_high = av_rmsd_high / (iframe-1)
      av_rmsd_all = av_rmsd_all /  (iframe-1)
    end if
  
    if ( mapfrac ) then
      write(*,"( f11.4, 4(tr2,f19.9) )") frac, av_rmsd_low, av_rmsd_high, av_rmsd_all,&
              (av_rmsd_low - av_rmsd_low_last)/mapstep
      av_rmsd_low_last = av_rmsd_low
      cycle map_fractions
    end if
  
    if ( printperframe ) then
      write(*,"(a,f12.5)") "# Average RMSD of least mobile atoms: ", av_rmsd_low
      write(*,"(a,f12.5)") "# Average RMSD of most mobile atoms:  ", av_rmsd_high
      write(*,"(a,f12.5)") "# Average RMSD of all atoms:          ", av_rmsd_all
    end if
  
    if ( print_rmsf ) then
      open(10,file=rmsffile)
      do i = 1, nc
        write(10,*) resnum(i), dsqrt(rmsf(i) / (iframe-1))
      end do
      close(10)
      write(*,"(a,a)") "# Wrote RMSF output file:             ", trim(rmsffile)
    end if
  
    ! If the user wanted the best alignment for each frame, nothing will be done from now on.
    ! However, probably it is easier to intepret the result using the average best bijection
    ! obtained along the trajectory. In this case, the alignment of each frame using this
    ! average, fixed, bijection to the reference is done, and the RMSD is computed again. 
           
    ! Redoing the alignment for the average best bijection
    
    if ( .not. printperframe ) then
    
      av_rmsd_low = 0.d0
      av_rmsd_high = 0.d0
      av_rmsd_all = 0.d0
      if ( print_rmsf ) then
        do i = 1, nc
          rmsf(i) = 0.d0
        end do
      end if
      
      call flashsort(av_best_bije, nc, lflash, mflash, bije_new)
      
      ! Saving inverse assignment in bije_old
      
      do i = 1, nc
        ic = bije_new(i)
        bije_old(ic) = i
      end do
      
      ! Computing center of mass of reference atoms
      
      xref_cm = 0.d0
      yref_cm = 0.d0
      zref_cm = 0.d0
      do i = 1, n_consider
        xref_cm = xref_cm + xref(bije_new(i))
        yref_cm = yref_cm + yref(bije_new(i))
        zref_cm = zref_cm + zref(bije_new(i))
      end do
      xref_cm = xref_cm / n_consider
      yref_cm = yref_cm / n_consider
      zref_cm = zref_cm / n_consider
      
      ! Move all CA atoms of the reference according to this center of mass
      
      do i = 1, nc
        xref_consider(i) = xref(bije_new(i)) - xref_cm
        yref_consider(i) = yref(bije_new(i)) - yref_cm
        zref_consider(i) = zref(bije_new(i)) - zref_cm
      end do
   
      iframe = 0
      do ipdb = 1, npdb
        open(10,file=pdbfile(ipdb),status="old",action="read",iostat=ioerror)
        if ( ioerror /= 0 ) then
          write(*,"( a )") " Error opening PDB file: ", trim(adjustl(pdbfile(ipdb)))
          stop
        end if
      
        ic = 0
        iall = 0
        do 
          read(10,"( a200 )",iostat=ioerror) record
          if ( ioerror /= 0 ) exit
      
          ! If this frame has ended, cycle to next frame
      
          if ( record(1:3) == "END" ) then
            ic = 0
            iall = 0
            cycle
          end if
  
          ! Read all coordinates of this frame (for later rotation/translation)
      
          if ( record(1:4) == "ATOM" .or. &
               record(1:6) == "HETATM" ) then
             if ( iall == 0 ) iframe = iframe + 1
             read(record(31:54),"( 3(f8.3) )",iostat=ioerror) xread, yread, zread
             if ( ioerror /= 0 ) then
               write(*,"( a,a )") " ERROR: Problem reading coordinates in file: ", trim(adjustl(pdbfile(ipdb)))
               stop
             end if
             iall = iall + 1
             if ( iall > nall ) cycle
             xall(iall) = xread
             yall(iall) = yread
             zall(iall) = zread

          ! If this is not a coordinate line, cycle
          else
            cycle
          end if
      
          ! Read CA coordinates of this frame
      
          if ( consider(iall) ) then
             read(record(31:54),"( 3(f8.3) )",iostat=ioerror) xread, yread, zread
             ic = ic + 1
             x(ic) = xread
             y(ic) = yread
             z(ic) = zread
          end if
      
          ! If the frame is completely read, align 
      
          if ( iall == nall ) then
      
            ! Compute center of mass of considered CA atoms in this frame
      
            xcm = 0.d0
            ycm = 0.d0
            zcm = 0.d0
            do i = 1, n_consider
              xcm = xcm + x(bije_new(i))
              ycm = ycm + y(bije_new(i))
              zcm = zcm + z(bije_new(i))
            end do
            xcm = xcm / n_consider
            ycm = ycm / n_consider
            zcm = zcm / n_consider
            do i = 1, n_consider
              x_consider(i) = x(bije_new(i)) - xcm    
              y_consider(i) = y(bije_new(i)) - ycm    
              z_consider(i) = z(bije_new(i)) - zcm    
            end do
      
            ! Align using standard Procrustes method
      
            call align(n_consider,ialign,0,x_consider,y_consider,z_consider,&
                       xref_consider,yref_consider,zref_consider,&
                       u,&
                       xm,ym,zm,xp,yp,zp)
      
            ! Compute scores per atom
      
            do i = 1, nc
              xtemp = x(i) - xcm
              ytemp = y(i) - ycm
              ztemp = z(i) - zcm
              call rotate(xtemp,ytemp,ztemp,u,xread,yread,zread)
              scores(i) = ( xread - xref_consider(bije_old(i)) )**2 + &
                          ( yread - yref_consider(bije_old(i)) )**2 + &
                          ( zread - zref_consider(bije_old(i)) )**2
            end do
      
            ! Compute RMSDs of final alignment of this frame
      
            rmsd_low = 0.d0
            rmsd_high = 0.d0
            do i = 1, n_consider
              rmsd_low = rmsd_low + scores(bije_new(i))
            end do
            if ( n_not_consider > 0 ) then 
              do i = n_consider + 1, nc
                rmsd_high = rmsd_high + scores(bije_new(i))
              end do
              rmsd_high = dsqrt( rmsd_high / n_not_consider )
            else
              rmsd_high = 0.d0
            end if
            rmsd_low = dsqrt( rmsd_low / n_consider )
            rmsd_all = (rmsd_low*n_consider + rmsd_high*n_not_consider)/nc
            
            av_rmsd_low = av_rmsd_low + rmsd_low
            av_rmsd_high = av_rmsd_high + rmsd_high
            av_rmsd_all = av_rmsd_all + rmsd_all
            write(*,"( i11, 3(tr2,f19.9) )") iframe, rmsd_low, rmsd_high, rmsd_all
  
            if ( print_rmsf ) then 
              do i = 1, nc
                rmsf(i) = rmsf(i) + scores(i)
              end do
            end if
        
            ! Apply transformation to all atoms and print aligned trajectory file
    
            write(20,"( a,i8 )") "REMARK FRAME: ", iframe
            ic = 0
            do i = 1, nall
              xtemp = xall(i) - xcm 
              ytemp = yall(i) - ycm 
              ztemp = zall(i) - zcm 
              call rotate(xtemp,ytemp,ztemp,u,xread,yread,zread)
              xread = xread + xref_cm
              yread = yread + yref_cm
              zread = zread + zref_cm
              record = pdbstring_left(i)
              occup = 0.d0
              beta = 0.d0
              if ( consider(i) ) then
                ic = ic + 1
                if ( bije_old(ic) <= n_consider ) occup = 1.d0
                beta = dmin1(dsqrt(scores(ic)),99.99d0)
              end if
              write(20,"( a30,3(f8.3),2(f6.2),a30 )") pdbstring_left(i),&
                                              xread, yread, zread,&
                                              occup, beta, pdbstring_right(i)
            end do
            write(20,"( a )") "END"
          end if 
    
        end do
        close(10)
      end do
      if ( iframe > 1 ) then
        av_rmsd_low = av_rmsd_low / (iframe-1)
        av_rmsd_high = av_rmsd_high / (iframe-1)
        av_rmsd_all = av_rmsd_all / (iframe-1)
      end if
      write(*,"(a,f12.5)") "# Average RMSD of least mobile atoms: ", av_rmsd_low
      write(*,"(a,f12.5)") "# Average RMSD of most mobile atoms:  ", av_rmsd_high
      write(*,"(a,f12.5)") "# Average RMSD of all atoms:          ", av_rmsd_all
      if ( iframe == 1 ) then
        write(*,"(a)") "# WARNING: only one frame found. Missing END/ENDMDL flags between frames?"
      end if
  
      if ( print_rmsf ) then
        open(10,file=rmsffile)
        do i = 1, nc
          write(10,*) resnum(i), dsqrt(rmsf(i) / (iframe-1))
        end do
        close(10)
        write(*,"(a,a)") "# Wrote RMSF output file:             ", trim(rmsffile)
      end if
  
    end if
  
    close(20)
    exit map_fractions
  
  end do map_fractions

end program mdlovofit

!
! Subroutine that prints instructions when there are input argument errors
!

subroutine arg_error

  write(*,"(a)") " "
  write(*,"(a)") " Run with: mdlovofit -f [fraction] -t align.pdb file1.pdb file2.pdb ..."
  write(*,"(a)") " "
  write(*,"(a)") " where:"
  write(*,"(a)") " "
  write(*,"(a)") " [fraction] is (real number) the fraction of the atoms that will be considered"
  write(*,"(a)") "            explicitly on the fit (that is, that will automatically be chosen"
  write(*,"(a)") "            by the method as the best aligned atoms). Use any real number"
  write(*,"(a)") "            between 0 and 1. For example: [fraction] = 0.7 for 70% of the atoms."
  write(*,"(a)") " "
  write(*,"(a)") " align.pdb : This is the PDB file which will containt the aligned trajectory."
  write(*,"(a)") " "
  write(*,"(a)") " file1.pdb etc. : These are the PDB files that contain the trajectory"
  write(*,"(a)") "                  to the aligned. Each PDB file may contain more than one"
  write(*,"(a)") "                  frame of the trajectory, and the files will be considered"
  write(*,"(a)") "                  as a sequential trajectory in the input order."
  write(*,"(a)") " "
  write(*,"(a)") " Optional parameters are available. Look for detailed instructions at:"
  write(*,"(a)") " "
  write(*,"(a)") "                 http://leandro.iqm.unicamp.br/mdlovofit"
  write(*,"(a)") " "
  write(*,"(a)") " "
  stop

end subroutine arg_error
             

!
! Subroutine align: given two sets of vectors, finds the best
!                   rotation that to align the two sets. Returns
!                   the second vector aligned to the first vector
!
!                 Method: S. K. Kearsley, 
!                         "On the orthogonal transformation used for
!                          structural comparisons"
!                         Acta Cryst. (1989) A45, 208-210 
!
!  Author: Leandro Martinez, IQ-UNICAMP, 26/10/2005
!
!  On input: n_align: the number of atoms of the group to be aligned
!            i_align: the indexes in vectors x-, y-, z- and mass of the atoms
!                     of the atoms to be aligned (1,...,n_align)
!            iatom: index of atom 0 (first atom-1) for xdcd, ydcd and zdcd
!            xdcd: The vector containing current x coordinates
!            ydcd: The vector containing current y coordinates
!            zdcd: The vector containing current z coordinates
!            xref: The vector containing referece x coordinates
!            yref: The vector containing referece y coordinates
!            zref: The vector containing referece z coordinates
!            Auxiliar arrays: xm, ym, zm, xp, yp, zp are vectors
!                             that have n_align positions. 
!
!  On return: u(3,3): rotation matrix to be applied
!
!  Obs: It is assumed that the atoms of both sets are already moved
!       in such a way that their center of mass is in the origin.
!

subroutine align(n_align,i_align,iatom,xdcd,ydcd,zdcd,&
                 xref,yref,zref,u,xm,ym,zm,xp,yp,zp)

  implicit none
  integer :: i, j, iq, n_align, i_align(*), iatom, ii
  double precision :: xdcd(*), ydcd(*), zdcd(*), xref(*), yref(*), zref(*)
  double precision :: u(3,3), a(4,4), q(4,4), qmin,&
                      xm(*), ym(*), zm(*), xp(*), yp(*), zp(*)
   
  ! Interface to jacobi subroutine, which is in Fortran77
  
  interface 
    subroutine jacobi(a,v,n)
    integer :: n
    double precision :: a(4,4), v(4,4)
    end subroutine jacobi
  end interface
  
  ! If the number of atoms of this group is 1, return the identity matrix
  
  if( n_align == 1 ) then
    do i = 1, 3
      do j = 1, 3
        if(i == j) then
          u(i,j) = 1.d0
        else
          u(i,j) = 0.d0
        end if
      end do
    end do
    return
  end if
   
  ! Computing the quaternion matrix
  
  do i = 1, n_align
    ii = iatom + i_align(i)
    xm(i) = xref(i) - xdcd(ii)
    ym(i) = yref(i) - ydcd(ii)
    zm(i) = zref(i) - zdcd(ii)
    xp(i) = xref(i) + xdcd(ii)
    yp(i) = yref(i) + ydcd(ii)
    zp(i) = zref(i) + zdcd(ii)
  end do
  
  do i = 1, 4
    do j = 1, 4
      q(i,j) = 0.d0
    end do
  end do
   
  do i = 1, n_align
    q(1,1) = q(1,1) + xm(i)**2 + ym(i)**2 + zm(i)**2
    q(1,2) = q(1,2) + yp(i)*zm(i) - ym(i)*zp(i)
    q(1,3) = q(1,3) + xm(i)*zp(i) - xp(i)*zm(i)
    q(1,4) = q(1,4) + xp(i)*ym(i) - xm(i)*yp(i)
    q(2,2) = q(2,2) + yp(i)**2 + zp(i)**2 + xm(i)**2
    q(2,3) = q(2,3) + xm(i)*ym(i) - xp(i)*yp(i)
    q(2,4) = q(2,4) + xm(i)*zm(i) - xp(i)*zp(i)
    q(3,3) = q(3,3) + xp(i)**2 + zp(i)**2 + ym(i)**2
    q(3,4) = q(3,4) + ym(i)*zm(i) - yp(i)*zp(i)
    q(4,4) = q(4,4) + xp(i)**2 + yp(i)**2 + zm(i)**2
  end do
  q(2,1) = q(1,2)
  q(3,1) = q(1,3)
  q(3,2) = q(2,3)
  q(4,1) = q(1,4)
  q(4,2) = q(2,4)
  q(4,3) = q(3,4)  
     
  ! Computing the eigenvectors 'a' and eigenvalues 'q' of the q matrix
  
  call jacobi(q,a,4)
   
  ! Choosing the quaternion that corresponds to the minimum
  
  iq = 1
  qmin = q(1,1)
  do i = 2, 4
    if(q(i,i).lt.qmin) then
      iq = i
      qmin = q(i,i)
    end if
  end do
  
  ! Computing the rotation matrix
  
  u(1,1) = a(1,iq)**2 + a(2,iq)**2 - a(3,iq)**2 - a(4,iq)**2
  u(1,2) = 2. * ( a(2,iq)*a(3,iq) + a(1,iq)*a(4,iq) )
  u(1,3) = 2. * ( a(2,iq)*a(4,iq) - a(1,iq)*a(3,iq) )  
  u(2,1) = 2. * ( a(2,iq)*a(3,iq) - a(1,iq)*a(4,iq) )  
  u(2,2) = a(1,iq)**2 + a(3,iq)**2 - a(2,iq)**2 - a(4,iq)**2 
  u(2,3) = 2. * ( a(3,iq)*a(4,iq) + a(1,iq)*a(2,iq) )  
  u(3,1) = 2. * ( a(2,iq)*a(4,iq) + a(1,iq)*a(3,iq) )  
  u(3,2) = 2. * ( a(3,iq)*a(4,iq) - a(1,iq)*a(2,iq) )  
  u(3,3) = a(1,iq)**2 + a(4,iq)**2 - a(2,iq)**2 - a(3,iq)**2 
  
  return

end subroutine align

!
! Subroutine rotate: Rotates the coordinates of a vector
!                    given the rotation matrix
! 
! On input:
!   x, y, z: Vector containing atoms' coordinates.
!   u: rotation matrix (3x3)
!
! On output:
!   xmove, ymove, zmove: Coordinates rotated and translated.
!

subroutine rotate(x,y,z,u,xmove,ymove,zmove)
  
  implicit none
  double precision :: x, y, z, u(3,3), xmove, ymove, zmove
  
  xmove = u(1,1)*x + u(1,2)*y + u(1,3)*z
  ymove = u(2,1)*x + u(2,2)*y + u(2,3)*z
  zmove = u(3,1)*x + u(3,2)*y + u(3,3)*z
  
  return

end subroutine rotate

!
! Random number generator
!

double precision function random(sem)
  implicit none
  integer sem, mult
  sem = mod( mult( sem, 3141581) + 1, 100000000)
  random = sem/100000000.0d0
  return
end function random

integer function mult( p, q)
  implicit none
  integer p, q, p0, p1, q0, q1
  p1 = p/10000
  p0 = mod(p,10000)
  q1 = q/10000
  q0 = mod(q,10000)
  mult = mod( mod( p0*q1+p1*q0,10000)*10000+p0*q0,100000000)
  return
end function mult
