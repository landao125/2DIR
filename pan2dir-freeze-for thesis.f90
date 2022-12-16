program evv2dir
  use kinds,  only : dp,i4b
  use constants
  use cal_x3
  use parameters,only : lowerw1,upperw1,stepw1,lowerw2,upperw2, &
                        stepw2, w1,w2,g16filename,inputfilename,&
                        lreadchk,g16chkname,lreadpes,g16pesname,&
                        gamma,scaling,gauss_revision_main,      &
                        gauss_revision_minor 
                       
  use io,only : open_file,close_file,io_file_unit,findkline,    &
                findkword,stdout,stdout_name
  use readinput,only : get_inputfile
  implicit none

  integer :: istat,itmp
  real (kind=dp) :: w1cc, w2cc,xhhh,xhhv,xhvh,evvhhh,evvhhv,evvxhvh
  integer (i4b) :: i,j,k,l,num_error,i_,j_,k_
  integer :: g16fileunit,g16chkunit,g16pesunit
  integer :: natoms,iatom,ipol,jpol,imode,jmode
  integer :: nqm,nqmf
  real(kind=dp),allocatable :: frc_consts(:),ir_inten(:)
  real(kind=dp),allocatable :: normal_coord(:,:,:),mass(:)                                                            
  real(kind=dp),allocatable :: car_coordinate(:,:,:,:),&
                               car_coordinate0(:,:)
  real(kind=dp),allocatable :: normal_coord_r(:),normal_coord_X(:,:)
  real(kind=dp),allocatable :: normal_e(:,:,:),normal_e1(:,:,:),&
                               normal_e2(:,:,:)
  real(kind=dp),allocatable :: normal_e_org(:,:,:), &
                               normal_e1_org(:,:,:),&
                               normal_e2_org(:,:,:)
  real(kind=dp),allocatable :: freqs_tmp(:),FI(:)
  integer,allocatable       :: mode_index(:)
  character(len=7),allocatable :: mode_char(:)
  character(len=maxlen) :: format_tmp,ctmp
  real(kind=dp) :: mu0(3),alpha0(3,3)
  real(kind=dp),allocatable :: mu(:,:,:),alpha(:,:,:,:)
  real(kind=dp),allocatable :: DDip0(:,:,:),DDip(:,:,:,:,:)
  integer :: nline,iline
  logical :: lout
  character(len=7) :: chari,charj,chark
  real(kind=dp) :: rtmp,tmp1,tmp2
  character(len=maxlen) :: specname="2dirspec.dat"
  integer :: specunit
  
  stdout = io_file_unit()
  open(unit=stdout,file=stdout_name,status='REPLACE')

  
  write (stdout,*)
  write (stdout,*),'*********************************************' 
  write (stdout,*),'( ___)( \/ )( \/ )  (__ \ (  _ \(_  _)(  _ \ '
  write (stdout,*),' )__)  \  /  \  /    / _/  )(_) )_)(_  )   / '
  write (stdout,*),'(____)  \/    \/    (____)(____/(____)(_)\_) '
  write (stdout,*),'                                             '
  write (stdout,*),'                   Fengqin Long    May 2022  '
  write (stdout,*),'*********************************************'
  write (stdout,*)  
  
  
   
  call get_inputfile(inputfilename)
   
  ! get natoms
  g16fileunit = io_file_unit()
  call open_file(g16filename,g16fileunit)
  call findkline(g16fileunit,"Revision",15,22)
  read(g16fileunit,"(10X,A2,11X,A4)") gauss_revision_main,&
                                      gauss_revision_minor
  write(stdout,*) "The Gaussian is calculate with Version:",&
                  gauss_revision_main,gauss_revision_minor
  
  ! get atom numbers
  call findkword(g16fileunit,"NAtoms=")
  !NAtoms=     78 NQM=       70 
  !NQMF=       8 NMMI=      0 NMMIF=      0
  read(g16fileunit,"(8X,I7,5X,I9,6X,I8)") natoms,nqm,nqmf
  if(nqmf==0) then
    nmode = 3*natoms - 6
  else
    nmode = 3*nqm
  endif

 
  g16pesunit = io_file_unit()
  call open_file(g16pesname,g16pesunit)
  read(g16pesunit,"(8X,I5)") itmp
  if(itmp /= natoms) then
    write(stdout,"(A)") &
      "Error!!! The number of atoms in *.log /= atoms in PES.out"
    write(stdout,"(A38,A30)") &
     "The PES.out is not the debuge file of ",&
     trim(adjustl(g16filename))
  endif
  

  allocate(freqs(nmode), redmass(nmode))
  allocate(frc_consts(nmode),ir_inten(nmode))
  allocate(normal_coord(3,natoms,nmode))
  freqs = 0.0
  redmass = 0.0
  normal_coord   = 0.0

  !read Dipole and Polarizability for equilibrium structure
  call findkline(g16fileunit," Dipole        =",1,16)
  read(g16fileunit,"(16X,3E15.8)")(mu0(ipol),ipol=1,3)
  read(g16fileunit,"(16X,3E15.8)")alpha0(1,1),alpha0(1,2),alpha0(2,2)
  read(g16fileunit,"(16X,3E15.8)")alpha0(1,3),alpha0(2,3),alpha0(3,3)    

  ! get Harmonic frequencies(cm**-1)
  call findkline(g16fileunit," Harmonic frequencies (cm**-1)",1,30)
  do i=1,4
    read(g16fileunit,*)
  enddo
  do imode=1,nmode
    if(MOD(imode,3)==1) then
      read(g16fileunit,*)
      read(g16fileunit,*)
      read(g16fileunit,"(15X,F12.4)") freqs(imode)
      read(g16fileunit,"(15X,F12.4)") redmass(imode)
      read(g16fileunit,"(15X,F12.4)") frc_consts(imode)
      read(g16fileunit,"(15X,F12.4)") ir_inten(imode)
      read(g16fileunit,*)
      do iatom=1,nqm !natoms
        read(g16fileunit,"(12X,3(F7.2))") &
            (normal_coord(ipol,iatom,imode),ipol=1,3)
      enddo
      do iatom=1,nqm+7!natoms+7
        backspace(g16fileunit)
      enddo
    elseif(MOD(imode,3)==2) then
      read(g16fileunit,*)
      read(g16fileunit,*)
      read(g16fileunit,"(38X,F12.4)") freqs(imode)
      read(g16fileunit,"(38X,F12.4)") redmass(imode)
      read(g16fileunit,"(38X,F12.4)") frc_consts(imode)
      read(g16fileunit,"(38X,F12.4)") ir_inten(imode)
      read(g16fileunit,*)
      do iatom=1,nqm !natoms
        read(g16fileunit,"(12X,23X,3(F7.2))") &
            (normal_coord(ipol,iatom,imode),ipol=1,3)
      enddo
      do iatom=1,nqm+7 !natoms+7
        backspace(g16fileunit)
      enddo      
    else
      read(g16fileunit,*)
      read(g16fileunit,*)
      read(g16fileunit,"(61X,F12.4)") freqs(imode)
      read(g16fileunit,"(61X,F12.4)") redmass(imode)
      read(g16fileunit,"(61X,F12.4)") frc_consts(imode)
      read(g16fileunit,"(61X,F12.4)") ir_inten(imode)
      read(g16fileunit,*)
      do iatom=1,nqm!natoms
        read(g16fileunit,"(58X,3(F7.2))") &
            (normal_coord(ipol,iatom,imode),ipol=1,3)
      enddo      
    endif
  enddo
          
  
  if(lreadchk) then
   
  ! read atom mass (amu.) in logfile
    call findkline(g16fileunit," - Thermochemistry -",1,20)
    allocate(mass(natoms))
    mass = 0.0
    
    read(g16fileunit,*)
    read(g16fileunit,*)
    read(g16fileunit,*)
    do iatom=1,natoms
      read(g16fileunit,"(41X,F10.5)") mass(iatom)
    enddo   
    
  ! read normal_coord in chkfile with high precision.
    g16chkunit = io_file_unit()
    call open_file(g16chkname,g16chkunit)
    
    allocate(normal_coord_r(nmode))
    allocate(normal_coord_X(nmode,nmode))
    normal_coord_r = 0.0
    normal_coord_X = 0.0      
    
    call findkword(g16chkunit,"Vib-Modes")
    read(g16chkunit,*)
    read(g16chkunit,*) (((normal_coord(ipol,iatom,imode),ipol=1,3),&
                       iatom=1,natoms),imode=1,nmode)
    !归一不正交，e''
      
    do imode=1,nmode
      normal_coord_r(imode) = sqrt(SUM(normal_coord(:,:,imode)**2))
      do jmode=1,nmode
        normal_coord_X(jmode,imode) = SUM(normal_coord(:,:,imode)*&
                                      normal_coord(:,:,jmode))
      enddo
    enddo
    
    allocate(normal_e(3,natoms,nmode),normal_e1(3,natoms,nmode))
    allocate(normal_e2(3,natoms,nmode),normal_e_org(3,natoms,nmode))
    allocate(normal_e1_org(3,natoms,nmode))
    allocate(normal_e2_org(3,natoms,nmode))
    ! change e2 to e and e1 
    ! e 正交归一,e1 既不正交也不归一，e2 归一不正交
    normal_e2 = normal_coord
    do imode=1,nmode
      do iatom=1,natoms
        normal_e(:,iatom,imode) = normal_e2(:,iatom,imode)*&
                                  sqrt(mass(iatom)/redmass(imode))
        normal_e1(:,iatom,imode)= normal_e2(:,iatom,imode)*&
                                  sqrt(mass_H/redmass(imode))
      enddo
    enddo    
     
    do imode=1,nmode
      normal_coord_r(imode) = sqrt(SUM(normal_e(:,:,imode)**2))
      do jmode=1,nmode
        normal_coord_X(jmode,imode) = SUM(normal_e(:,:,imode)*&
                                      normal_e(:,:,jmode))
      enddo
    enddo
    
    do imode=1,nmode
      normal_coord_r(imode) = sqrt(SUM(normal_e1(:,:,imode)**2))
      do jmode=1,nmode
        normal_coord_X(jmode,imode) = SUM(normal_e1(:,:,imode)*&
                                      normal_e1(:,:,jmode))
      enddo
    enddo    
    
    do imode=1,nmode
      normal_coord_r(imode) = sqrt(SUM(normal_e2(:,:,imode)**2))
      do jmode=1,nmode
        normal_coord_X(jmode,imode) = SUM(normal_e2(:,:,imode)*&
                                      normal_e2(:,:,jmode))
      enddo
    enddo   
    
    !read Cartesian coordinates in PES.out 
    allocate(car_coordinate0(3,natoms))
    allocate(car_coordinate(3,natoms,2,nmode))

    call findkline(g16pesunit," Cartesian Coordinates:",1,23)
    read(g16pesunit,*)
    do iatom=1,natoms
      read(g16pesunit,*) (car_coordinate0(ipol,iatom),ipol=1,3)
    enddo
    
    do imode=1,nmode
      do i=1,2
        call findkline(g16pesunit," Cartesian Coordinates:",1,23)
        read(g16pesunit,*)
        do iatom=1,natoms
        read(g16pesunit,*)&
                      (car_coordinate(ipol,iatom,i,imode),ipol=1,3)
        enddo      
      enddo
    enddo
            
    do imode=1,nmode
      write(stdout,"(10X,I5)") imode
      do iatom=1,natoms
        do ipol=1,3
          write(stdout,"(I7,5(1X,E13.6))") (iatom-1)*3+ipol,&
          (car_coordinate(ipol,iatom,1,imode)-&
           car_coordinate0(ipol,iatom)),&
           car_coordinate(ipol,iatom,2,imode)-&
           car_coordinate0(ipol,iatom),&
           car_coordinate(ipol,iatom,2,imode)-&
           car_coordinate(ipol,iatom,1,imode),&
           normal_e1(ipol,iatom,imode),&
          (car_coordinate(ipol,iatom,2,imode)-&
           car_coordinate(ipol,iatom,1,imode))/&
          (normal_e1(ipol,iatom,imode))
        enddo
      enddo
    enddo      
       
  endif   
      
  !Dipole(mu) and Polarizability(alpha) for original &
  !normal mode shift(0.01A for e1) structure.
  allocate(mu(3,2,nmode))
  allocate(alpha(3,3,2,nmode))
  mu =0.0
  alpha = 0.0
  
  do imode=1,nmode
    do i =1 ,2
      call findkline(g16fileunit," Dipole        =",1,16)
      read(g16fileunit,"(16X,3E15.8)") (mu(ipol,i,imode),ipol=1,3)
      read(g16fileunit,"(16X,3E15.8)") alpha(1,1,i,imode),&
                                       alpha(1,2,i,imode),&
                                       alpha(2,2,i,imode)
      read(g16fileunit,"(16X,3E15.8)") alpha(1,3,i,imode),&
                                       alpha(2,3,i,imode),&
                                       alpha(3,3,i,imode)
    enddo
  enddo
    
    
  ! get the index in "QUADRATIC FORCE CONSTANTS IN NORMAL MODES"
  ! get the Frequency [cm-1] in new index
  allocate(freqs_tmp(nmode),FI(nmode),mode_index(nmode))
  allocate(mode_char(nmode))
  freqs_tmp =0.0
  FI =0.0
  mode_index = 0
  mode_char = ' '
  
  call findkline(g16fileunit,&
                "QUADRATIC FORCE CONSTANTS IN NORMAL MODES",9,49)
  do i = 1,9
    read(g16fileunit,*)
  enddo
  do imode=1,nmode
    read(g16fileunit,"(A7,22X,F12.5)") mode_char(imode),FI(imode)
  enddo
  ! get the reflection for the two index
  freqs_tmp = freqs 
  mode_index = 0
  do imode=1,nmode
    inner: do jmode=1,nmode
      if(ABS(FI(imode)-freqs_tmp(jmode))<0.01) then
        mode_index(imode) = jmode
        freqs_tmp(jmode) = 0.0
        exit inner 
      endif
    enddo inner
  enddo  
   
  if(lreadchk) then
  
    allocate(freqs_org(nmode),redmass_org(nmode))
    do imode=1,nmode
      freqs_org(imode) = FI(imode)
      redmass_org(imode)= redmass(mode_index(imode))
      normal_e_org(:,:,imode) = normal_e(:,:,mode_index(imode))
      normal_e1_org(:,:,imode)= normal_e1(:,:,mode_index(imode)) 
      normal_e2_org(:,:,imode)= normal_e2(:,:,mode_index(imode)) 
    enddo  
  
    deallocate(mass)
    deallocate(normal_coord_r,normal_coord_X)
    deallocate(normal_e,normal_e2)
    deallocate(normal_e_org,normal_e1_org,normal_e2_org)      
    deallocate(freqs_org,redmass_org) 
    deallocate(car_coordinate0,car_coordinate)
    
    call close_file(g16chkname,g16chkunit)  
    
  endif
   
 
  !get CUBIC FORCE CONSTANTS IN NORMAL MODES  
  allocate(vijk(nmode,nmode,nmode))
  vijk =0.0
  
  call findkline(g16fileunit,&
                 "CUBIC FORCE CONSTANTS IN NORMAL MODES",11,47)    
  do i=1,9
    read(g16fileunit,*)
  enddo
  
  nline = 0
  lout  = .FALSE.
  do while(.NOT. lout)
    read(g16fileunit,"(A)") ctmp
    nline = nline +1
    if(trim(adjustl(ctmp))=='') lout = .TRUE. 
  enddo
  
  do iline = 1 ,nline
    backspace(g16fileunit)
  enddo
  nline = nline -1
  
  do iline=1,nline
    !vijk_org(i,j,k)
    read(g16fileunit,"(3A7,33X,f13.5)") chari,charj,chark,rtmp 
    do imode=1,nmode
      if(mode_char(imode)==chari) i = imode
      if(mode_char(imode)==charj) j = imode
      if(mode_char(imode)==chark) k = imode
    enddo
   
    i_=mode_index(i)
    j_=mode_index(j)
    k_=mode_index(k)
    vijk(i_,j_,k_) = rtmp/(redmass(i_)*redmass(j_)*redmass(k_))
    vijk(i_,k_,j_) = vijk(i_,j_,k_)
    vijk(j_,i_,k_) = vijk(i_,j_,k_)
    vijk(j_,k_,i_) = vijk(i_,j_,k_)
    vijk(k_,i_,j_) = vijk(i_,j_,k_)
    vijk(k_,j_,i_) = vijk(i_,j_,k_)
    
  enddo
  
  ! read  Cartesian Dipole derivatives in PES.out  (dmudqq)
  allocate(DDip0(3,3,natoms))
  DDip0 = 0.0
  rewind(g16pesunit)
  call findkline(g16pesunit," Cartesian Dipole derivatives:",1,30)
  read(g16pesunit,*)
  read(g16pesunit,"(4(E20.10))")(((DDip0(jpol,ipol,iatom),jpol=1,3),&
                                 ipol=1,3),iatom=1,natoms)
  
  allocate(DDip(3,3,natoms,2,nmode))
  DDip = 0.0
  do imode=1,nmode
    do i=1,2
      call findkline(g16pesunit,&
      " Cartesian Dipole derivatives:",1,30)
      read(g16pesunit,*)
      read(g16pesunit,"(4(E20.10))") &
          (((DDip(jpol,ipol,iatom,i,imode),jpol=1,3),ipol=1,3),&
          iatom=1,natoms)       
    enddo
  enddo 

  ! calculate the first dipole derivatives.
  ! dmudq
  allocate(dmudq(3,nmode))
  dmudq = 0.0
  
  do imode =1 ,nmode
    dmudq(:,imode) = (mu(:,1,imode) - mu(:,2,imode))/&
                     (sqrt(redmass(imode))*2.0*step)
  enddo
  
  ! calculate the second dipole derivatives.
  ! d2mudq2
  allocate(d2mudq2(3,nmode,nmode))
  d2mudq2 = 0.0

  do imode =1, nmode
    tmp1 = sqrt(redmass(imode))
    do jmode = 1,nmode
      tmp2 = sqrt(redmass(jmode))
      do jpol = 1,3
        rtmp = 0.0        
        rtmp = SUM((DDip(jpol,:,:,1,jmode)-DDip(jpol,:,:,2,jmode))*&
               normal_e1(:,:,imode))
        d2mudq2(jpol,jmode,imode) =rtmp /(tmp1*tmp2*step*2.0)
      enddo
    enddo
  enddo

  
  num_error = 0
  do imode = 1, nmode
    do jmode = 1, nmode
      if (jmode <= imode) then
        if ( SUM(abs(abs(d2mudq2(:,jmode,imode))-&
             abs(d2mudq2(:,imode,jmode)))) < 0.1 ) then
          write (stdout, *) jmode, imode, '   ----------'
        else if ( SUM(abs((abs(d2mudq2(:,jmode,imode)) - &
                  abs(d2mudq2(:,imode,jmode)))/ &
                  abs(d2mudq2(:,jmode,imode)))) < 0.01 ) then 
          write (stdout, *) jmode, imode, '   ----------'
        else
        write (stdout, '(2i4, 6f9.3)') jmode, imode, &
        ((d2mudq2(ipol,jmode,imode),d2mudq2(ipol,imode,jmode)),&
        ipol=1,3)
        num_error = num_error + 1
        end if
      end if
    end do
  end do

  write(stdout,*)'Second dipole derivatives done!'
  if (num_error /= 0) write (stdout, '(a, i8, a)') &
	'Warning:', num_error, 'Please check cross-derivatives!'    
  
  
  ! Calculate the first derivatives of polarizability.
  ! dadq
  allocate(dadq(3,3,nmode))
  dadq = 0.0
  
  do imode=1,nmode
    tmp1 = sqrt(redmass(imode)) * (step * 2.0)
    dadq(1,1,imode) = (alpha(1,1,1,imode)- alpha(1,1,2,imode) ) /tmp1
    dadq(1,2,imode) = (alpha(1,2,1,imode)- alpha(1,2,2,imode) ) /tmp1
    dadq(2,1,imode) = dadq(1,2,i)
    dadq(2,2,imode) = (alpha(2,2,1,imode)- alpha(2,2,2,imode) ) /tmp1
    dadq(1,3,imode) = (alpha(1,3,1,imode)- alpha(1,3,2,imode) ) /tmp1
    dadq(3,1,imode) = dadq(1,3,i)
    dadq(2,3,imode) = (alpha(2,3,1,imode)- alpha(2,3,2,imode) ) /tmp1
    dadq(3,2,imode) = dadq(2,3,i)
    dadq(3,3,imode) = (alpha(3,3,1,imode)- alpha(3,3,2,imode) ) /tmp1
  enddo
  write(stdout,*) 'First derivatives of polarizability done!'    
  
  deallocate(DDip0,DDip)
  deallocate(frc_consts,ir_inten,normal_coord)
  deallocate(normal_e1)

  call close_file(g16pesname,g16pesunit)
  call close_file(g16filename,g16fileunit)

  allocate(freqscc(nmode), w_Q1(nmode), w_Q12(nmode,nmode))
  freqscc = 0.0
  w_Q1 = 0.0
  w_Q12 = 0.0
  gamma = gamma * cc
  sigma = gamma / (2._dp * sqrt( 2._dp * log(2._dp) ))
  freqs = freqs * scaling
  
  do imode = 1,nmode
    w_Q1(imode) = freqs(imode)
    do jmode = 1, nmode
      w_Q12(imode,jmode) = freqs(imode)+freqs(jmode)
    enddo
  enddo
  
  freqscc = freqs * cc
  redmass = redmass * amu
  dmudq   = dmudq * e0
  dadq    = dadq * e0 * e0 * a0 / eh
  d2mudq2 = d2mudq2 * e0 / a0
  vijk    = vijk * eh / (a0 ** 3)    

  specunit = io_file_unit()
  call open_file(specname,specunit)
  
  do w2 = lowerw2, upperw2, stepw2
    w2cc = w2 * cc
    do w1 = lowerw1, upperw1, stepw1
      w1cc = w1 * cc
      evvhhh  = 0.0
      evvhhv  = 0.0
      evvxhvh = 0.0  
      do jmode = 1, nmode
        if ( abs(w_Q1(jmode) - w2) > threshold ) cycle
        do imode = 1, nmode
          if ( abs(w_Q12(imode,jmode) - w1) > threshold ) cycle
          call x3ac(jmode,imode,w2cc,w1cc,xhhh,xhhv,xhvh)
          evvhhh  = evvhhh  + xhhh
          evvhhv  = evvhhv  + xhhv
          evvxhvh = evvxhvh + xhvh
        end do
      end do
      evvhhh  = evvhhh  * evvhhh
      evvhhv  = evvhhv  * evvhhv
      evvxhvh = evvxhvh * evvxhvh
      if (abs(evvhhh) < 1.d-307) evvhhh = 0.0
      if (abs(evvhhv) < 1.d-307) evvhhv = 0.0
      write (specunit, '(2f10.2, 3e14.4e3)') w2, w1, evvhhh,&
                                             evvhhv,evvxhvh
    end do
  end do   

  deallocate(freqs, redmass)
  deallocate(mu,alpha)
  deallocate(freqs_tmp,FI,mode_index)
  deallocate(vijk,dmudq,d2mudq2,dadq)
  deallocate(freqscc, w_Q1, w_Q12) 
  call close_file(specname,specunit)
    

end program evv2dir

module kinds
  implicit none
  save
  !kinds definition
  integer ,parameter :: dp    = kind(1.0d0)
  integer ,parameter :: dpc   = kind((1.0d0,1.0d0))
  integer ,parameter :: i4b   = selected_int_kind(9)  
  integer ,parameter :: i8b   = selected_int_kind(18)
  integer ,parameter :: sgl   = selected_real_kind(6,30)
  integer ,parameter :: dpreal= selected_real_kind(14,200)  
end module kinds

module constants
  use kinds,only : dp
  implicit none
  integer,parameter :: maxlen = 86
  real(kind=dp),parameter :: mass_H = 1.00794
  real(kind=dp),parameter, public :: bohr2ang=0.52917721092
  real(kind=dp),parameter :: step =  0.01/bohr2ang
  real(kind=dp),parameter :: pi=3.14159265358979323846264338327950288
  real(kind=dp),parameter :: e0 = 1.602d-19 * 2.998d9
  real(kind=dp),parameter :: a0 = 52.92d-10
  real(kind=dp),parameter :: eh = 4.360d-11
  real(kind=dp),parameter :: amu = 1.661d-24
  real(kind=dp),parameter :: threshold = 200.0
  real(kind=dp),parameter :: cc = 2.998d10 * 2 * pi 
  real(kind=dp),parameter :: van = 0.0
  real(kind=dp),parameter :: NN = 1.d22
  real(kind=dp),parameter :: FF = 1.0
  
end module constants

module parameters
  use kinds, only : dp
  use constants,only:maxlen
  implicit none
  
  character(len = maxlen) :: inputfilename="evv2dir.in"
  character(len = maxlen) :: g16filename
  character(len = maxlen) :: g16chkname,g16pesname
  character(len = 2) :: gauss_revision_main
  character(len = 4) :: gauss_revision_minor
  real(kind=dp) :: lowerw1,upperw1,stepw1,lowerw2,upperw2,stepw2
  real(kind=dp) :: w1,w2
  real(kind=dp) :: gamma,scaling
  logical :: lreadchk,lreadpes

  namelist /shinput/ &
  lowerw1,upperw1,stepw1,lowerw2,upperw2,stepw2,gamma,scaling,&
  g16filename,g16chkname,lreadchk,g16pesname,lreadpes
  
end module parameters

module cal_x3
	use kinds, only : dp,i4b
	implicit none
	
	real(kind=dp), allocatable :: dadq(:,:,:)
	real(kind=dp), allocatable :: dmudq(:,:),w_Q12(:,:)
	real(kind=dp), allocatable :: freqs(:),freqscc(:)
	real(kind=dp), allocatable :: redmass(:), w_Q1(:)
  real(kind=dp), allocatable :: freqs_org(:),redmass_org(:)
  real(kind=dp), allocatable :: d2mudq2(:,:,:)
  real(kind=dp), allocatable :: vijk(:,:,:)
  real(kind=dp), allocatable :: sigma
	integer (i4b) :: nmode

  contains
 
  subroutine x3ac(c,a,walpha,wbeta,x3_zzzz,x3_xxzz,x3_xzxz)
  
  	use kinds, only  : dp, i4b
    use parameters,only : gamma
    use constants, only : van, NN, FF, pi
  	implicit none
  	real (kind=dp), intent (out) :: x3_zzzz, x3_xxzz, x3_xzxz
  	real (kind=dp), intent (in)  :: walpha, wbeta
  	integer (i4b), intent (in) :: c, a
  	real (kind=dp) :: mazzzz(nmode), maxxzz(nmode), maxzxz(nmode)
    real (kind=dp) :: eazzzz,eaxxzz, eaxzxz,pzzzz, pxxzz, pxzxz
  	real (kind=dp) :: dadQa(3,3)
  	real (kind=dp) :: tmp(6)
  	real (kind=dp) :: dudQc(3), dudQi(3), d2udQcdQa(3)
  	real (kind=dp) :: v_c, v_a, m_c, m_a, vcg, vca, vga, vcia, &
                      gamma_ca,gamma_ga,sigma_ca, sigma_ga, vi, mi
  	integer (i4b) :: n, m, i
  		
  	mazzzz = 0.0
  	maxxzz = 0.0
  	maxzxz = 0.0
  	v_c = freqscc(c)
  	m_c = redmass(c)
  	dudQc = dmudq(:,c)
  	v_a = freqscc(a)
  	m_a = redmass(a)
  	dadQa = dadq(:,:,a) !dadq(3,3,nmode)
  	vcg = v_c + v_a + van
  	vca = v_a + van
  	vga = - v_c
  	gamma_ca = gamma
  	gamma_ga = gamma
  	sigma_ca = sigma
  	sigma_ga = sigma
    
  	tmp = 0.0
  	do n = 1, 3
  		tmp(5) = tmp(5) + dadQa(n,n)
  	end do
  
  	tmp(6)= sqrt(exp(-(walpha+vga)*(walpha+vga)/2/sigma_ga/sigma_ga)&
            /sigma_ga/sqrt(2*pi) * exp(-(wbeta-walpha-vca)*         &
            (wbeta-walpha-vca)/2/sigma_ca/sigma_ca)/sigma_ca/       &
            sqrt(2*pi)) * pi / sqrt(gamma_ca * gamma_ga)
      
  	pzzzz = 0.0
  	pxxzz = 0.0
  	pxzxz = 0.0
    
  	tmp(1) = - NN * FF / (192 * m_c * m_a * v_c * v_a)
  	do i = 1, nmode
  		vi = freqscc(i)
  		mi = redmass(i)
  		dudQi = dmudq(:,i)   !dmudq(3,nmode)
  		vcia =  vijk(c,i,a)
  		tmp(2:4) = 0.0
  		do n = 1, 3
  			tmp(2) = tmp(2) + dudQi(n) * dudQc(n)
  			do m = 1, 3
  				tmp(3) = tmp(3) + dadQa(m, n) * dudQi(m) * dudQc(n)
  			end do
  		end do
      
  		tmp(4) = tmp(1) * vcia / ( mi  * (v_c + v_a - vi) * &
               (v_c + v_a + vi))
  		mazzzz(i) = tmp(4)*((tmp(5) * tmp(2) + tmp(3) * 2)/15)  *tmp(6)
  		maxxzz(i) = tmp(4)*((2 * tmp(5) * tmp(2) - tmp(3))/15)  *tmp(6)
  		maxzxz(i) = tmp(4)*((- tmp(5) * tmp(2) + tmp(3) * 3)/30)*tmp(6)
  		
      if (i /= c .and. i /= a) then
  			pzzzz = pzzzz + mazzzz(i)
  			pxxzz = pxxzz + maxxzz(i)
  			pxzxz = pxzxz + maxzxz(i)
  		end if
  	end do
  
  	d2udQcdQa = d2mudq2(:,c,a)  !d2mudq2(3,nmode,nmode)
  	tmp(2:3) = 0.0
  	do n = 1, 3
  		tmp(2) = tmp(2) + d2udQcdQa(n) * dudQc(n)
  		do m = 1, 3
  			tmp(3) = tmp(3) + dadQa(m, n) * d2udQcdQa(m) * dudQc(n)
  		end do
  	end do
  	eazzzz = tmp(1) * (tmp(5) * tmp(2) + tmp(3) * 2) / 15  * tmp(6)
  	eaxxzz = tmp(1) * (2 * tmp(5) * tmp(2) - tmp(3)) / 15 * tmp(6)
  	eaxzxz = tmp(1) * (- tmp(5) * tmp(2) + tmp(3) * 3) / 30* tmp(6)
  
  	if (c == a) then
  		x3_zzzz = mazzzz(c) + eazzzz
  		x3_xxzz = maxxzz(c) + eaxxzz
  		x3_xzxz = maxzxz(c) + eaxzxz
  	else
  		x3_zzzz = mazzzz(c) + mazzzz(a) + eazzzz
  		x3_xxzz = maxxzz(c) + maxxzz(a) + eaxxzz
  		x3_xzxz = maxzxz(c) + maxzxz(a) + eaxzxz
  	end if
  
  end subroutine x3ac

end module cal_x3

module io
  !! Module to handle operations related to file input and output.
  use kinds     ,only : dp
  use constants ,only : maxlen
  implicit none
  
  integer, public, save           :: stdin,stdout
  character(len=maxlen)           :: stdin_name
  character(len=maxlen)           :: stdout_name = "evv2dir.out"
  !! Unit on which * is written
  !integer, parameter, public      :: maxlen = 120  
  !! Max column width of input file
  integer                         :: ierr
  character(len=maxlen)           :: msg
  character(len=maxlen)           :: home_dir
  ! For parallel execution: I/O within an image
  ! These are set at startup by calling mp_world_start
  integer :: ionode_id= 0      ! index of the i/o node for this image
  Logical :: Lionode  = .True. ! true if this processor is a i/o node
                               ! for this image  
contains
   
  !========================================
  subroutine io_error ( error_msg )
  !========================================
  !! Abort the code giving an error message 
  !========================================
    implicit none
    character(len=*), intent(in) :: error_msg

    write(stdout,*)  'Exiting.......' 
    write(stdout, '(1x,a)') trim(error_msg)       
    write(stdout, '(1x,a)') trim(error_msg)
    write(stdout,'(A)') "Error: examine the output/error file" 
    STOP      
  end subroutine io_error
     
  function io_file_unit() !得到一个当前未使用的unit，用于打开文件
  !===========================================                                     
  !! Returns an unused unit number
  !! so we can later open a file on that unit.                                       
  !===========================================
    implicit none
    integer :: io_file_unit,unit_index
    logical :: file_open

    unit_index = 9
    file_open  = .true.
    do while ( file_open )
      unit_index = unit_index + 1
      inquire( unit=unit_index, OPENED = file_open ) 
    end do
    io_file_unit = unit_index
    return
  end function io_file_unit
  
  subroutine open_file(file_name,file_unit)
    implicit none 
    character(len=*),intent(in) :: file_name
    integer,intent(in)          :: file_unit
    open(unit=file_unit, file=file_name,iostat=ierr,iomsg=msg)
    if(ierr /= 0 ) then
      call io_error('Error: Problem opening "'//&
                     trim(adjustl(file_name))//' " file')
      call io_error(msg)
    endif
  end subroutine open_file
  
  subroutine close_file(file_name,file_unit)
    implicit none
    integer,intent(in)          :: file_unit
    character(len=*),intent(in) :: file_name
    close(file_unit,iostat=ierr,iomsg=msg)
    if(ierr /= 0 ) then
      call io_error('Error: Problem close "'//&
                     trim(adjustl(file_name))//' " file')
      call io_error(msg)
    endif   
  end subroutine close_file
    
  subroutine findkword(funit,kword)
    implicit none
    integer,intent(in):: funit
    character(len=*),intent(in) :: kword
    logical :: lfindkword 
    character(len=maxlen) :: ctmp
    
    lfindkword = .FAlSE.
    do while(.NOT. lfindkword)
      read(funit,*) ctmp
      if (trim(adjustl(ctmp))==trim(adjustl(kword))) then
        lfindkword = .TRUE.
        backspace(unit=funit)
      endif
    enddo  
  end subroutine findkword

  subroutine findkline(funit,kline,indexi,indexe)
    implicit none
    integer,intent(in):: funit
    character(len=*),intent(in) :: kline
    integer,intent(in)::indexi,indexe
    logical :: lfindkline
    character(len=maxlen) :: ctmp
    
    lfindkline = .FAlSE.
    do while(.NOT. lfindkline)
      read(funit,"(A)") ctmp
      if (ctmp(indexi:indexe)==kline) then
        lfindkline = .TRUE.
        backspace(unit=funit)
      endif
    enddo  
  end subroutine findkline
   
end module io

module readinput
  use kinds,only : dp
  use constants,only : maxlen,step
  use parameters,only :lowerw1, upperw1, stepw1, lowerw2, upperw2,& 
                       stepw2, gamma, scaling,g16filename,        &
                       g16chkname,lreadchk,g16pesname,lreadpes
  implicit none
  integer               :: in_unit,tot_num_lines,ierr,loop,in1,in2
  integer               :: num_lines,line_counter
  character(len=maxlen),allocatable:: in_data(:) 
  character(len=maxlen) :: dummy,ctmp
  integer               :: ipos
  character, parameter  :: TABCHAR = char(9) !char(9)为制表符TAB    
  
  contains
  
  subroutine get_inputfile(filename)
    implicit none
    character(len=*),intent(in) :: filename
    call treat_inputfile(filename)
    call read_namelist()
  end subroutine

  subroutine treat_inputfile(filename)  
  !用于将输入文件中的每一条写入字符串文件
  !并且(将非文件路径)改为小写，去除注释          
  !=======================================!
  !! Load the shin file into a character  
  !! array in_file, ignoring comments and  
  !! blank lines and converting everything 
  !! to lowercase characters               
  !=======================================!
    use io,only : io_file_unit,io_error
    implicit none
    character(*),intent(in):: filename
    logical :: alive
    character(len=maxlen) :: msg
    
    in_unit=io_file_unit( )
    inquire(file=trim(adjustl(filename)),exist=alive)
    if(.NOT. alive) then
      call io_error("Error:Input file "//trim(adjustl(filename))&
                    //" doesn't exist.")
    else
      open (unit=in_unit, file=trim(adjustl(filename)),&
            form='formatted',status='old',iostat=ierr)
      if(ierr /= 0) then
        call io_error('Error: Problem opening input file LVCSH.in')
      endif
    endif
    
    num_lines=0;tot_num_lines=0;ierr=0
    do while( ierr == 0 )
      read(in_unit, '(a)', iostat = ierr ,iomsg=msg) dummy 
      if(ierr > 0 ) then
        call io_error('Error: Problem reading input file SHIN')
        call io_error(msg)
      elseif(ierr == 0 )then
        ! convert all tabulation characters to spaces
        ipos = index(dummy,TABCHAR) 
        !查询字符串在字符串中出现的位置,并将制表符改为空格
        do while (ipos /= 0)
          dummy(ipos:ipos) = ' '
          ipos = index(dummy,TABCHAR)
        end do
        dummy=adjustl(dummy)
        tot_num_lines=tot_num_lines+1
        if( dummy(1:1)/='!'  .and. dummy(1:1)/='#' ) then
          if(len_trim(adjustl(dummy)) > 0 ) num_lines=num_lines+1
        endif
      endif
    end do
    !得到SHIN文件中总的行数tot_num_lines以及非注释和空行 num_lines

    rewind(in_unit)
    allocate(in_data(num_lines),stat=ierr)  
    !字符串数组，内部文件 line=449
    if (ierr/=0) call io_error('Error allocate &
                       in_data in para_in_file')
    line_counter=0
    do loop=1,tot_num_lines
      read(in_unit, '(a)', iostat = ierr ,iomsg=msg) dummy
      if(ierr /= 0) then
        call io_error('Error: Problem opening input file SHIN')
        call io_error(msg)
      endif
      !I convert all tabulation characters to spaces
      ipos = index(dummy,TABCHAR)
      do while (ipos /= 0)
        dummy(ipos:ipos) = ' '
        ipos = index(dummy,TABCHAR)
      end do
      !if(index(dummy,"dir")<=0) then 
      !  dummy=utility_lowercase(dummy) 
      !将(不表示文件路径的)字符串中大写字母全部改为小写
      !endif
      dummy=trim(adjustl(dummy))     !
      if( dummy(1:1)=='!' .or.  dummy(1:1)=='#' ) cycle
      if(len(trim(dummy)) == 0 ) cycle
      if(index(dummy,'=') <=1 )  cycle  
      !当该行中没有‘=’ 或‘=’前没有内容则跳过该行
      line_counter=line_counter+1
      
      !去除有效行信息中的注释部分，注释可以采用 ！或者 #
      in1=index(dummy,'!')
      in2=index(dummy,'#')
      if(in1==0 .and. in2==0)  in_data(line_counter)=dummy
      !不存在'!'与'#'
      if(in1==0 .and. in2>0 )  in_data(line_counter)=dummy(:in2-1)
      if(in2==0 .and. in1>0 )  in_data(line_counter)=dummy(:in1-1)
      if(in2> 0 .and. in1>0 )  in_data(line_counter)=&
                               dummy(:min(in1,in2)-1)  
    end do
    !得到包含有效信息的行数line_counter,和相应的数据in_data(line_counter)
    close(in_unit)
  end subroutine treat_inputfile 

  subroutine read_namelist()
    use parameters,only : shinput
    use io
    implicit none
    integer::incar_unit,i
    
    !!write input file to namelist input file
    incar_unit = io_file_unit()
    open(unit=incar_unit,status='SCRATCH',iostat=ierr,iomsg=msg)
    if(ierr > 0 ) then
      call io_error('Error: Problem reading SCRATCH input namelist')
      call io_error(msg)
    elseif(ierr == 0 )then    
      write(incar_unit,*)"&shinput" 
      do i=1,line_counter
        write(incar_unit,*) trim(adjustl(in_data(i)))
      enddo
      write(incar_unit,*) "/"
    endif
    rewind(incar_unit)
    
    !set default values for variables in namelist
    lowerw1 = 500
    upperw1 = 2000
    stepw1  = 5
    lowerw2 = 500
    upperw2 = 2000
    stepw2  = 5
    gamma   = 20
    scaling     = 0.96
    g16filename = "gaussian.log"
    g16chkname  = "gaussian.chk"
    lreadchk    = .true.
    g16pesname  = "PES.out"
    lreadpes    = .TRUE.
 
    write(stdout,"(/,1X,A67)")   repeat("=",67)
    write(stdout,"(1X,10X,A)") "The namelist file as follows"
    write(stdout,"(1X,A67)")   repeat("=",67)
    do i=1,line_counter+2
      read(incar_unit,"(A80)") ctmp
      write(stdout,"(A80)") ctmp
    enddo
    rewind(incar_unit)
    read(UNIT=incar_unit,nml=shinput,iostat=ierr,iomsg=msg)
    if(ierr /= 0) then
      call io_error('Error: Problem reading namelist file SHIN')
      call io_error(msg)
    endif  
    close(incar_unit)
    
  end subroutine read_namelist

end module readinput

