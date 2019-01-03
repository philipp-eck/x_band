	Program x_band_PROCAR 
	! Program is written to extract the orbital character from the PROCAR file, writes the output band ordered. K values in A/2pi
	! compile with ifort -free -heap-arrays 1 -o extract_char_band.x extract_char_band.f
	Implicit None

	integer*8 :: ikpt, i, j, n, ispin, norb, iorb, iband, bandi, iatom, nat, l
	integer*8 :: Nbands, Nkpt, count_max
	integer*8 :: Natoms
	integer*8, allocatable :: iband_sort(:,:)
	integer*8, allocatable :: AtomArray(:), dummyarray(:)
	real*8 :: dummy, SurfChar, Chari, rec_vec(3,3), k_car(3), k_car_prev(3), scal, real_vec(3,3)
	real*8 :: dist, kpt(3)
	

	character(len=4), allocatable :: header(:)
	character(len=80) :: TXT, input, filetype, output, vaspversion, orb 
	character(len=40) :: fs
	real*8, allocatable :: SquaredAmplitude(:,:), kdist(:), eval(:,:)
	real*8, allocatable :: OrbArray(:),tot(:,:), m(:,:,:)
	real*8, allocatable :: OrbChar(:,:,:)
	logical :: LSORBIT,ok_flag, ISORT
	

	allocate(header(16))
	header = [character(len=65) :: 's','py','pz','px','dxy','dyz','dz2',&
	      'dxz','dx2','f-3','f-2','f-1','f0','f1','f2','f3']

	! Read input x_input.dat
	open( unit=12,file='x_input.dat', status='old', form='formatted' )
	read( unit=12,nml=FLAGS )
	close( unit=12 )

	Namelist/FLAGS/ LSORBIT, &   !SOC flag
	                VASPVERSION, & !VASP.x.x.x
	                ORB, &         ! "d" or "f" electron system
	                INPUT, &
	                OUTPUT, &
	                ISORT         ! Option, switch on eigenvalue sorting for hybrid comps

	! Read Natoms and allocate AtomArray
	open(unit=13,file=input, status='old')
	read(13,*)
	read(13,*) txt, txt, txt, Nkpt, txt, txt, txt, Nbands, txt, txt, txt, Natoms
	close(13)

	Namelist/ATARRAY/ ATOMARRAY !Array containing atoms

	! Allocate and read Atomarray
	allocate(AtomArray(Natoms))
	AtomArray = 0.0
	open( unit=12,file='x_input.dat', status='old', form='formatted' )
	read( unit=12,nml=ATARRAY )
	close( unit=12 )

	! Specify the atoms 
	! Compute reciprocal lattice vectors from the CONTCAR
	open(unit=11,file='CONTCAR',status='old')
	
	read(11,*)
	read(11,*) scal
	do i=1,3
	 read(11,*) (real_vec(i,j), j=1,3)
	end do
	
	real_vec=scal*real_vec
	call M33INV(real_vec,rec_vec,ok_flag)

	! Read PROCAR
	open(unit=10,file=input,status='old')
	read(10,'(A)') filetype
	read(10,*) txt, txt, txt, Nkpt, txt, txt, txt, Nbands, txt, txt, txt, Natoms
	! allocate  arrays
	if (orb .EQ. "d") then
	 norb = 10
	else if (orb .EQ. "f") then
	 norb = 17
	else
	 write(*,*) "Unknown ORB, choose 'd' or 'f'"
	end if
	allocate(m(Nkpt,Nbands,3))
	allocate(SquaredAmplitude(Natoms,norb))
	allocate(dummyarray(norb-1))
	allocate(tot(Nkpt,Nbands))
	allocate(OrbChar(Nkpt,Nbands,norb))
	allocate(kdist(Nkpt))
	allocate(eval(Nkpt,Nbands))
	allocate(iband_sort(Nkpt,Nbands))
	dist = 0.0000000000d0
	do ikpt=1,Nkpt
	!write(*,'(A,1i6)') 'Reading k-point: ', ikpt
	 read(10,*)
	!read(10,*) txt, dummy, txt, (kpt(i), i=1,3)
	 read(10,'(A19,3f11.8)') txt, (kpt(i), i=1,3) 
	!write(*,*) kpt(:)
	 k_car=matmul(rec_vec,kpt(:))
	 if ( ikpt > 1) then
	  dist = dist + sqrt( (k_car(1)-k_car_prev(1))**2 + (k_car(2)-k_car_prev(2))**2 + (k_car(3)-k_car_prev(3))**2 )
	! write(*,*) "Compute dist"
	 end if
	!write(*,*) dist, kpt
	 kdist(ikpt) = dist 
	 read(10,*)
	 do iband=1,Nbands
	  read(10,*) txt, dummy, txt, txt, eval(ikpt,iband)
	! write(*,*)'eval',  eval(ikpt,iband)
	  read(10,*)
	  read(10,*)
	  do iatom=1,Natoms
	   read(10,*) dummy, (SquaredAmplitude(iatom,iorb), iorb=1,norb)
          end do
	  do i = 1,10
	   Chari = 0.d0
	   do j = 1,SIZE((AtomArray))
	    if (AtomArray(j) .NE. 0) then
	     Chari = Chari + SquaredAmplitude(AtomArray(j),i)
	    end if
	   end do
	   OrbChar(ikpt,iband,i) = Chari
	  end do
	  read(10,*) txt, (dummyarray(iorb), iorb=1,norb-1), tot(ikpt,iband)
	  if (vaspversion =='5.4.1') then
	   read(10,*)
	  end if
	  if (LSORBIT ==.TRUE.) then
	   do l=1,3
	    do iatom=1,Natoms
	     read(10,*)
	    end do
 	    read(10,*) txt, (dummyarray(iorb), iorb=1,norb-1), m(ikpt,iband,l)
	    if (vaspversion =='5.4.1') then
	     read(10,*)
	    end if
	   end do
	  end if
	  if (filetype =='PROCAR lm decomposed + phase') then
	   if (vaspversion =='5.4.1') then
	    count_max = 2*Natoms+1
	   else if (vaspversion == '5.4.4') then
	    count_max = Natoms+2
	   end if
 	   do iatom=1,count_max 
	    read(10,*)
	   end do
	  end if
	! write(*,*) "finished reading loop"
	 end do
	 read(10,*)
	 k_car_prev(1) = k_car(1)
	 k_car_prev(2) = k_car(2)
	 k_car_prev(3) = k_car(3)
	end do
	close(10)
	
	!Sorting the eigenvalues with respect to their energies
	if (ISORT==.TRUE.) then
	call SORTIBAND(eval, Nkpt, Nbands, iband_sort)
	end if

	! Start writing the output
	
	! Write header
	open(unit=200,file=output)
	if (LSORBIT ==.TRUE.) then
	 write(fs,'(A11,i2,A4)') '(A3,A9,A10,',norb+5,'A7)'
	 write(200,fs)'i','kdist','eval','mx','my','mz',(trim(header(iorb)),iorb=1,norb-1),'tot','Abs'
	else
	 write(fs,'(A11,i2,A4)') '(A3,A9,A10,',norb+2,'A7)'
	 write(200,fs) 'i', 'kdist', 'eval',(trim(header(iorb)),iorb=1,norb-1), 'tot','Abs'
	end if
	
	!Write projections
	do iband=1, Nbands
	 do ikpt=1, Nkpt
	  if (ISORT==.TRUE.) then
	   j = iband_sort(ikpt,iband)
	  else
	   j = iband
	  end if
	  if (LSORBIT ==.TRUE.) then
	   write(fs,'(A15,i2,A6)') '(i3,f9.5,f10.5,',norb+4,'f7.3)'
	   write(200,fs) iband, kdist(ikpt), eval(ikpt,j), (m(ikpt,j,i), i=1,3), (OrbChar(ikpt,j,i), i= 1,norb), tot(ikpt,j)
	  else
	   write(fs,'(A15,i2,A6)') '(i3,f9.5,f10.5,',norb+1,'f7.3)'
	   write(200,fs) iband, kdist(ikpt), eval(ikpt,j), (OrbChar(ikpt,j,i), i= 1,norb),tot(ikpt,j)
	  end if
	 end do
	 write(200,*)
	end do
	close(200)


	CONTAINS

	 !Sort the band index wrt to the energy
	SUBROUTINE SORTIBAND (eval, Nkpt, Nbands, sort_iband)
	
	implicit none
	integer*8, intent(in) :: Nkpt, Nbands
	integer*8 :: ikpt, iband, s
	integer*8, intent(out) :: sort_iband(Nkpt,Nbands)
	real*8, intent(in) :: eval(Nkpt,Nbands)
	logical :: SORT
	
	do ikpt=1,Nkpt
	 sort_iband(ikpt,1)=1
	 do iband=2,Nbands
	! write(*,*) 'sort kpoint band', ikpt, iband
	  if (eval(ikpt,iband) .lt. eval(ikpt,sort_iband(ikpt,iband-1))) then
	   SORT = .FALSE.
	!  write(*,*)'start sorting'
	   s = 0
	   do while (SORT == .FALSE.)
	    sort_iband(ikpt,iband-s)=sort_iband(ikpt,iband-s-1)
	    if (iband-s-1==1) then
	     SORT = .TRUE.
	    else if (eval(ikpt,sort_iband(ikpt,iband-s-2))<eval(ikpt,iband)) then
	     SORT = .TRUE.
	    end if
	    s = s + 1
	   end do
	   sort_iband(ikpt,iband-s) = iband
	  else
	   sort_iband(ikpt,iband) = iband
	  end if
	 end do
	end do
	return
	END SUBROUTINE SORTIBAND

	SUBROUTINE M33INV (A, AINV, OK_FLAG)
	
	IMPLICIT NONE
	
	DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
	DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
	LOGICAL, INTENT(OUT) :: OK_FLAG
	
	DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
	DOUBLE PRECISION :: DET
	DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR
	
	
	DET =   A(1,1)*A(2,2)*A(3,3)  &
	      - A(1,1)*A(2,3)*A(3,2)  &
	      - A(1,2)*A(2,1)*A(3,3)  &
	      + A(1,2)*A(2,3)*A(3,1)  &
	      + A(1,3)*A(2,1)*A(3,2)  &
	      - A(1,3)*A(2,2)*A(3,1)
	
	IF (ABS(DET) .LE. EPS) THEN
	   AINV = 0.0D0
	   OK_FLAG = .FALSE.
	   RETURN
	END IF
	
	COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
	COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
	COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
	COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
	COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
	COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
	COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
	COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
	COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))
	
	AINV = TRANSPOSE(COFACTOR) / DET
	
	OK_FLAG = .TRUE.
	
	RETURN
	
	END SUBROUTINE M33INV

	End Program x_band_PROCAR
