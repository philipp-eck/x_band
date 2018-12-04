	Program x_band_PROCAR 
	! Program is written to extract the orbital character from the PROCAR file, writes the output band ordered. K values in A/2pi
	! compile with ifort -free -heap-arrays 1 -o extract_char_band.x extract_char_band.f
	Implicit None

	integer*8 :: ikpt, i, j, n, ispin, iorb, iband, bandi, iatom, nat, l
	integer*8 :: Nbands, Nkpt
	integer*8 :: Natoms, AtomArray(100)
	integer*8, allocatable :: iband_sort(:,:)
	!integer*8, allocatable :: AtomArray(:)
	real*8 :: dummy, SurfChar, Chari, rec_vec(3,3), k_car(3), k_car_prev(3), scal, real_vec(3,3)
	real*8 :: dist, kpt(3)
	

	CHARACTER(LEN=80) :: TXT, input, filetype, output
	real*8, allocatable :: SquaredAmplitude(:,:), kdist(:), eval(:,:)
	real*8, allocatable :: OrbArray(:),tot(:,:), m(:,:,:)
	real*8, allocatable :: OrbChar(:,:,:)
	logical :: LSORBIT,ok_flag, ISORT
	

	! Read Natoms and allocate AtomArray
	open(unit=13,file='PROCAR', status='old')
	read(13,*)
	read(13,*) txt, txt, txt, Nkpt, txt, txt, txt, Nbands, txt, txt, txt, Natoms
	close(13)
	!allocate(AtomArray(Natoms))

	Namelist/FLAGS/ LSORBIT, &   !SOC flag
	                ATOMARRAY, & !Array containing atoms
	                INPUT, &
	                OUTPUT, &
	                ISORT         ! Option, switch on eigenvalue sorting for hybrid comps

	! Read input x_input.dat
	open( unit=12,file='x_input.dat', status='old', form='formatted' )
	read( unit=12,nml=FLAGS )
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
	allocate(m(Nkpt,Nbands,3))
	allocate(SquaredAmplitude(Natoms,10))
	allocate(tot(Nkpt,Nbands))
	allocate(OrbChar(Nkpt,Nbands,10))
	allocate(kdist(Nkpt))
	allocate(eval(Nkpt,Nbands))
	allocate(iband_sort(Nkpt,Nbands))
	dist = 0.0000000000d0
	do ikpt=1,Nkpt
	!write(*,'(A,1i6)') 'Reading k-point: ', ikpt
	 read(10,*)
	!read(10,*) txt, dummy, txt, (kpt(i), i=1,3)
	 read(10,'(A18,3f11.8)') txt, (kpt(i), i=1,3) 
	 write(*,*) kpt(:)
	 k_car=matmul(rec_vec,kpt(:))
	 if ( ikpt > 1) then
	  dist = dist + sqrt( (k_car(1)-k_car_prev(1))**2 + (k_car(2)-k_car_prev(2))**2 + (k_car(3)-k_car_prev(3))**2 )
	  write(*,*) "Compute dist"
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
	   read(10,*) dummy, (SquaredAmplitude(iatom,iorb), iorb=1,10)
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
	  read(10,*) txt, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, tot(ikpt,iband)
	  if (LSORBIT ==.TRUE.) then
	   do l=1,3
	    do iatom=1,Natoms
	     read(10,*)
	    end do
 	    read(10,*) txt, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, m(ikpt,iband,l)
	   end do
	  end if
	  if (filetype =='PROCAR lm decomposed + phase') then
 	!   do iatom=1,2*Natoms+1 for VASP.5.4.1
	   do iatom=1,Natoms+2 ! for VASP.5.4.4
	    read(10,*)
	   end do
	  end if
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
	 write(200,'(A6,A12,1A15,15A10)')'iband','kdist','eval','mx','my','mz','s','py','pz','px','dxy', 'dyz','dz2','dxz','dx2','tot','Abs'
	else
	 write(200,'(A6,A12,1A15,15A10)') 'iband', 'kdist', 'eval', 's', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'tot','Abs'
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
	   write(200,'(i6,f12.5,f15.8,14f10.3)') iband, kdist(ikpt), eval(ikpt,j), (m(ikpt,j,i), i=1,3), (OrbChar(ikpt,j,i), i= 1,10), tot(ikpt,j)
	  else
           write(200,'(i6,f12.5,f15.8,11f10.3)') iband, kdist(ikpt), eval(ikpt,j), (OrbChar(ikpt,j,i), i= 1,10),tot(ikpt,j)
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
