program orotron
   use, intrinsic :: iso_c_binding
   use types

   implicit none

   interface
      subroutine oro(INP, OUTP) bind(c, name='oro')
         use types
         type(Input), target, intent(in) :: INP
         type(Output), target, intent(inout) :: OUTP
      end subroutine
   end interface

   type(Input) INP
   type(Output) OUTP

   namelist /param/ Ne, Lz, Tend, Delta, Ic, dz, dt

   real(c_double) Lz, Tend, Delta, dz, dt, Ic
   real(c_double), allocatable :: ZAxis(:), TAxis(:), ZBEG, ZEND

   integer(c_int) i, Ne, Nz, Nt, err_alloc, INTERVALT, INTERVALZ, OUTNz, OUTNt

   complex(c_double_complex), allocatable :: InitialField(:)

   logical(c_bool) IND1

   open (unit=1, file='input_fortran.in', status='old', err=101)
   read (unit=1, nml=param, err=102)
   close (unit=1)

   write (*, nml=param)

   Nz = Lz/dz ! Razmer vektora ZAxis (vvod)
   Nt = Tend/dt ! Razmer vektora TAxis (vvod)
   allocate (ZAxis(0:Nz), TAxis(0:Nt), InitialField(0:Nz), INP%ZAxis(0:Nz), INP%TAxis(0:Nt), INP%InitialField(0:Nz), stat=err_alloc)
   if (err_alloc /= 0) then
      print *, "allocation error"
      pause
      stop
   end if

   ! Vychislit' razmery massivov dlya vyvoda
   INTERVALT = Nt/500
   INTERVALZ = Nz/500
   !INTERVALT = 1
   !INTERVALZ = 1

   if (INTERVALT < 1) then
      print *, 'Too small "Tend"'
      pause
      stop
   end if

   if (INTERVALZ < 1) then
      print *, 'Too small "Lz"'
      pause
      stop
   end if

   if (INTERVALT > 1 .and. INTERVALZ > 1) then
      OUTNt = Nt/INTERVALT
      OUTNz = Nz/INTERVALZ
   elseif (INTERVALT == 1 .and. INTERVALZ > 1) then
      OUTNt = Nt
      OUTNz = Nz/INTERVALZ
   elseif (INTERVALT > 1 .and. INTERVALZ == 1) then
      OUTNt = Nt/INTERVALT
      OUTNz = Nz
   else ! (INTERVALT == 1) .and. (INTERVALZ == 1) then
      OUTNt = Nt
      OUTNz = Nz
   end if

   ! Vydelit' pamyat' pod massivy vyvoda
   allocate (OUTP%OUTB(0:OUTNz, 0:OUTNt), OUTP%OUTCu(0:OUTNz, 0:OUTNt), stat=err_alloc)
   if (err_alloc /= 0) then
      print *, "allocation error"
      pause
      stop
   end if

   do i = 0, Nz
      ZAxis(i) = i*dz
   end do

   do i = 0, Nt
      TAxis(i) = i*dt
   end do

   ZBEG = 0
   ZEND = 0.5
   where ((ZAxis .GE. ZBEG) .AND. (ZAxis .LE. ZEND)) InitialField = 0.001*sin(pi*(ZAxis - ZBEG)/(ZEND - ZBEG))**2; 
   INP%Ne = Ne
   INP%ZAxis(:) = ZAxis(:)
   INP%TAxis(:) = TAxis(:)
   INP%Delta = Delta
   INP%Ic = Ic
   INP%INTERVALT = INTERVALT
   INP%INTERVALZ = INTERVALZ
   INP%InitialField(:) = InitialField
   INP%OUTNz = OUTNz
   INP%Nz = Nz
   INP%Nt = Nt
   INP%dz = dz
   INP%dt = dt

   call oro(INP, OUTP)

   !open(1, file = 'test.dat', err = 101)
   !write(1, '(4I)', err = 103)  size(TAxis,1), INTERVALT, size(ZAxis,1), INTERVALZ
   !do i=1,Nz
   !write(1, '(3F17.8)', err = 103) ZAxis(i), real(InitialField(i)), imag(InitialField(i))
   !enddo
   !close(1)

   pause
   stop
101 print *, "error of file open"; pause; stop
102 print *, 'error of file read'; pause; stop
103 print *, 'error of file write'; pause; stop
end
