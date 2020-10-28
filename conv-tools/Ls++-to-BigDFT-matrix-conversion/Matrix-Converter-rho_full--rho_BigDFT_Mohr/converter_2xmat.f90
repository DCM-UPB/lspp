program convert_densitymatrix
  implicit none
  integer :: irow, icol, iat, nat, nspin, nrow, ncol, ndata
  character(len=128) :: cc, filename
  real(kind=8),dimension(:,:),allocatable :: matrix

  write(*,*) 'enter filename, number of atoms and number of spin components'
  read(*,*) filename, nat, nspin
  open(unit=1,file=trim(filename))
  read(1,*) cc, cc
  read(1,*) cc, ncol
  read(1,*) cc, nrow
  read(1,*) cc, ndata
  if (nrow/=ncol) stop 'nrow/=ncol'
  if (nrow*ncol/=ndata) stop 'nrow*ndata/=ndata'
  allocate(matrix(nrow,ncol))
  do irow=1,nrow
      read(1,*) (matrix(irow,icol), icol=1,ncol)
  end do
  close(unit=1)

  open(unit=1,file=trim(filename)//'_converted')
  write(1,'(a,3i9,a)') '#', ncol, nat, nspin, '  number of basis functions, number of atoms, number of spins'
  do iat=1,nat
      write(1,'(a,3es16.8,a)') '#', 0.0d0, 0.d0, 0.0d0, '  fake atomic positions'
  end do
  do irow=1,nrow
      do icol=1,ncol
          write(1,'(2i8,es24.15,2i8)') irow, icol, 2*matrix(irow,icol), 0, 0
      end do
  end do
  close(unit=1)
  deallocate(matrix)

end program convert_densitymatrix

