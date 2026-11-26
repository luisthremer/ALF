! Unit test for the Hubbard potential energy Pot_scal on bilayer lattices (KNOWN_ISSUES.md H7).
! Usage: 44-Hubbard-pot-weights <Lattice_type>   (Bilayer_square or Bilayer_honeycomb)
! With a Green function G = 0 every Grc(i,i) is 1, so Pot_scal = N_cells * sum_orbitals weight(orbital):
! the first Norb/2 orbitals carry Ham_U, the others Ham_U2 (as the interaction vertices in Ham_V).
! Ham_U = 1, Ham_U2 = 3, 2x2 unit cells: Bilayer_square (Norb 2) 4*(1+3) = 16, Bilayer_honeycomb (Norb 4) 4*(2+6) = 32.
Program Test_Hubbard_pot_weights
#ifdef MPI
  Use mpi
#endif
  Use Hamiltonian_main
  Implicit None

  Character (len=64) :: name
  Integer :: unit_para
  Real (Kind=Kind(0.d0)) :: expected, pot
  Complex (Kind=Kind(0.d0)), allocatable :: GR(:,:,:)
  Complex (Kind=Kind(0.d0)) :: z
  Integer :: n_col
  Real (Kind=Kind(0.d0)) :: sgn
  Character (len=64) :: ham_name
#ifdef MPI
  Integer :: ierr
  Call MPI_Init(ierr)
  Call Set_Group_Comm(MPI_COMM_WORLD)
#endif
#ifdef HDF5
  Write(*,*) 'SKIP: needs the legacy text output (built without HDF5)'
#else

  Call get_command_argument(1, name)
  Select Case (trim(name))
  Case ("Bilayer_square")
     expected = 16.d0
  Case ("Bilayer_honeycomb")
     expected = 32.d0
  Case Default
     Write(*,*) 'usage: 44-Hubbard-pot-weights Bilayer_square|Bilayer_honeycomb'
     error stop 1
  end Select

  Open(newunit=unit_para, file='parameters', status='replace', action='write')
  Write(unit_para,'(A)') '&VAR_lattice'
  Write(unit_para,'(A)') ' Lattice_type = "' // trim(name) // '", L1 = 2, L2 = 2'
  Write(unit_para,'(A)') '/'
  Write(unit_para,'(A)') '&VAR_Model_Generic /'
  Write(unit_para,'(A)') '&VAR_Hubbard'
  Write(unit_para,'(A)') ' Ham_U = 1.d0, Ham_U2 = 3.d0'
  Write(unit_para,'(A)') '/'
  Close(unit_para)

  ham_name = 'Hubbard'
  Call Alloc_Ham(ham_name)
  Call ham%Ham_set()
  Call Validate_Ham_Variables()
  Call Setup_flavor_map()
  Call ham%Alloc_obs(0)
  Call ham%Init_obs(0)          ! zero the accumulators, as main does before every bin

  allocate(GR(Ndim,Ndim,N_FL))
  GR = (0.d0, 0.d0)
  Call ham%Obser(GR, (1.d0, 0.d0), 1, 1.d0)
  ! Pot_scal is not accessible from outside the model: read it from the legacy text output
  Call ham%Pr_obs(0)
  Open(newunit=unit_para, file='Pot_scal', status='old', action='read')
  Read(unit_para,*) n_col, z, sgn
  Close(unit_para)
  pot = real(z, kind(0.d0))

  if (abs(pot - expected) > 1.d-12) then
     Write(*,'(A,A,A,ES14.6,A,ES14.6)') 'FAIL: ', trim(name), ' Pot_scal = ', pot, ' expected ', expected
     error stop 1
  endif
  Write(*,'(A,A,A,ES14.6)') 'PASS: ', trim(name), ' Pot_scal = ', pot
#endif
#ifdef MPI
  Call MPI_Finalize(ierr)
#endif
end Program Test_Hubbard_pot_weights
