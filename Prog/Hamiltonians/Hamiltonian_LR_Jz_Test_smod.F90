!  Copyright (C) 2016 - 2020 The ALF project
!
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
!
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version
!
!  New file of the delayed/submatrix research branch (feature/delayed-submatrix, 2026). See docs/algorithm-notes.md.


!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> TEST-ONLY model of the low-rank (delayed/submatrix) validation suite, not a physical production model.
!> Mz-Hubbard model on the square lattice (N_FL = 2, N_SUN = 1) with an additional Ising-type spin-spin
!> interaction on every nearest-neighbour bond,
!>   H = -t sum_<ij>,s c^dag_{is} c_{js}  - mu N  + U sum_i (n_{i,up}-1/2)(n_{i,do}-1/2)
!>       + Jz sum_<ij> S^z_i S^z_j,      S^z_i = (n_{i,up} - n_{i,do})/2.
!> Note: the Predefined_Int_Jz decomposition simulates -(|Jz|/2) (S^z_i - sgn(Jz) S^z_j)^2 per bond, i.e. the Jz term
!> above plus on-site (S^z_i)^2 terms (an effective shift of U by about Jz per bond pair); Pot and Ener measure the H
!> written above only, not that extra term. Irrelevant for the scheme comparisons this model exists for.
!> The on-site term uses Predefined_Int_U_MZ (diagonal rank-1 vertices). The bond term uses Predefined_Int_Jz
!> (type 2, N_P = 2 multi-site vertices acting on sites i and j, one per nearest-neighbour bond;
!> O = diag(1, -sgn Jz), alpha = 0, g = +-sqrt(dtau |Jz| / 8) for the two flavours).
!> Operator order in Op_V(:,nf): the Ndim on-site vertices first, then for each site I the x-bond (I, I+a_1)
!> followed by the y-bond (I, I+a_2). No Symm decomposition, no checkerboard, no projector.
!> The purpose of the model is to give the delayed/submatrix schemes multi-site vertices of the
!> Predefined_Int_Jz class to be validated on.
!--------------------------------------------------------------------


    submodule (Hamiltonian_main) ham_LR_Jz_Test_smod

      Use Operator_mod
      Use WaveFunction_mod
      Use Lattices_v3
      Use MyMats
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      Use Observables
      Use Fields_mod
      Use Predefined_Hoppings
      use runtime_error_mod


      Implicit none

      type, extends(ham_base) :: ham_LR_Jz_Test
      contains
        ! Set Hamiltonian-specific procedures
        procedure, nopass :: Ham_Set
        procedure, nopass :: Alloc_obs
        procedure, nopass :: Obser
        procedure, nopass :: ObserT
#ifdef HDF5
        procedure, nopass :: write_parameters_hdf5
#endif
      end type ham_LR_Jz_Test

      !#PARAMETERS START# VAR_lattice
      Character (len=64) :: Model = ''  ! Value irrelevant
      Character (len=64) :: Lattice_type = 'Square'  ! Possible Values: 'Square'
      Integer            :: L1 = 4   ! Length in direction a_1
      Integer            :: L2 = 4   ! Length in direction a_2
      !#PARAMETERS END#

      !#PARAMETERS START# VAR_LR_Jz_Test
      real(Kind=Kind(0.d0)) :: ham_T    = 1.d0      ! Hopping parameter
      real(Kind=Kind(0.d0)) :: ham_U    = 4.d0      ! Hubbard interaction
      real(Kind=Kind(0.d0)) :: ham_Jz   = 1.d0      ! Ising spin-spin coupling on nearest-neighbour bonds (Predefined_Int_Jz)
      real(Kind=Kind(0.d0)) :: Ham_chem = 0.d0      ! Chemical potential
      real(Kind=Kind(0.d0)) :: Dtau     = 0.1d0     ! Thereby Ltrot=Beta/dtau
      real(Kind=Kind(0.d0)) :: Beta     = 2.d0      ! Inverse temperature
      !#PARAMETERS END#

      Type (Lattice),       target :: Latt
      Type (Unit_cell),     target :: Latt_unit

    contains

      module Subroutine Ham_Alloc_LR_Jz_Test
        allocate(ham_LR_Jz_Test::ham)
      end Subroutine Ham_Alloc_LR_Jz_Test

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_LR_Jz_Test_read_write_parameters.F90"

!--------------------------------------------------------------------
!> @brief
!> Sets the Hamiltonian
!--------------------------------------------------------------------
      Subroutine Ham_Set

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif
          Implicit none

          integer                :: unit_info
          Character (len=64)     :: file_info

#ifdef MPI
          Integer        :: ierr, Isize, Irank, irank_g, isize_g, igroup
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif

          ! From dynamically generated file "Hamiltonian_LR_Jz_Test_read_write_parameters.F90"
          call read_parameters()

          If (L1 < 2) then
             Write(6,*) 'LR_Jz_Test: L1 must be >= 2 (L2 = 1 gives a chain)'
             Call Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          endif
          If (abs(ham_Jz) < 1.d-12) then
             Write(6,*) 'LR_Jz_Test: ham_Jz must be non-zero'
             Call Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          endif

          Ltrot     = nint(beta/dtau)
          Thtrot    = 0
          Projector = .false.
          Symm      = .false.
          N_SUN     = 1
          N_FL      = 2

          ! Setup the Bravais lattice
          Call  Ham_Latt

          ! Setup the hopping / single-particle part
          Call  Ham_Hop

          ! Setup the interaction.
          call Ham_V

#ifdef MPI
          If (Irank_g == 0) then
#endif
             File_info = "info"
#if defined(TEMPERING)
             write(File_info,'(A,I0,A)') "Temp_",igroup,"/info"
#endif

             Open(newunit=unit_info, file=file_info, status="unknown", position="append")
             Write(unit_info,*) '====================================='
             Write(unit_info,*) 'Model is      : LR_Jz_Test (test-only model)'
             Write(unit_info,*) 'Lattice is    : ', Lattice_type
             Write(unit_info,*) 'L1            : ', L1
             Write(unit_info,*) 'L2            : ', L2
             Write(unit_info,*) '# of orbitals : ', Ndim
             Write(unit_info,*) 'Symm. decomp  : ', Symm
             Write(unit_info,*) 'Finite temperture version'
             Write(unit_info,*) 'Beta          : ', Beta
             Write(unit_info,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
             Write(unit_info,*) 't             : ', Ham_T
             Write(unit_info,*) 'Ham_U         : ', Ham_U
             Write(unit_info,*) 'Ham_Jz        : ', Ham_Jz
             Write(unit_info,*) 'Ham_chem      : ', Ham_chem
             Close(unit_info)
#ifdef MPI
          Endif
#endif

        End Subroutine Ham_Set

!--------------------------------------------------------------------
!> @brief
!> Sets  the  Lattice
!--------------------------------------------------------------------
        Subroutine Ham_Latt

          Implicit none
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)

          If (Lattice_Type /=  "Square")  then
             Write(6,*) 'LR_Jz_Test is only defined for the square lattice'
             Call Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          Endif

          Latt_unit%Norb    = 1
          Latt_unit%N_coord = 2
          allocate(Latt_unit%Orb_pos_p(Latt_unit%Norb,2))
          Latt_unit%Orb_pos_p(1, :) = [0.d0, 0.d0]

          a1_p(1) =  1.0  ; a1_p(2) =  0.d0
          a2_p(1) =  0.0  ; a2_p(2) =  1.d0
          L1_p    =  dble(L1)*a1_p
          L2_p    =  dble(L2)*a2_p
          Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )

          Ndim = Latt%N*Latt_unit%Norb

        end Subroutine Ham_Latt
!--------------------------------------------------------------------
!> @brief
!> Sets  the Hopping (one dense operator, no checkerboard)
!--------------------------------------------------------------------
        Subroutine Ham_Hop

          Implicit none

          Integer :: nf , I, Ix, Iy
          allocate(Op_T(1,N_FL))
          do nf = 1,N_FL
             Call Op_make(Op_T(1,nf),Ndim)
             Do I = 1,Latt%N
                Ix = Latt%nnlist(I,1,0)
                Op_T(1,nf)%O(I,  Ix) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                Op_T(1,nf)%O(Ix, I ) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                If ( L2 > 1 ) then
                   Iy = Latt%nnlist(I,0,1)
                   Op_T(1,nf)%O(I,  Iy) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                   Op_T(1,nf)%O(Iy, I ) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                endif
                Op_T(1,nf)%O(I,  I ) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                Op_T(1,nf)%P(i) = i
             Enddo
             Op_T(1,nf)%g      = -Dtau
             Op_T(1,nf)%alpha  =  cmplx(0.d0,0.d0, kind(0.D0))
             Call Op_set(Op_T(1,nf))
          enddo

        end Subroutine Ham_Hop

!--------------------------------------------------------------------
!> @brief
!> Sets the interaction: Ndim on-site Mz vertices, then the Jz bond vertices (x-bond, y-bond per site).
!--------------------------------------------------------------------
        Subroutine Ham_V

          Use Predefined_Int
          Implicit none

          Integer :: I, Ix, Iy, nc, N_bonds

          N_bonds = Latt%N
          If (L2 > 1) N_bonds = 2*Latt%N
          Allocate(Op_V(Ndim + N_bonds, N_FL))

          nc = 0
          Do I = 1, Ndim
             nc = nc + 1
             Call Predefined_Int_U_MZ( Op_V(nc,1), Op_V(nc,2), I, DTAU, Ham_U )
          Enddo
          Do I = 1, Latt%N
             Ix = Latt%nnlist(I,1,0)
             nc = nc + 1
             Call Predefined_Int_Jz( Op_V(nc,1), Op_V(nc,2), I, Ix, DTAU, Ham_Jz )
             If ( L2 > 1 ) then
                Iy = Latt%nnlist(I,0,1)
                nc = nc + 1
                Call Predefined_Int_Jz( Op_V(nc,1), Op_V(nc,2), I, Iy, DTAU, Ham_Jz )
             Endif
          Enddo

        end Subroutine Ham_V

!--------------------------------------------------------------------
!> @brief
!> Specifiy the equal time and time displaced observables
!--------------------------------------------------------------------
        Subroutine  Alloc_obs(Ltau)

          Implicit none
          !>  Ltau=1 if time displaced correlations are considered.
          Integer, Intent(In) :: Ltau
          Integer    ::  i, N, Nt
          Character (len=64) ::  Filename
          Character (len=:), allocatable ::  Channel


          ! Scalar observables
          Allocate ( Obs_scal(4) )
          Do I = 1,Size(Obs_scal,1)
             select case (I)
             case (1)
                N = 1;   Filename ="Kin"
             case (2)
                N = 1;   Filename ="Pot"
             case (3)
                N = 1;   Filename ="Part"
             case (4)
                N = 1;   Filename ="Ener"
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
          enddo

          ! Equal time correlators
          Allocate ( Obs_eq(5) )
          Do I = 1,Size(Obs_eq,1)
             select case (I)
             case (1)
                Filename = "Green"
             case (2)
                Filename = "SpinZ"
             case (3)
                Filename = "SpinXY"
             case (4)
                Filename = "SpinT"
             case (5)
                Filename = "Den"
             case default
                Write(6,*) "Error in Alloc_obs"
             end select
             Nt = 1
             Channel = "--"
             Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
          enddo

          If (Ltau == 1) then
             ! time-displaced correlators
             Allocate ( Obs_tau(5) )
             Do I = 1,Size(Obs_tau,1)
                select case (I)
                case (1)
                   Channel = 'P' ; Filename = "Green"
                case (2)
                   Channel = 'PH'; Filename = "SpinZ"
                case (3)
                   Channel = 'PH'; Filename = "SpinXY"
                case (4)
                   Channel = 'PH'; Filename = "SpinT"
                case (5)
                   Channel = 'PH'; Filename = "Den"
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = Ltrot+1-2*Thtrot
                If(Projector) Channel = 'T0'
                Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             enddo
          endif

        End Subroutine Alloc_obs

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Computes equal time observables (Pot includes U sum n_up n_do and the Jz bond energy)
!> @details
!> @param [IN] Gr   Complex(:,:,:)
!> \verbatim
!>  Green function: Gr(I,J,nf) = <c_{I,nf } c^{dagger}_{J,nf } > on time slice ntau
!> \endverbatim
!> @param [IN] Phase   Complex
!> \verbatim
!>  Phase
!> \endverbatim
!> @param [IN] Ntau Integer
!> \verbatim
!>  Time slice
!> \endverbatim
!-------------------------------------------------------------------
        subroutine Obser(GR,Phase,Ntau, Mc_step_weight)

          Use Predefined_Obs

          Implicit none

          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

          !Local
          Complex (Kind=Kind(0.d0)), allocatable :: GRC(:,:,:)
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, ZP,ZS, ZZ, ZXY, ZDen
          Integer :: I,J, imj, nf,  Ix, Iy

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))

          ZS = ZS*Mc_step_weight

          allocate(GRC(Ndim,Ndim,N_FL))

          Do nf = 1,N_FL
             Do I = 1,Ndim
                Do J = 1,Ndim
                   GRC(I, J, nf) = -GR(J, I, nf)
                Enddo
                GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
             Enddo
          Enddo
          ! GRC(i,j,nf) = < c^{dagger}_{i,nf } c_{j,nf } >

          ! Compute scalar observables.
          Do I = 1,Size(Obs_scal,1)
             Obs_scal(I)%N         =  Obs_scal(I)%N + 1
             Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo


          Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
          Do I = 1,Latt%N
             Ix = Latt%nnlist(I,1,0)
             Zkin = Zkin  + GRC(I,Ix,1)  + GRC(Ix,I,1)  &
                  &       + GRC(I,Ix,2)  + GRC(Ix,I,2)
          Enddo
          If (L2 > 1) then
             Do I = 1,Latt%N
                Iy = Latt%nnlist(I,0,1)
                Zkin = Zkin + GRC(I,Iy,2)  + GRC(Iy,I,2)   &
                     &      + GRC(I,Iy,1)  + GRC(Iy,I,1)
             Enddo
          Endif
          Zkin = Zkin*cmplx(-Ham_T,0.d0,Kind(0.d0))
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *ZP* ZS


          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          Do I = 1,Ndim
             ZPot = ZPot + Grc(i,i,1) * Grc(i,i, 2)
          Enddo
          Zpot = Zpot*ham_U
          ! Jz sum_<ij> S^z_i S^z_j,  S^z = (n_up - n_do)/2  (Wick: connected + disconnected parts)
          Do I = 1,Latt%N
             Ix = Latt%nnlist(I,1,0)
             Zpot = Zpot + ham_Jz * 0.25d0 * ( GRC(I,Ix,1)*GR(I,Ix,1) + GRC(I,Ix,2)*GR(I,Ix,2) + &
                  &   (GRC(I,I,2) - GRC(I,I,1))*(GRC(Ix,Ix,2) - GRC(Ix,Ix,1)) )
             If (L2 > 1) then
                Iy = Latt%nnlist(I,0,1)
                Zpot = Zpot + ham_Jz * 0.25d0 * ( GRC(I,Iy,1)*GR(I,Iy,1) + GRC(I,Iy,2)*GR(I,Iy,2) + &
                     &   (GRC(I,I,2) - GRC(I,I,1))*(GRC(Iy,Iy,2) - GRC(Iy,Iy,1)) )
             Endif
          Enddo
          Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + Zpot * ZP*ZS


          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Do I = 1,Ndim
             Zrho = Zrho + Grc(i,i,1) +  Grc(i,i,2)
          enddo
          Obs_scal(3)%Obs_vec(1)  =    Obs_scal(3)%Obs_vec(1) + Zrho * ZP*ZS
          Obs_scal(4)%Obs_vec(1)  =    Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*ZP*ZS

          Do I = 1,Size(Obs_eq,1)
             Obs_eq(I)%N         =  Obs_eq(I)%N + 1
             Obs_eq(I)%Ave_sign  =  Obs_eq(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo

          Do I = 1,Latt%N
             Do J = 1,Latt%N
                imj  = latt%imj(I,J)
                ZXY  = GRC(I,J,1) * GR(I,J,2) +  GRC(I,J,2) * GR(I,J,1)
                ZZ   = GRC(I,J,1) * GR(I,J,1) +  GRC(I,J,2) * GR(I,J,2)    + &
                       (GRC(I,I,2) - GRC(I,I,1))*(GRC(J,J,2) - GRC(J,J,1))

                ZDen = (GRC(I,I,1) + GRC(I,I,2)) * (GRC(J,J,1) + GRC(J,J,2)) + &
                     &  GRC(I,J,1) * GR(I,J,1)   +  GRC(I,J,2) * GR(I,J,2)
                
                Obs_eq(1)%Obs_Latt(imj,1,1,1) =  Obs_eq(1)%Obs_Latt(imj,1,1,1) + (GRC(I,J,1) + GRC(I,J,2))*ZP*ZS
                Obs_eq(2)%Obs_Latt(imj,1,1,1) =  Obs_eq(2)%Obs_Latt(imj,1,1,1) +  ZZ  *ZP*ZS
                Obs_eq(3)%Obs_Latt(imj,1,1,1) =  Obs_eq(3)%Obs_Latt(imj,1,1,1) +  ZXY *ZP*ZS
                Obs_eq(4)%Obs_Latt(imj,1,1,1) =  Obs_eq(4)%Obs_Latt(imj,1,1,1) + (2.d0*ZXY + ZZ)*ZP*ZS/3.d0
                Obs_eq(5)%Obs_Latt(imj,1,1,1) =  Obs_eq(5)%Obs_Latt(imj,1,1,1) +  ZDen * ZP * ZS


             enddo
             Obs_eq(5)%Obs_Latt0(1) = Obs_eq(5)%Obs_Latt0(1) + (GRC(I,I,1) + GRC(I,I,2)) *  ZP*ZS
          enddo

          deallocate(GRC)
        end Subroutine Obser
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Computes time displaced  observables
!> @details
!> @param [IN] NT, Integer
!> \verbatim
!>  Imaginary time
!> \endverbatim
!> @param [IN] GT0, GTT, G00, GTT,  Complex(:,:,:)
!> \verbatim
!>  Green functions:
!>  GT0(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(0  )>
!>  G0T(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(tau)>
!>  G00(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(0  )>
!>  GTT(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(tau)>
!> \endverbatim
!> @param [IN] Phase   Complex
!> \verbatim
!>  Phase
!> \endverbatim
!-------------------------------------------------------------------
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE, Mc_step_weight )

          Use Predefined_Obs

          Implicit none

          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

          !Locals
          Complex (Kind=Kind(0.d0)) :: ZP, ZS, ZZ, ZXY, ZDEN
          Integer :: IMJ, I, J

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          ZS = ZS*Mc_step_weight

          If (NT == 0 ) then
             Do I = 1,Size(Obs_tau,1)
                Obs_tau(I)%N         =  Obs_tau(I)%N + 1
                Obs_tau(I)%Ave_sign  =  Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
             Enddo
          Endif
          Do I = 1,Latt%N
             Do J = 1,Latt%N
                imj  = latt%imj(I,J)

                ZZ   =      (GTT(I,I,1) -  GTT(I,I,2) ) * ( G00(J,J,1)  -  G00(J,J,2) )   &
                     &    -  G0T(J,I,1) * GT0(I,J,1)  -  G0T(J,I,2) * GT0(I,J,2)
                ZXY  =    -  G0T(J,I,1) * GT0(I,J,2)  -  G0T(J,I,2) * GT0(I,J,1)


                ZDen =   (cmplx(2.d0,0.d0,kind(0.d0)) -  GTT(I,I,1) - GTT(I,I,2) ) * &
                     &   (cmplx(2.d0,0.d0,kind(0.d0)) -  G00(J,J,1) - G00(J,J,2) )   &
                     &   -G0T(J,I,1) * GT0(I,J,1)  -  G0T(J,I,2) * GT0(I,J,2)

                Obs_tau(1)%Obs_Latt(imj,NT+1,1,1) =  Obs_tau(1)%Obs_Latt(imj,NT+1,1,1) + (GT0(I,J,1) + GT0(I,J,2))*ZP*ZS
                Obs_tau(2)%Obs_Latt(imj,NT+1,1,1) =  Obs_tau(2)%Obs_Latt(imj,NT+1,1,1) +  ZZ  *ZP*ZS
                Obs_tau(3)%Obs_Latt(imj,NT+1,1,1) =  Obs_tau(3)%Obs_Latt(imj,NT+1,1,1) +  ZXY *ZP*ZS
                Obs_tau(4)%Obs_Latt(imj,NT+1,1,1) =  Obs_tau(4)%Obs_Latt(imj,NT+1,1,1) + (2.d0*ZXY + ZZ)*ZP*ZS/3.d0
                Obs_tau(5)%Obs_Latt(imj,NT+1,1,1) =  Obs_tau(5)%Obs_Latt(imj,NT+1,1,1) +  ZDen * ZP * ZS


             enddo
             Obs_tau(5)%Obs_Latt0(1) = Obs_tau(5)%Obs_Latt0(1) + &
                  &                   (cmplx(2.d0,0.d0,kind(0.d0)) -  GTT(I,I,1) - GTT(I,I,2))  *  ZP*ZS
          enddo



        end Subroutine OBSERT
    end submodule ham_LR_Jz_Test_smod
