!  Copyright (C) 2016-2026 The ALF project
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
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
!
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de .
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.


Program MaxEnt_Wrapper
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> General wrapper for maxent.  Handles particle, partile-hole, particle-particle,
!> channels as well as zero temperature.  See documentation for details.
!>
!
!--------------------------------------------------------------------
       Use runtime_error_mod
       Use Natural_Constants, only: pi, Eps_small
       Use MaxEnt_stoch_mod
       Use MaxEnt_mod
       use MaxEnt_Wrapper_mod
       use iso_fortran_env, only: output_unit, error_unit
       use Files_mod, only: str_to_upper
       Use Random_Wrap, only: lcg
#ifdef MPI
       Use mpi
#endif
#ifdef _OPENMP
       use check_omp_num_threads_mod
#endif

       Implicit  None

       Real (Kind=Kind(0.d0)), Dimension(:)  , allocatable :: XQMC, XQMC_st, XTAU, Xtau_st, &
            &                                                 Alpha_tot, om_bf, alp_bf, xom, A
       Real (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: XCOV, XCOV_st
       Real (Kind=Kind(0.d0))                              :: X_moments(2), Xerr_moments(2), ChiSq
       Character (Len=64)                                  :: command, File1, File2, file_boot
       Complex (Kind=Kind(0.d0))                           :: Z
       Logical, parameter                                  :: Test =.false.
       

       !----
       Integer                :: iboot, N_boot = 1 ! Set default to 1 in order to be compatible with original parameter files
       Real (Kind=Kind(0.d0)), Allocatable :: U_cov(:,:), Sigma_cov(:), Z_sample(:)
       Real (Kind=Kind(0.d0))              :: u, v, welford_delta, welford_delta2
       logical                             :: initial_checkpoint, N_boot_checkpoint = .false.
       logical                             :: Do_bootstrap, Reference_run
       integer                             :: N_runs
       integer                             :: ia
       Real (Kind=Kind(0.d0))              :: A_err
       Character (len=64)                  :: FA_boot_Name
       ! Welford accumulators for one bootstrap sample, flattened:
       !   Buf((ia-1)*Ndis + nw)   = A(om_nw) of Aom_ps_<ia>, ia = 1..N_alpha-10
       !   Buf(Ndis*N_alpha_1 + n) = (F,A) at alpha_n (QFI_ph.dat, DIDV.dat), n = 1..N_alpha
       Real (Kind=Kind(0.d0)), Allocatable :: Buf(:), W_mean(:), W_M2(:)
       Integer                             :: N_buf, N_acc, N_my_runs, k_buf
       ! MPI: rank 0 does the reference run and accumulates (Welford, online) the samples that
       ! ranks 1..Isize-1 compute in their own directories MaxSAC_rank_<r>/ with own seeds.
       Integer                             :: Isize = 1, Irank = 0, ierr, k_seed, i_seed
       Integer, Allocatable                :: seed_vec(:)
       Logical                             :: Use_MPI
#ifdef MPI
       Integer                             :: mpi_stat(MPI_STATUS_SIZE)
#endif
      !  Character (len=64)     :: temp, temp2

       !TODO: 22.06.2026
       !5. Ask for MPI, numerical stability -> of interpolation table
       !
       !9. Finish analysis of data: a) Extract Kappa(T) -> Plot; b) extract \tau(T) -> Plot; c) Plot GF; d) Extract \sigma(T)
       !10. Analytical continuation for electrical current. What could be the reason, why in Superlattice, I get such an answer, and here this answer. Did I do something wrong? 
            !=> I make a sanity check for the old simulation, whether I can reproduce the results I had with superlattice.




       !TODO: MPI -> seeds for random number generators
       !TODO: more robust implementation is not to generate the aom files until at last.
       !The important change is: do not write Aom_ps_*, Best_fit, energies, moments, dump_*, Max_stoch_log during every bootstrap sample. That avoids MPI file collisions and removes a lot of unnecessary I/O.

       !TODO: old max stoch mod + bootstrap -> new mod -> fallback to old if no bootstrap. -> new does not write the files, but explicitly calculates the cumulative means .> and only in final run writes the files.

       !Todo Numerical accuracy will go down with to many ntau -> dynamically check numerical accuracy -> dynamical increase Ndis_table size
       !Todo: better interpolation for kernel table -> sinh (in maxent_stoch_mod)
       !TODO: fallback if kernel goes out of table -> use to see, whether interpolation is actually good
       !TODO: How to use default model + Bootstrapping? Default model for smooth curves. => Best practice would be to start with T1 -> make one bootstrap sample -> T2 -> next bootstrap sample so on and on ..., but this will be slow AF, because parameters have to be read etc.
       !TODO: Average chi^2
       !TODO: Default model.
       !----

       Integer                :: Ngamma, Ndis,  NBins, NSweeps, Nwarm, N_alpha, N_cov
       Integer                :: N_skip, N_rebin, N_Back, N_auto, Norb
       Real (Kind=Kind(0.d0)) :: OM_st, OM_en,  alpha_st, R, Tolerance
       Logical                :: Checkpoint,  Stochastic, Default_model_exists, Particle_channel_PH
       Character (Len=:), allocatable :: Channel
       Character (Len=1)      :: Char, Char1
       Character (len=64)     :: str_temp, FA_Name
       ! Space  for classic MaxEnt
       Real (Kind=Kind(0.d0)), allocatable ::  Xker_classic(:,:),  A_classic(:),  Default(:)

       Integer                :: nt, nt1, io_error, n,nw, nwp, ntau, N_alpha_1, i,  nbin_qmc
       Integer                :: ntau_st, ntau_en, Ntau_old
       Real (Kind=Kind(0.d0)) :: dtau, xmom1, x,x1,x2, tau, omp, om, Beta,err, delta, Dom
       Real (Kind=Kind(0.d0)), parameter :: Zero = Eps_small
       Real (Kind=Kind(0.d0)) :: Alpha_classic_st=100000.d0
       Integer ::  N_BZ_Zones     =  1 
       Logical ::  Extended_Zone = .false.

       NAMELIST /VAR_Max_Stoch/ Ngamma, Ndis,  NBins, NSweeps, Nwarm, N_alpha, N_boot, N_boot_checkpoint, &
            &                   OM_st, OM_en,  alpha_st, R,  Checkpoint, Tolerance, &
            &                   Stochastic

       NAMELIST /VAR_errors/    N_skip, N_rebin, N_cov,  N_Back, N_auto,  N_BZ_Zones,  Extended_Zone

#ifdef MPI
       Call MPI_Init(ierr)
       Call MPI_Comm_size(MPI_COMM_WORLD, Isize, ierr)
       Call MPI_Comm_rank(MPI_COMM_WORLD, Irank, ierr)
#endif
#ifdef _OPENMP
       call check_omp_num_threads()
#endif

       Stochastic = .true. !  This is  the  default
       open(unit=30,file='parameters',status='old',action='read', iostat=io_error)
       if (io_error.eq.0) then
          READ(30,NML=VAR_errors)
          READ(30,NML=VAR_Max_Stoch)
       else
          write(error_unit,*) 'No file parameters '
          CALL Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
       endif
       close(30)

       ! MPI is only used for the bootstrap samples. Otherwise rank 0 does the usual single run.
       Use_MPI = Isize > 1 .and. Stochastic .and. N_boot > 1
       If (Isize > 1 .and. .not. Use_MPI .and. Irank /= 0) then
#ifdef MPI
          Call MPI_Finalize(ierr)
#endif
          Stop
       endif
       If (Use_MPI) then
          ! Own random number stream for every rank (otherwise all ranks draw the same samples)
          Call Random_seed(size=k_seed)
          Allocate(seed_vec(k_seed))
          i_seed = 8752143 + 104729*Irank
          Do i = 1, k_seed
             x = lcg(i_seed)
             seed_vec(i) = i_seed
          Enddo
          Call Random_seed(put=seed_vec)
          ! Own working directory for the files of MaxEnt_stoch on the worker ranks
          If (Irank /= 0) then
             write(MaxEnt_stoch_dir,'("MaxSAC_rank_",I0,"/")') Irank
             Call EXECUTE_COMMAND_LINE("mkdir -p " // trim(MaxEnt_stoch_dir))
          endif
       endif
       
       N = 20
       Call Set_Ra_ba(N)
       
       INQUIRE(FILE="Default", EXIST=Default_model_exists)

       open (unit=10,File="g_dat", status="unknown")
       read(10,*)  ntau, nbin_qmc, Beta, Norb,  str_temp
       Channel  = trim(str_temp)
       Allocate ( XCOV(NTAU,NTAU), XQMC(NTAU),XTAU(NTAU) )
       XCOV  = 0.d0
       Do nt = 1,NTAU
          read(10,*)  xtau(nt), xqmc(nt), err
          xcov(nt,nt) = err*err
       Enddo
       if (N_cov.eq.1) then
          do nt = 1,ntau
             do nt1 = 1,ntau
                read(10,*) xcov(nt,nt1)
             enddo
          enddo
       endif
       close(10)

       dtau = Xtau(2) - Xtau(1)
       11 format(A20, ': ', A)
       12 format(A20, ': ', I10)
       13 format(A20, ': ', *(E25.17E3))
       14 format(A20, ': ', (L1))
       15 format((E25.17E3))

       If (Stochastic) then 
          Open(unit=50,File=trim(MaxEnt_stoch_dir)//'Info_MaxEnt',Status="unknown")
       else
          Open(unit=50,File='Info_MaxEnt_cl',Status="unknown")
       endif
       write(50,11) 'Channel', Channel
       If (str_to_upper(Channel) == "PH" .or. str_to_upper(Channel) == "P_PH" .or. str_to_upper(Channel) ==  'PH_C')  then
          Write(50,"(A72)")  'Om_start is set to zero. PH  and  P_PH and PH_C channels corresponds to symmetric data'
          Om_st = 0.d0
       endif
       Write(50, 12) "Covariance", N_cov
       Write(50, 13) "Om_st", Om_st
       Write(50, 13) "Om_en", Om_en
       Write(50, 13) "Delta Om",  (Om_en - Om_st)/real(Ndis,kind(0.d0))
       Write(50, 14) "Default model exists",Default_model_exists 
       If (Stochastic) then
         Write(50, 14)  'Checkpoint' , Checkpoint
         Write(50, 12) "Bins",    NBins
         Write(50, 12) "Sweeps",  NSweeps
         Write(50, 12) "Warm",    Nwarm
         If (N_alpha <= 10 ) then
            Write(error_unit,*) 'Not enough temperatures: N_alpha has to be bigger than 10'
            CALL Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
         Endif
         Write(50,'(A54)')  "Tempertaure  set:  1/T(n=1..N_alpha) = Alpha_st*R**(n)"
         Write(50,12) "N_Alpha", N_alpha
         Write(50,13) "Alpha_st",alpha_st
         Write(50,13) "R", R
       endif
       Ntau_st = 1
       Ntau_en = Ntau
       Select Case (str_to_upper(Channel))
       Case ("PH", "PH_C")
          xmom1 = pi * xqmc(1)
       Case ("PP")
          xmom1 =  pi * (xqmc(1) + xqmc(ntau) )
       Case ("P")
          xmom1 =  pi * ( xqmc(1) + xqmc(ntau) )
          !  Remove the tau = beta point from the data since it is  correlated
          !  due to the sum rule with  the tau=0 data point. Also if the tau = 0
          !  data point has no fluctations (due to particle-hole symmetry for instance)
          !  it will be removed.
          Ntau_en = Ntau - 1
          Ntau_st = 1
          if ( xcov(1,1) < zero )  ntau_st = 2
       Case ("P_PH")
          xmom1 =  pi*xqmc(1)
          if ( xcov(1,1) < zero )  ntau_st = 2
          Particle_channel_PH = .true.
       Case ("T0")
          xmom1 =  pi*xqmc(1)
          Ntau_st = 1
          if ( xcov(1,1) < zero )  ntau_st = 2
       Case default
          Write(error_unit,*) "Channel '" // Channel // "' not yet implemented"
          CALL Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
       end Select
       Ntau_old = Ntau
       Call Rescale ( XCOV, XQMC,XTAU, Ntau_st, Ntau_en, Tolerance, NTAU)
       Write(50,"(A32, I4,A4, I4)") trim('Data has been rescaled from Ntau'), NTAU_old,trim("to"), Ntau
       If ( Ntau <= 4 ) then
          write(error_unit,*) 'Not enough data!'
          CALL Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
       Endif
       If (  nbin_qmc > 2*Ntau .and. N_cov == 0  )   Write(50,"(A72)") 'Consider using the covariance. You seem to have enough bins'
       If (  nbin_qmc < 2*Ntau .and. N_cov == 1  )   Write(50,"(A72)") 'Not enough bins for a reliable estimate of the covariance '

       !Store
       Allocate ( XCOV_st(NTAU,NTAU), XQMC_st(NTAU),XTAU_st(NTAU) )
       XCOV_st = XCOV
       XQMC_st = XQMC
       XTAU_st = XTAU
       Allocate (A_classic(Ndis), Default(Ndis), XKer_classic(size(Xqmc,1),Ndis))
       If (Default_model_exists) then
          Open(Unit=10,file="Default",status="unknown") 
          read (10,*) Char 
          rewind(10) 
          If (Char == "X" )   then
             do nw = 1,Ndis
               read(10,*) Char1,X, X1, Default(nw)
             enddo
          else 
             do nw = 1,Ndis
               read(10,*) X,Default(nw)
             enddo
          endif
          close(10)
       endif
       Call Set_default(Default,beta,Channel, OM_st, Om_en, xmom1,Default_model_exists,Stochastic)
      
       Allocate (Alpha_tot(N_alpha) )
       do nt = 1,N_alpha
          alpha_tot(nt) = alpha_st*(R**(nt-1))
       enddo
       write(50,13) "First Moment",  Xmom1
       write(50,13) "Beta", Beta

       !------Injection--------
       !Bootstrap
       N_alpha_1 = N_alpha - 10
       Allocate (xom(Ndis), A(Ndis))
       Allocate(U_cov(Ntau,Ntau), Sigma_cov(Ntau), Z_sample(Ntau)) !Allocate arrays for Gaussian resampling
       N_buf = Ndis*N_alpha_1 + N_alpha
       Allocate(Buf(N_buf), W_mean(N_buf), W_M2(N_buf)) !Allocate arrays for Welford online algorithm
       Buf = 0.d0; W_mean = 0.d0; W_M2 = 0.d0; N_acc = 0; N_my_runs = 0; A_err = 0.d0

       ! Bootstrap only for stochastic MaxEnt with N_boot > 1. The bootstrap samples are followed
       ! by one reference run on the unperturbed data, which provides energies, Best_fit, Data_out,
       ! Aom_ps_* and Green. N_boot <= 1 reproduces the original single run on the original data.
       Do_bootstrap = Stochastic .and. N_boot > 1
       If (Do_bootstrap) then
          N_runs = N_boot + 1
       else
          N_runs = 1
       endif
       If (.not. Stochastic .and. N_boot > 1) &
            & Write(50,"(A)") 'N_boot > 1 ignored: bootstrap is only implemented for stochastic MaxEnt'
       If (Do_bootstrap) Write(50,12) "N_boot", N_boot
       If (Use_MPI)      Write(50,12) "MPI ranks", Isize

       If (Do_bootstrap) Call Diag(XCOV_st, U_cov, Sigma_cov)  !Considering the usual sizes of Ntau this will mostlikely not be a bottlneck

       If (Do_bootstrap) then
       Do nt1 = 1, Ntau

         !! ---------------------------------------------------------
         !! SANITY CHECK: Ensure covariance is positive semi-definite
         !! We tolerate tiny negative values from numerical round-off
         !! (relative to the largest eigenvalue) and set them to zero,
         !! but crash if the matrix is fundamentally broken.
         !! ---------------------------------------------------------
           If (Sigma_cov(nt1) < -1.d-10*maxval(abs(Sigma_cov))) then
               write(error_unit,*) 'FATAL ERROR: Resampled Covariance matrix has a significantly negative eigenvalue!'
               write(error_unit,*) 'Eigenvalue index: ', nt1, ' Value: ', Sigma_cov(nt1)
               write(error_unit,*) 'Your QMC data may be corrupted, or you have too few bins.'
               CALL Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
           Endif
           Sigma_cov(nt1) = max(Sigma_cov(nt1), 0.d0)
       End Do
       End If

       initial_checkpoint = Checkpoint
       ! MaxEnt_stoch resets the random number generator to a fixed seed by default. Since the
       ! number of random draws per call does not depend on the data, every bootstrap sample
       ! would then receive identical Gaussian noise. Continue the global stream instead.
       If (N_boot > 1) MaxEnt_stoch_reset_seed = .false.

      ! 1. Box-Mueller Covariance Resampling
       Do iboot = 1, N_runs
            Reference_run = (iboot == N_runs)
            If (Use_MPI) then
               ! Rank 0: only the reference run (first, so that it can receive afterwards).
               ! Ranks 1..Isize-1: the bootstrap samples, round robin.
               If (Reference_run) then
                  If (Irank /= 0) cycle
               else
                  If (Irank /= 1 + mod(iboot-1, Isize-1)) cycle
               endif
            endif
            N_my_runs = N_my_runs + 1
            XQMC = XQMC_st
            XCOV = XCOV_st
            If (.not. Reference_run) then
            Do nt = 1, Ntau !For each tau point we get a different noise
               Do! U \in (0,1]
                   u = ranf_wrap()
                   If (u > tiny(0.d0)) Exit !Tiny, in order to not slow down the code, if a number smaller than 10^(-308) was sampled
               End Do

               Do ! v \in [0,1)
                  v = ranf_wrap()
                  If (v< 1.d0) Exit
               End Do 
               Z_sample(nt) = sqrt(-2.d0 * log(u)) * cos(2.d0 * pi * v)  !Box mueller of the two uniformly sampled u, v to a Gaussian N(\mu =0, \sigma=1)
            End Do

            !Apply noise to data
            Do nt = 1, Ntau
               Do nt1 = 1, Ntau
                  XQMC(nt) = XQMC(nt) + U_cov(nt, nt1) * sqrt(Sigma_cov(nt1)) * Z_sample(nt1) ! Sigma_cov >= 0 ensured above
               End Do
            End Do
            End If
         
            !2. Warmup & File Caching Logic
            !   A dump of a previous run (Checkpoint at iboot=1, or the previous bootstrap
            !   sample if N_boot_checkpoint) is only used as starting configuration: the data
            !   differ, so accumulated histograms must not be reused and the warmup is kept.
            If (Do_bootstrap) then
               If (N_my_runs == 1) then
                  If (.not. Initial_Checkpoint) &
                       & Call EXECUTE_COMMAND_LINE("rm -f " // trim(MaxEnt_stoch_dir) // "dump*")
               Else If (.not. N_boot_checkpoint) then
                  Call EXECUTE_COMMAND_LINE("rm -f " // trim(MaxEnt_stoch_dir) // "dump*")
               End If
               MaxEnt_stoch_warm_start_only = .true.
               Checkpoint = Initial_Checkpoint .or. N_boot_checkpoint
            End If
       !------Injection_END----

       FA_Name = ""   ! Set below if the channel writes an (F,A) file
       Select Case (str_to_upper(Channel))
       Case ("PH")
          If  (Stochastic)  then
             FA_Name = "QFI_ph.dat"
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_ph, Back_Trans_ph, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm, F_QFI_ph, &
                  &            Filename_F=FA_Name, Default_provided=Default)
             ! Beware: Xqmc and cov are modified in the MaxEnt_stoch call.
          else
             Call Set_Ker_classic(Xker_ph,Xker_classic,Om_st,Om_en,beta,xtau_st)
             Call  MaxEnt( XQMC, XCOV, A_classic, XKER_classic, Alpha_classic_st, CHISQ ,DEFAULT)
          endif
       Case ("PH_C")
          If  (Stochastic)  then
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_ph, Back_trans_pp, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm, F, Default_provided=Default)
             ! Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_ph_c, Back_Trans_ph_c, Beta, &
                  ! &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm, F_QFI_ph_c, Default)
             ! Beware: Xqmc and cov are modified in the MaxEnt_stoch call.
          else
             ! Call Set_Ker_classic(Xker_ph_c,Xker_classic,Om_st,Om_en,beta,xtau_st)
             Call Set_Ker_classic(XKER_ph,Xker_classic,Om_st,Om_en,beta,xtau_st)
             Call  MaxEnt( XQMC, XCOV, A_classic, XKER_classic, Alpha_classic_st, CHISQ ,DEFAULT)
          endif       
       Case ("PP")
          If  (Stochastic) then
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_pp, Back_Trans_pp, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm, F , Default_provided=Default)
             ! Beware: Xqmc and cov are modified in the MaxEnt_stoch call.
          else
             Call Set_Ker_classic(Xker_pp,Xker_classic,Om_st,Om_en,beta,xtau_st)
             Call  MaxEnt( XQMC, XCOV, A_classic, XKER_classic, Alpha_classic_st, CHISQ ,DEFAULT)
          endif
       Case ("P")
          If  (Stochastic)  then
             FA_Name = "DIDV.dat"
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_p, Back_Trans_p, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm , F_DIDV, & 
                  &            Filename_F=FA_Name, Default_provided=Default)
             ! Beware: Xqmc and cov are modified in the MaxEnt_stoch call.
          else  ! Classic
             Call Set_Ker_classic(Xker_p,Xker_classic,Om_st,Om_en,beta,xtau_st)
             Call  MaxEnt( XQMC, XCOV, A_classic, XKER_classic, Alpha_classic_st, CHISQ ,DEFAULT)
          endif  
       Case ("P_PH")
          If  (Stochastic)  then
             FA_Name = "DIDV.dat"
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_p_ph, Back_Trans_p, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm ,F_DIDV_PH, & 
                  &            Filename_F=FA_Name, Default_provided=Default)
          else  ! Classic
             Call Set_Ker_classic(Xker_p_ph,Xker_classic,Om_st,Om_en,beta,xtau_st)
             Call  MaxEnt( XQMC, XCOV, A_classic, XKER_classic, Alpha_classic_st, CHISQ ,DEFAULT)
          endif  
       Case ("T0")
          If (Stochastic)  then
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_T0, Back_Trans_T0, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm,F,Default_provided=Default)
             ! Beware: Xqmc and cov are modified in the MaxEnt_stoch call.
          else
             Call Set_Ker_classic(Xker_T0,Xker_classic,Om_st,Om_en,beta,xtau_st)
             Call  MaxEnt( XQMC, XCOV, A_classic, XKER_classic, Alpha_classic_st, CHISQ ,DEFAULT)
          endif
       Case default
          Write(error_unit,*) "Channel '" // Channel // "' not yet implemented"
          CALL Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
       end Select

       !Injection Part 2: Welford Online Algorithm (No data is written yet!)
       !TODO: Implement welford in MaxEnt_stpcj directly. -> Write files only after last N_boot
       !TODO: Also bootstrap A_om
       If (Do_bootstrap .and. .not. Reference_run) then
         Call Read_sample(Buf)
         If (Use_MPI) then
#ifdef MPI
            Call MPI_Send(Buf, N_buf, MPI_DOUBLE_PRECISION, 0, 77, MPI_COMM_WORLD, ierr)
#endif
         else
            Call Welford_add(Buf)
         endif
         !INJECTION_PART_2_END
      End If
   End Do

#ifdef MPI
       If (Use_MPI) then
          If (Irank == 0) then
             ! Online Welford update in the order in which the samples arrive.
             Do n = 1, N_boot
                Call MPI_Recv(Buf, N_buf, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 77, &
                     &        MPI_COMM_WORLD, mpi_stat, ierr)
                Call Welford_add(Buf)
             Enddo
          else
             ! Workers are done.
             If (.not. Checkpoint) &
                  & Call EXECUTE_COMMAND_LINE("rm -f " // trim(MaxEnt_stoch_dir) // "dump*")
             close(50)
             Call MPI_Finalize(ierr)
             Stop
          endif
       endif
#endif


         !INJECTION 3
         !Calculate Bootstrap error from variance.
         If (Do_bootstrap) Then
         ! Frequency grid from the reference run (current directory)
         Open(Unit=66, File="Aom_ps_1", status="old")
         Do nw = 1, Ndis
            read(66,*) xom(nw)
         Enddo
         close(66)
         Do ia = 1, N_alpha_1
            write(file_boot, '(A,"_",I0)') "Aom_ps_boot", ia
            open (Unit=12, File=file_boot, Status="unknown", action="write") 
            ! Write(12,"(A14,2x,A16,2x,A16)") "# omega", "A_mean", "A_error"
            Do nw = 1, Ndis
               k_buf = (ia-1)*Ndis + nw
               A_err = sqrt(W_M2(k_buf) / dble(max(1, N_acc - 1)))
               Write(12,"(3(1X,E25.17E3))") xom(nw), W_mean(k_buf), A_err
            End Do
            Close(12)
         End Do
         ! (F,A) observable: <name>_boot.dat with alpha, bootstrap mean, bootstrap error.
         If (len_trim(FA_Name) > 0) then
            i = index(FA_Name, ".dat", back=.true.)
            If (i > 1) then
               FA_boot_Name = FA_Name(1:i-1) // "_boot.dat"
            else
               FA_boot_Name = trim(FA_Name) // "_boot"
            endif
            open (Unit=12, File=trim(FA_boot_Name), Status="unknown", action="write")
            Do ia = 1, N_alpha
               k_buf = Ndis*N_alpha_1 + ia
               Write(12,"(3(1X,E25.17E3))") Alpha_tot(ia), W_mean(k_buf), &
                    & sqrt(W_M2(k_buf) / dble(max(1, N_acc - 1)))
            End Do
            Close(12)
         End If
         End If

         !INJECTION 3 END
      
         If  ( Stochastic )   then
         If ( .not.  Checkpoint  ) then
           Command = "rm dump*"
           Call EXECUTE_COMMAND_LINE(Command)
           Command = "ls"
           Call EXECUTE_COMMAND_LINE(Command)
         endif
          
          Open (Unit=10,File="energies",status="unknown")

          Do n = 1,N_alpha
             Read(10,*) X,X1,X2
          enddo
          Write(50,13) "Best Chisq", X1
          close(50)

          Open (Unit = 10,File="Best_fit", Status ="unknown")
          Allocate (om_bf(Ngamma), alp_bf(Ngamma) )
          DO i = 1, Ngamma
            read(10,*)  om_bf(i), alp_bf(i)
          Enddo
          close(10)

          Open (Unit = 11,File="Data_out", Status ="unknown")
          do nt = 1,Ntau
             X = 0.d0
             tau = xtau_st(nt)
             Select Case (str_to_upper(Channel))
                Case ("PH")
                   do i = 1,Ngamma
                      X = X + alp_bf(i)*Xker_ph(tau,om_bf(i), beta)
                   enddo
                Case ("PH_C")
                   do i = 1,Ngamma
                      ! X = X + alp_bf(i)*Xker_ph_c(tau,om_bf(i), beta)
                      X = X + alp_bf(i)*Xker_ph(tau,om_bf(i), beta)
                   enddo
                Case ("PP")
                   do i = 1,Ngamma
                      X = X + alp_bf(i)*Xker_pp(tau,om_bf(i), beta)
                   enddo
                Case ("P_PH")
                   do i = 1,Ngamma
                      X = X + alp_bf(i)*Xker_p_ph(tau,om_bf(i), beta)
                   enddo
                Case ("P")
                   do i = 1,Ngamma
                      X = X + alp_bf(i)*Xker_p(tau,om_bf(i), beta)
                   enddo
                Case ("T0")
                   do i = 1,Ngamma
                      X = X + alp_bf(i)*Xker_T0(tau,om_bf(i), beta)
                   enddo
                Case default
                   Write(error_unit,*) "Channel '" // Channel // "' not yet implemented"
                   CALL Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
             end Select
             Write(11,"(E25.17E3,2x,E25.17E3,2x,E25.17E3,2x,E25.17E3)")  xtau_st(nt), xqmc_st(nt),  sqrt(xcov_st(nt,nt)), xmom1*X
          enddo
          close(11)
          N_alpha_1 = N_alpha - 10
          File1 ="Aom_ps"
          write(file2, '(A,"_",I0)') trim(File1), N_alpha_1
          Open(Unit=66,File=file2,status="unknown")
          do nw = 1,Ndis
             read(66,*) xom(nw), A(nw), x, x1, x2
          enddo
          close(66)
          Dom = xom(2) - xom(1)
          Open (Unit=43,File="Green", Status="unknown", action="write")
       else 
          DOM  =  (OM_En -  OM_St)/dble(Ndis)
          A =  A_classic/Dom
          do  nw  = 1,Ndis
             xom(nw) =  OM_St +  dble(nw-1)*dom
          enddo
          Open (Unit = 11,File="Data_out_cl", Status ="unknown")
          Do Nt = 1,Ntau
             X = 0.d0
             Do nw = 1,Ndis
                X = X  +  A_classic(nw)*Xker_classic(nt,nw)
             enddo
             Write(11,"(E25.17E3,2x,E25.17E3,2x,E25.17E3,2x,E25.17E3)")  xtau_st(nt), xqmc_st(nt),  sqrt(xcov_st(nt,nt)), X
          enddo
          close(11)
          Write(50,13) "Final value of  alpha ", 1.d0/Alpha_classic_st
          Write(50,13) "CHISQ" , CHISQ
          close(50)
          Select Case (str_to_upper(Channel))
             Case ("PH")
                do  nw  = 1,Ndis
                   A(nw) =  Back_trans_ph(A(nw), xom(nw), beta)
                enddo
             Case ("PH_C")
                do  nw  = 1,Ndis
                   ! A(nw) =  Back_trans_ph_c(A(nw), xom(nw), beta)
                   A(nw) =  Back_trans_pp(A(nw), xom(nw), beta)
                enddo
             Case ("PP")
                do  nw  = 1,Ndis
                   A(nw) =  Back_trans_pp(A(nw), xom(nw), beta)
                enddo
             Case ("P")
                do  nw  = 1,Ndis
                   A(nw) =  Back_trans_p(A(nw), xom(nw), beta)
                enddo
             Case ("P_PH")
                do  nw  = 1,Ndis
                   A(nw) =  Back_trans_p(A(nw), xom(nw), beta)
                enddo
             Case ("T0")
                do  nw  = 1,Ndis
                   A(nw) =  Back_trans_T0(A(nw), xom(nw), beta)
                enddo
             Case default
                Write(error_unit,*) "Channel '" // Channel // "' not yet implemented"
                CALL Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
             end Select
          Open (Unit=43,File="Green_cl", Status="unknown", action="write")
       endif 

       ! Compute the real frequency Green function.
       delta = Dom
       x  = 0.d0
       x1 = 0.d0
       x2 = 0.d0
       do nw = 1,Ndis
          Z = cmplx(0.d0,0.d0,Kind(0.d0))
          om = xom(nw)
          do nwp = 1,Ndis
             omp = xom(nwp)
             If  (((str_to_upper(Channel) == "P_PH") .or. (str_to_upper(Channel) == "PH_C")) .and.  omp > 0.00001d0  )  then 
               ! In this case A(om)  = A (-om) and om > 0
               Z = Z + A(nwp)/cmplx(  om -  omp, delta, kind(0.d0)) &
                   & + A(nwp)/cmplx(  om +  omp, delta, kind(0.d0)) 
            elseif ( str_to_upper(Channel) == "PH" .and.  omp > 0.00001d0  ) then 
               ! In this case S(q,-om)=  exp(-beta*om) S(q,om)   and om > 0 
               Z = Z + A(nwp)/cmplx(  om -  omp, delta, kind(0.d0)) &
                   & + A(nwp)*exp(-beta*om)/cmplx(  om +  omp, delta, kind(0.d0)) 
            else
               ! No symmetry is used
               Z = Z + A(nwp)/cmplx( om -  omp, delta, kind(0.d0))
            endif
          enddo
          Z = -Z * dom/pi
          If (Test)   then 
            x  = x  + Aimag(Z)
            x1 = x1 + om*Aimag(Z)
            x2 = x2 + om*om*Aimag(Z)
            write(43,"('X',2x,E25.17E3,2x,E25.17E3,2x,E25.17E3,2x,E25.17E3,2x,E25.17E3,2x,E25.17E3)")  & 
               & xom(nw), dble(Z), Aimag(Z),  X*dom, x1*dom, x2*dom
          else
            write(43,"('X',2x,E25.17E3,2x,E25.17E3,2x,E25.17E3)")  & 
               & xom(nw), dble(Z), Aimag(Z)
          endif   
       enddo
       close(43)

       call clean_Set_Ra_ba()
#ifdef MPI
       Call MPI_Finalize(ierr)
#endif

     contains

       !> Reads the result of one bootstrap sample (files in MaxEnt_stoch_dir) into B.
       Subroutine Read_sample(B)
         Implicit none
         Real (Kind=Kind(0.d0)), Intent(Out) :: B(:)
         Integer :: ia_l, nw_l
         Real (Kind=Kind(0.d0)) :: y0, y2, y3, y4
         Character (len=320) :: fname
         B = 0.d0
         Do ia_l = 1, N_alpha_1
            write(fname, '(A,"_",I0)') trim(MaxEnt_stoch_dir) // "Aom_ps", ia_l
            Open(Unit=66, File=trim(fname), status="old")
            Do nw_l = 1, Ndis
               read(66, *) y0, B((ia_l-1)*Ndis + nw_l), y2, y3, y4
            Enddo
            close(66)
         Enddo
         ! (F,A) observable: one line (alpha, F, MC error) per alpha
         If (len_trim(FA_Name) > 0) then
            Open(Unit=67, File=trim(MaxEnt_stoch_dir) // trim(FA_Name), status="old")
            Do ia_l = 1, N_alpha
               read(67, *) y0, B(Ndis*N_alpha_1 + ia_l), y2
            Enddo
            close(67)
         endif
       end Subroutine Read_sample

       !> Welford online update of mean and M2 with one more sample B.
       !> mean_n = mean_{n-1} + (x_n - mean_{n-1})/n
       !> M2_n   = M2_{n-1} + (x_n - mean_{n-1}) (x_n - mean_n);  variance = M2_n/(n-1)
       Subroutine Welford_add(B)
         Implicit none
         Real (Kind=Kind(0.d0)), Intent(In) :: B(:)
         Integer :: k
         N_acc = N_acc + 1
         Do k = 1, N_buf
            welford_delta  = B(k) - W_mean(k)
            W_mean(k)      = W_mean(k) + welford_delta/dble(N_acc)
            welford_delta2 = B(k) - W_mean(k)
            W_M2(k)        = W_M2(k) + welford_delta*welford_delta2
         Enddo
       end Subroutine Welford_add

     end Program MaxEnt_Wrapper
