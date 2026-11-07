s=open("main.F90").read()
old="                          CALL ham%Obser( GR_Tilde, PHASE, Ntau1, Mc_step_weight )"
new="""                          if (Ntau1 >= 17 .and. Ntau1 <= 24 .and. Ndim == 16) write(*,'(A,I4,4F14.8)') 'DBG up symm raw1 raw2 tilde1 tilde2', Ntau1, real(GR(6,6,1)), real(GR(6,6,2)), real(GR_Tilde(6,6,1)), real(GR_Tilde(6,6,2))
                          CALL ham%Obser( GR_Tilde, PHASE, Ntau1, Mc_step_weight )"""
assert old in s
s=s.replace(old,new,1)
old2="                          CALL ham%Obser( GR, PHASE, Ntau1, Mc_step_weight  )"
new2="""                          if (Ntau1 >= 17 .and. Ntau1 <= 24 .and. Ndim == 16) write(*,'(A,I4,2F14.8)') 'DBG up nosymm raw1 raw2', Ntau1, real(GR(6,6,1)), real(GR(6,6,2))
                          CALL ham%Obser( GR, PHASE, Ntau1, Mc_step_weight  )"""
assert old2 in s
s=s.replace(old2,new2,1)
old3="                          CALL ham%Obser( GR, PHASE, Ntau1, Mc_step_weight  )"
new3="""                          if (Ntau1 == 17 .and. Ndim == 16) write(*,'(A,I4,8F8.3)') 'DBGg Ltrot', Ltrot, real(Op_T(2,1)%g_t(17:24))
                          if (Ntau1 == 17 .and. Ndim == 16) then
                             do ii = 1, 16
                                do jj = 1, 16
                                   write(*,'(A,2I3,2F16.10)') 'DBGO', ii, jj, real(Op_T(1,1)%O(ii,jj)), aimag(Op_T(1,1)%O(ii,jj))
                                enddo
                             enddo
                          endif
                          CALL ham%Obser( GR, PHASE, Ntau1, Mc_step_weight  )"""
s=s.replace(old3,new3,1)
s=s.replace("        Integer :: NTAU, NTAU1, nt_obs\n","        Integer :: NTAU, NTAU1, nt_obs, ii, jj\n",1)
open("main.F90","w").write(s)
