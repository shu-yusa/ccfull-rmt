!======================================================================!
!     Title    : scat_noncoll.f90                                !
!                                                                      !
!     *** Structure of This Program ***                                !
!                                                                      !
!     module    global_constant                                        !
!     module    input_data                                             !
!     module    potentials                                             !
!     module    relative_potential -> potentials (inheritance)         !
!     module    coupling_matrix    -> potentials                       !
!     module    coupled_channels                                       !
!     module    calc_profile                                           !
!                                                                      !
!     program   main                                                   !
!       subroutine  comm                                               !
!       subroutine  output                                             !
!======================================================================!
  program main
    use, intrinsic :: iso_fortran_env
    use global_constant, only : PI
    use input_data
    use relative_potential
    use coupling_matrix
    use coupled_channels
    use calc_profile
    use MPI
!   use my_fitting
    use print_array
    use angular_coup
    implicit none
    logical :: b
    integer :: i, n, m, ios, getcwd, chdir, idum, nr
    integer :: ierr, nprocs, myrank, ist, ied, proc_len, rem   ! MPI
    integer, allocatable, dimension(:) :: displs, recvcnt
    real(8), allocatable, dimension(:) :: Ea, Esig, angl
    real(8), allocatable, dimension(:) :: sig_iel, sig_fus, spin
    real(8), allocatable, dimension(:) :: sig_fus_sum
    real(8), allocatable, dimension(:,:) :: sig_iel_n
    real(8), allocatable, dimension(:,:) :: dsig_iel, dsig_qel, dsig_R
    real(8), allocatable, dimension(:,:,:) :: dsig_qel_n
    real(8), allocatable, dimension(:,:,:) :: dsig_qel_n_sum
    real(8) :: w_l(2:10)
    complex(8), allocatable, dimension(:,:) :: fc
    complex(8), allocatable, dimension(:,:,:) :: fN
    character(len=500) :: dir, cwd, wd
    character(len=500) :: dir2='angular'
    character(len=500) :: dir3='angular_dist'
    character(len=500) :: dir4="Q_val_dist"
    character(len=20) :: outdir = 'results'
    character(len=8) :: proc_name
    character(len=100), parameter :: INPF="input_scat_noncoll"
    character(len=100) :: noncoll_inp
    character(len=3) :: rand
    type(inp) :: ip
    type(cc_scat) :: cc_calc
    type(rel_pot) :: vrel_pot
    type(coup_mat) :: cm
    type(profiler) :: prof
!   type(my_fit) :: mf
!   type(my_fit2) :: mf2
    integer :: I1, L
    real(8) :: e1, s, f, wg, r
    real(8) :: Delta, sig, alpha, rp

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(mpi_comm_world, nprocs, ierr)
    call MPI_COMM_RANK(mpi_comm_world, myrank, ierr)
    call MPI_GET_PROCESSOR_NAME(proc_name, proc_len, ierr)
    call ip%read_input(trim(INPF),noncoll_inp, dir)

    ios = getcwd(cwd)
    if (myrank == 0) then
      call check_directory(trim(cwd)//'/'//trim(outdir))
      call check_directory(trim(cwd)//'/'//trim(outdir)//'/'//dir)
      call check_directory(trim(cwd)//'/'//trim(outdir)//'/'//trim(dir)//'/'//dir2)
      call check_directory(trim(cwd)//'/'//trim(outdir)//'/'//trim(dir)//'/'//dir3)
      call check_directory(trim(cwd)//'/'//trim(outdir)//'/'//trim(dir)//'/'//dir4)
      call make_dirs_below(trim(cwd)//'/'//trim(outdir)//'/'//trim(dir), ip%Ncomp)
    end if
    call mpi_barrier(mpi_comm_world, ierr)

    idum = -123456789
    call vrel_pot%rel_pot_(ip)
    call cm%coup_mat_(ip, "input_rmt", trim(noncoll_inp), nprocs, myrank,idum)
    call cc_calc%cc_scat_(ip, vrel_pot, cm)
    call prof%profiler_(ip, vrel_pot, cm)

    allocate(Ea(ip%Egrid+1))
    allocate(Esig(ip%Egrid+1))
    allocate(angl(ip%nth))
    allocate(fN(ip%Nch,ip%nth,ip%Egrid+1))
    allocate(fc(ip%nth,ip%Egrid+1))
    allocate(dsig_iel(ip%nth,ip%Egrid+1))
    allocate(dsig_qel(ip%nth,ip%Egrid+1))
    allocate(dsig_qel_n(ip%Nch,ip%nth,ip%Egrid+1))
    allocate(dsig_R(ip%nth,ip%Egrid+1))
    allocate(spin(ip%Egrid+1))
    allocate(sig_fus(ip%Egrid+1))
    allocate(sig_iel_n(ip%Nch,ip%Egrid+1))
    allocate(sig_iel(ip%Egrid+1))

    allocate(sig_fus_sum(ip%Egrid+1))
    allocate(dsig_qel_n_sum(ip%Nch,ip%nth,ip%Egrid+1))

    allocate(displs(0:nprocs-1), recvcnt(0:nprocs-1))

    forall (i=1:ip%Egrid+1) Ea(i) = ip%Emin + dble(i-1) * ip%dE
    forall (i=1:ip%nth) angl(i) = ip%thmin + dble(i-1) * ip%dth
    angl = angl * PI / 180.0d0

    rem = mod(ip%Egrid+1, nprocs)
    if (rem == 0 .or. myrank < rem) then
      ist = (ip%Egrid / nprocs + 1) * myrank + 1
      ied = (ip%Egrid / nprocs + 1) * (myrank + 1)
    else
      ist = (ip%Egrid / nprocs) * myrank + 1 + rem 
      ied = (ip%Egrid / nprocs) * (myrank + 1) + rem
    end if

    b = cm%find_rmin(0)
    call vrel_pot%make_Vrel()
    if (myrank == 0) then
        ios = chdir(trim(cwd)//'/'//trim(outdir)//'/'//trim(dir))
      call prof%pot_prof()
      call prof%Reaction_prof("calc_info", trim(noncoll_inp))
      ios = chdir(trim(cwd))
    end if
!   rp = cm%rb
!   write(6,*) '  r    =', rp, ' fm'
!   write(6,*) '  V(r) =', Vn(ip,rp) + Vc(ip,rp), ' MeV'
!   stop

    sig_fus_sum = 0.0d0
    dsig_qel_n_sum = 0.0d0
    do nr=1, ip%Ncomp
      if (myrank == 0) then
        write(output_unit,*) nr, ' th calculation'
        write(output_unit,*) 'make Vcp...'
      end if
      call cm%noise_regenerate()
      call cm%make_Vcp()

      if (nprocs > 1) then
        call mpi_barrier(mpi_comm_world, ierr)
        call comm2()
!       write(output_unit,*) myrank,'@',proc_name,':comm2. finished '
      end if

      do i=ist, ied
        write(output_unit,'(1x,2i6,a,i6)') myrank,i," /",ied
        call cc_calc%cc_scattering(Ea(i), spin(i), sig_fus(i), sig_iel_n(:,i), sig_iel(i))
        call cc_calc%scat_amp_Coul(Ea(i), angl, fc(:,i))
        call cc_calc%scat_amp_nucl(Ea(i), angl, fN(:,:,i))
      end do
      write(output_unit,*) myrank, '@',proc_name,':finished'

      if (nprocs > 1) then
        call mpi_barrier(mpi_comm_world, ierr)
        call comm()
        write(output_unit,*) myrank,'@',proc_name,':comm. finished '
      end if

      if (myrank == 0) then
        Esig = Ea * sig_fus
        dsig_R = abs(fc) ** 2 * 10.0d0
        dsig_qel_n(1,:,:) = abs(fc + fN(1,:,:)) ** 2 * 10.0d0
        do n=2, ip%Nch
          dsig_qel_n(n,:,:) = abs(fN(n,:,:)) ** 2 * 10.0d0
        end do
        dsig_iel = sum(dsig_qel_n(2:ip%Nch,:,:),dim=1) 
        dsig_qel = dsig_iel + dsig_qel_n(1,:,:)

        sig_fus_sum = sig_fus_sum + sig_fus
        dsig_qel_n_sum = dsig_qel_n_sum + dsig_qel_n

        write(rand,'(i3)') nr
        wd = trim(dir)//'/rand_'//trim(adjustl(rand))
        ios = chdir(trim(cwd)//'/'//trim(outdir)//'/'//trim(dir))
        call output()
        call qel_corrected_angle()

        ! averaged cross sections
        sig_fus = sig_fus_sum / dble(nr)
        dsig_qel_n = dsig_qel_n_sum / dble(nr)

        Esig = Ea * sig_fus
        dsig_iel = sum(dsig_qel_n(2:ip%Nch,:,:),dim=1) 
        dsig_qel = dsig_iel + dsig_qel_n(1,:,:)
        ios = chdir(trim(cwd)//'/'//trim(dir))
        call output()
        call qel_corrected_angle()
        open(7,file="last_cacl")
          write(7,*) '# of comps =',nr
        close(7)

        ios = chdir(trim(cwd))

      end if
    end do

    call cc_calc%destruct_cc_scat()
    call vrel_pot%destruct_rel_pot()
    call cm%destruct_coup_mat()
    if (myrank == 0) then
      write(output_unit,*) "calculation finished."
      write(output_unit,*) 'directory : '//trim(outdir)//'/'//trim(dir)
    end if

    deallocate(Ea,  Esig, angl)
    deallocate(sig_fus, sig_iel, sig_iel_n, spin)
    deallocate(fc, fN)
    deallocate(dsig_iel, dsig_qel, dsig_qel_n, dsig_R)
    call MPI_FINALIZE(ierr)

    contains
!********************************************************************!
    subroutine comm2()
      implicit none
      integer :: ist, ied, rem
      integer :: n, m, k, step, nst, nlen
      real(8), pointer :: p
      real(8), allocatable :: tmp(:)


!     step = nint(ip%rmax / ip%dr) + 2
      step = nint(ip%rcut / ip%dr)
      rem = mod(step, nprocs-1)
      if (rem == 0 .or. myrank < rem) then
        ist = ((step-1)/ (nprocs-1)+1) * myrank + 1
        ied = ((step-1)/ (nprocs-1)+1) * (myrank + 1)
      else
        ist = ((step-1)/ (nprocs-1)) * myrank + 1 + rem 
        ied = ((step-1)/ (nprocs-1)) * (myrank + 1) + rem
      end if
      if (myrank == nprocs-1) then
        ist = step+1
        ied = ip%rgrid+2
      end if


      do m=0, nprocs-2
        if (rem == 0 .or. m < rem) then
          recvcnt(m) = (step-1) / (nprocs-1) + 1
        else
          recvcnt(m) = (step-1) / (nprocs-1)
        end if
        displs(m) = sum(recvcnt(:m-1))
      end do
      recvcnt(nprocs-1) = ip%rgrid+2 - step
      displs(nprocs-1) = sum(recvcnt(:nprocs-2))

!     if (myrank == 0) then
!       write(6,*) 'n=',n!,'nlen=',nlen
!       call print_vec(dble(recvcnt))
!       call print_vec(dble(displs))
!       write(6,*) 'rgrid+2=',ip%rgrid+2
!       write(6,*) '------'
!     end if
!     stop

      n = ied - ist + 1
      nlen = ip%Nch * (ip%Nch + 1) / 2
      allocate(tmp(n))
      nst = 1
      do k=1, nlen
        tmp = cm%Vcp_linear(k,ist:ied)
        call mpi_allgatherv(tmp,n,mpi_real8,           &
                        cm%Vcp_linear(k,:),recvcnt,displs,&
                        mpi_real8, mpi_comm_world,ierr)
      end do

      deallocate(tmp)


    end subroutine
!********************************************************************!
    subroutine comm()
      implicit none
      integer, allocatable, dimension(:) :: displs, recvcnt
      real(8), allocatable, dimension(:) :: tmp1
      complex(8), allocatable, dimension(:) :: tmp2

      n = ied - ist + 1
      allocate(tmp1(n))
      allocate(displs(0:nprocs-1), recvcnt(0:nprocs-1))


      rem = mod(ip%Egrid+1, nprocs)

      do m=0, nprocs-1
        if (rem == 0 .or. m < rem) then
          recvcnt(m) = ip%Egrid / nprocs + 1
        else
          recvcnt(m) = ip%Egrid / nprocs
        end if
        displs(m) = sum(recvcnt(:m-1))
      end do


      tmp1 = sig_fus(ist:ied)
      call mpi_gatherv(tmp1,n,mpi_real8,sig_fus,recvcnt,displs, &
                       mpi_real8,0, mpi_comm_world,ierr)
      tmp1 = spin(ist:ied)
      call mpi_gatherv(tmp1,n,mpi_real8,spin,recvcnt,displs, &
                       mpi_real8,0,mpi_comm_world,ierr)
      tmp1 = sig_iel(ist:ied)
      call mpi_gatherv(tmp1,n,mpi_real8,sig_iel,recvcnt,displs, &
                       mpi_real8,0,mpi_comm_world,ierr)

      do m=0, nprocs-1
        if (rem == 0 .or. m < rem) then
          recvcnt(m) = (ip%Egrid/nprocs + 1) * ip%Nch
        else
          recvcnt(m) = ip%Egrid/nprocs * ip%Nch
        end if
        displs(m) = sum(recvcnt(:m-1))
      end do

      deallocate(tmp1)
      allocate(tmp1(n*ip%Nch))

      tmp1 = pack(sig_iel_n(:,ist:ied),.true.)
      call mpi_gatherv(tmp1,recvcnt(myrank),mpi_real8,       &
                       sig_iel_n,recvcnt,displs,mpi_real8,0, &
                       mpi_comm_world,ierr)

      do m=0, nprocs-1
        if (rem == 0 .or. m < rem) then
          recvcnt(m) = (ip%Egrid/nprocs + 1) * ip%nth
        else
          recvcnt(m) = ip%Egrid/nprocs * ip%nth
        end if
        displs(m) = sum(recvcnt(:m-1))
      end do

      allocate(tmp2(ip%nth*n))

      tmp2 = pack(fc(:,ist:ied),.true.)
      call mpi_gatherv(tmp2,recvcnt(myrank),mpi_double_complex, &
                       fc,recvcnt,displs,mpi_double_complex,    &
                       0,mpi_comm_world,ierr)


      do m=0, nprocs-1
        if (rem == 0 .or. m < rem) then
          recvcnt(m) = (ip%Egrid/nprocs + 1) * ip%nth * ip%Nch
        else
          recvcnt(m) = ip%Egrid/nprocs * ip%nth * ip%Nch
        end if
        displs(m) = sum(recvcnt(:m-1))
      end do

      deallocate(tmp2)
      allocate(tmp2(ip%nth*n*ip%Nch))

      tmp2 = pack(fN(:,:,ist:ied),.true.)
      call mpi_gatherv(tmp2,recvcnt(myrank),mpi_double_complex, &
                       fN,recvcnt,displs,mpi_double_complex,    &
                       0,mpi_comm_world,ierr)

      deallocate(tmp1,tmp2)
    end subroutine
!**********************************************************************!
    subroutine record_partial_results(id)
      implicit none
      integer, intent(in) :: id
      integer :: i, n, m
      real(8) :: E, tmp1
      character(len=50), parameter :: FM10='(1x,f10.3,3es11.3)'

      Esig = Ea * sig_fus
      dsig_R = abs(fc) ** 2 * 10.0d0
      dsig_qel_n(1,:,:) = abs(fc + fN(1,:,:)) ** 2 * 10.0d0
      do n=2, ip%Nch
        dsig_qel_n(n,:,:) = abs(fN(n,:,:)) ** 2 * 10.0d0
      end do
      dsig_iel = sum(dsig_qel_n(2:ip%Nch,:,:),dim=1) 
      dsig_qel = dsig_iel + dsig_qel_n(1,:,:)

      open(7,file=trim(dir)//'/tmp_fusion.dat')
      open(8,file=trim(dir)//'/tmp_qel.dat')
      do i=ist, id
        write(7,*) Ea(i), sig_fus(i)
        write(8,FM10) Ea(i), dsig_qel(ip%nth,i)/dsig_R(ip%nth,i)
!       write(8,FM10) Ea(i), dsig_qel(ip%nth-2,i)/dsig_R(ip%nth-2,i), &
!                         dsig_qel(ip%nth-1,i)/dsig_R(ip%nth-1,i), &
!                         dsig_qel(ip%nth,i)/dsig_R(ip%nth,i)
      end do
      close(7)
      close(8)

    end subroutine
!**********************************************************************!
    subroutine output()
      implicit none
      integer :: t, di_u, di_l
      real(8) :: E, Dfus, Dqel
      character(len=45), parameter :: FM1='(1x,a,f8.3,a)'
      character(len=45), parameter :: FM2='(1x,a,es13.4,a)'
      character(len=45), parameter :: FM3='(1x,a,i3,a,es13.4,a)'
      character(len=50), parameter :: FM4='(1x,f8.3,3es13.4)'
      character(len=50), parameter :: FM5='(1x,f7.3,f8.2,3x,a,es14.4)'
      character(len=50), parameter :: FM6='(19x,i3,1x,2es14.4)'
      character(len=50), parameter :: FM7='(8x,f8.2,3x,a,es14.4)'
      character(len=50), parameter :: FM8='(1x,2f10.3,es14.4)'
      character(len=50), parameter :: FM9='(1x,2f10.3,es11.3)'
      character(len=50), parameter :: FM10='(1x,f10.3,es11.3)'
      character(len=500) :: c, c1, c2, c3, c4

!-- Fusion cross section -------------------------------
!     open(7,file=trim(dir)//'/fusion.dat')
      open(7,file='fusion.dat')
      do i=1, ip%Egrid+1
        write(7,*) Ea(i), sig_fus(i)
      end do
      close(7)

!-- Fusion barrier distribution -----------------------------
!     open(7,file=trim(dir)//'/fus_bar_dist.dat')
      open(7,file='fus_bar_dist.dat')
      do i=ip%di+1, ip%Egrid-ip%di+1
        Dfus = 2.0d0 * ((Esig(i+ip%di)-Esig(i)) / (Ea(i+ip%di)-Ea(i))  &
                     -  (Esig(i)-Esig(i-ip%di)) / (Ea(i)-Ea(i-ip%di))) &
                     / (Ea(i+ip%di) - Ea(i-ip%di))
        write(7,FM4) Ea(i), Dfus
      end do
      close(7)

!-- Energy dependence of differential scattering cross sections -----
      do n=1, ip%nth
        t = nint(angl(n) * 1800.0d0 / PI)   ! per 10 degree
        if (mod(t,100) == 0) then
          write(c,'(f6.2)') angl(n) * 180.0d0 / PI
          c1 = trim(dir2)//'/qel_'//trim(adjustl(c))//'deg.dat'
!         c2 = trim(dir2)//'/el_'//trim(adjustl(c))//'deg.dat'
!         c3 = trim(dir2)//'/iel_'//trim(adjustl(c))//'deg.dat'
          open(7,file=trim(c1))
!         open(8,file=trim(c2))
!         open(9,file=trim(c3))
          do i=1, ip%Egrid+1
            E = 2.0d0 * Ea(i) * sin(0.5d0 * angl(n)) &
                     / (1.0d0 + sin(0.5d0 * angl(n)))
            write(7,FM9) Ea(i), E, dsig_qel(n,i)/dsig_R(n,i)
!           write(8,FM10) Ea(i), dsig_qel_n(1,n,i)/dsig_R(n,i)
!           write(9,FM10) Ea(i), dsig_iel(n,i)/dsig_R(n,i)
          end do
          close(7)
!         close(8)
!         close(9)
        end if
      end do

!-- qel. barrier distribution -----------------------------
      if (mod(ip%di,2) == 0) then
        di_u = ip%di / 2
        di_l = di_u
      else
        di_u = ceiling(dble(ip%di)/2.0d0)
        di_l = ip%di / 2
      end if
      do n=1, ip%nth
        t = nint(angl(n) * 1800.0d0 / PI)   ! per 10 degree
        if (t > 1100 .and. mod(t,100) == 0) then
          write(c,'(f6.2)') angl(n) * 180.0d0 / PI
         c1 = trim(dir2)//'/qel_bar'//trim(adjustl(c))//'deg.dat'
          open(7,file=trim(c1))
          do i=di_l+1, ip%Egrid-di_u+1
            E = 0.5d0 * (Ea(i+di_u) + Ea(i-di_l))
            E = 2.0d0 * E * sin(0.5d0 * angl(n))            &
                     / (1.0d0 + sin(0.5d0 * angl(n)))
            Dqel = - (dsig_qel(n,i+di_u) / dsig_R(n,i+di_u)   &
                   -  dsig_qel(n,i-di_l) / dsig_R(n,i-di_l))  &
                   / (Ea(i+di_u) - Ea(i-di_l))
            write(7,FM8) Ea(i), E, Dqel
          end do
          close(7)
        end if
      end do


!-- angular distribution of elastic differential scattering cross section
!     do i=1, ip%Egrid+1
!       if (mod(nint(Ea(i)*10.0d0),50) == 0) then     ! per 5 MeV
!         write(c,'(f6.2)') Ea(i)
!         open(7,file=trim(dir3)//'/diff_el_'//trim(adjustl(c))//'MeV.dat')
!         do n=1, ip%nth
!           write(7,*) angl(n)*180.0d0/PI, dsig_qel_n(1,n,i)/dsig_R(n,i)
!         end do
!         close(7)
!       end if
!     end do

!   return
!-- Q-value distribution for scattering
      do i=1, ip%Egrid+1
!       if (mod(nint(Ea(i)*10.0d0),50) == 0) then   ! per 5 MeV
          do n=1, ip%nth
            t = nint(angl(n) * 1800.0d0 / PI)
            if (t > 1100 .and. mod(t,50) == 0) then    ! per 5 degree
!           if (t > 1100 .and. mod(t,100) == 0) then   ! per 10 degree
              write(c1,'(f6.2)') Ea(i)
              write(c2,'(f6.2)') angl(n) * 180.0d0 / PI
              c = trim(dir4)//'/Qdist_'//trim(adjustl(c1))//'MeV_'// &
                  trim(adjustl(c2))//'deg.dat'
              open(7,file=trim(c))
                do m=1, ip%Nch
                  write(7,*) cm%e_n(m), dsig_qel_n(m,n,i), &
                             dsig_qel_n(m,n,i)/dsig_R(n,i)
                end do
              close(7)
            end if
          end do
!       end if
      end do


    end subroutine
!**********************************************************************!
    subroutine qel_corrected_angle()
      implicit none
      integer :: i, n, m
      integer :: s, di_u, di_l
      integer :: t(ip%Nch)
      integer, dimension(ip%Nch) :: k
      integer, dimension(ip%nth,ip%Egrid+1) :: kcm2
      real(8), allocatable, dimension(:,:,:) :: dsig_qel_n_lab
      real(8), allocatable, dimension(:,:) :: dsig_R_lab
      real(8), dimension(ip%Nch) :: Kcm, vcm, tcm, tlab
      real(8) :: E, Dqel, Elab, V, E2
      character(len=50), parameter :: FM8='(1x,2f10.3,2es14.4)'
      character(len=50), parameter :: FM9='(1x,2f10.3,2es11.3)'
      character(len=500) :: c, c1, c2, c3, c4

      allocate(dsig_qel_n_lab(ip%Nch,ip%nth,ip%Egrid+1))
      allocate(dsig_R_lab(ip%nth,ip%Egrid+1))
      dsig_qel_n_lab = 0.0d0
!     write(6,*) sum(dsig_qel(1,:))*ip%dE
      do i=1, ip%Egrid+1
        Kcm = Ea(i) - cm%e_n
        Elab = (ip%Ap + ip%At) / ip%At * Ea(i)
        vcm = sqrt(2.0d0*ip%rmass*Kcm) / (ip%Ap*mass)
        V = sqrt(2.0d0*(ip%Ap*mass)*Elab) / ((ip%Ap+ip%At)*mass)
        do n=1, ip%nth
!         tlab = vcm*sin((PI-angl(n))) / (vcm*cos((PI-angl(n))) - V)
!         tlab = (PI-atan(tlab)) * 180.0d0 / PI
!         tlab = nint(tlab * 10.0d0) / 10.0d0
!         k = nint((tlab - ip%thmin) / ip%dth) + 1

          tcm = -V/vcm*sin(angl(n))**2 + cos(angl(n))*sqrt(1.0d0-(V/vcm*sin(angl(n)))**2)
          tcm = acos(tcm) * 180.0d0 / PI
          tcm = nint(tcm * 10.0d0) / 10.0d0
          k = nint((tcm - ip%thmin) / ip%dth) + 1
          kcm2(n,i) = k(1)
!         write(6,'(1x,3f8.2)') angl(n)*180.0d0/PI, angl(k(1))*180.0d0/PI,(angl(n)-angl(k(1)))*180.0d0/PI
          do m=1, ip%Nch
            if (k(m) > 0 .and. k(m) <= ip%nth) then
              dsig_qel_n_lab(m,n,i) = dsig_qel_n(m,k(m),i)
!             dsig_qel_n_lab(m,k(m),i) = dsig_qel_n(m,n,i)
            end if
          end do
          dsig_R_lab(n,i) = dsig_R(k(1),i)
        end do
!       stop
      end do
      dsig_qel = sum(dsig_qel_n_lab(:,:,:),dim=1) 
!     write(6,*) sum(dsig_qel(1,:))*ip%dE

!     write(6,*) 'angle arranged.'

!-- Energy dependence of differential scattering cross sections -----
      do n=1, ip%nth
        s = nint(angl(n) * 1800.0d0 / PI)   ! per 10 degree
        if (mod(s,100) == 0) then
          write(c,'(f6.2)') angl(n) * 180.0d0 / PI
          c1 = trim(dir2)//'/qel_'//trim(adjustl(c))//'deg_lab.dat'
          open(7,file=trim(c1))
          do i=1, ip%Egrid+1
            E = 2.0d0 * Ea(i) * sin(0.5d0 * angl(n)) &
                     / (1.0d0 + sin(0.5d0 * angl(n)))
            E2= 2.0d0 * Ea(i) * sin(0.5d0 * angl(kcm2(n,i))) &
                     / (1.0d0 + sin(0.5d0 * angl(kcm2(n,i))))
            write(7,FM9) Ea(i), E, dsig_qel(n,i)/dsig_R_lab(n,i), E2
          end do
          close(7)
        end if
      end do

!-- qel. barrier distribution -----------------------------
      if (mod(ip%di,2) == 0) then
        di_u = ip%di / 2
        di_l = di_u
      else
        di_u = ceiling(dble(ip%di)/2.0d0)
        di_l = ip%di / 2
      end if
      do n=1, ip%nth
        s = nint(angl(n) * 1800.0d0 / PI)   ! per 10 degree
        if (s > 1100 .and. mod(s,100) == 0) then
          write(c,'(f6.2)') angl(n) * 180.0d0 / PI
         c1 = trim(dir2)//'/qel_bar'//trim(adjustl(c))//'deg_lab.dat'
          open(7,file=trim(c1))
          do i=di_l+1, ip%Egrid-di_u+1
            E = 0.5d0 * (Ea(i+di_u) + Ea(i-di_l))
            E2= 2.0d0 * E * sin(0.5d0 * angl(kcm2(n,i))) &
                 / (1.0d0 + sin(0.5d0 * angl(kcm2(n,i))))
            E = 2.0d0 * E * sin(0.5d0 * angl(n))            &
                     / (1.0d0 + sin(0.5d0 * angl(n)))
            Dqel = - (dsig_qel(n,i+di_u) / dsig_R_lab(n,i+di_u)   &
                   -  dsig_qel(n,i-di_l) / dsig_R_lab(n,i-di_l))  &
                   / (Ea(i+di_u) - Ea(i-di_l))
            write(7,FM8) Ea(i), E, Dqel, E2
          end do
          close(7)
        end if
      end do

      deallocate(dsig_qel_n_lab)
      deallocate(dsig_R_lab)

    end subroutine
!**********************************************************************!
  end program

