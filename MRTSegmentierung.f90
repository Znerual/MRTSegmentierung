    !  MRTSegmentierung.f90
    !
    !  FUNCTIONS:
    !  MRTSegmentierung - Entry point of console application.
    !

    !****************************************************************************
    !
    !  PROGRAM: MRTSegmentierung
    !
    !  PURPOSE:  Entry point for the console application.
    !
    !****************************************************************************

    program MRTSegmentierung
    use IFPORT
    implicit none

    !dimension of the image
    integer Lx,Ly
    parameter (Lx=254,Ly=333)
    
    integer, parameter :: NN=4

    integer Inn(NN)                    !/* Nearest neighbor array I */
    parameter (Inn=(/1,-1,0,0/))
    integer Jnn(NN)                       !/* Nearest neighbor array J */
    parameter (Jnn= (/0,0,1,-1/))
    real*8 :: JJ = 1d0

    integer :: i
    real*8 :: acceptRate
    real*8 :: startTime, endTime
    real*8 :: midRR
        
    !Sweep amount
    integer, parameter :: warm_up = 500, step=50
    !mean values and standard deviation
    integer MW(5)
    !               BG  WM GM  CSF  SB
    parameter (MW=(/30,426,602,1223,167/))
    integer SIGMA(5)
    parameter (SIGMA=(/30,59,102,307,69/))

    
    !temperatures
    real*8 :: t_i = 4d0, T, t_f= 1d-1, beta
    real*8 :: lambda = 1.5d0
    integer :: confMax(Lx,Ly), confMin(Lx,Ly)
    
    !MR Picture
    integer :: mag(Lx,Ly), conf(Lx,Ly)
    integer :: MR_correct(Lx, Ly)
    real*8 :: dif(5)
    
    !Open files
    open(unit=33, file="SimMRimage.dat", action="read")
    open(unit=34, file="TAcceptRate.dat", action="write")
    open(unit=36, file="error.dat", action="write")
    open(unit=37, file="CorrectSegImage.dat", action="read")
    
    call RANDOM_SEED()
    
    call CPU_TIME(startTime)
    
    call input(33, mag, Lx,Ly)
    close(33)
    
    call input(37, MR_correct, Lx, Ly)
    close(37)
    call CPU_TIME(endTime)
    print*,"Input done in ", endTime - startTime, " seconds"
    !Initialise Configuration
    
    call CPU_TIME(startTime)
    call random_init(conf, Lx,Ly)
    call CPU_TIME(endTime)
    print*,"Initialising done in ", endTime - startTime, " seconds"
    
    call CPU_TIME(startTime)
    
    beta=1/t_i
    do i = 1, warm_up
        call do_sweep(mag, Inn,Jnn, JJ,NN, Lx, Ly, beta,  acceptRate, midRR, MW, Sigma, conf)       
    end do
    call CPU_TIME(endTime)
    print*,"Warm up done in ", endTime - startTime, " seconds with <r>: ", midRR
    call CPU_TIME(startTime)
    T=t_i
    do while(T > t_f)
        
        beta = 1d0/T
        do i = 1, step
            call do_sweep(mag, Inn,Jnn, JJ,NN, Lx, Ly, beta,  acceptRate, midRR, MW, Sigma, conf)
            
        end do
        write(34, *) T, acceptRate
        if (acceptRate < 1d-4) then
            exit
        endif
        !print*, "T: ", T, "rate: ", acceptRate
        T = T/lambda
    end do
    call CPU_TIME(endTime)
    print*,"Calculation done in ", endTime - startTime, " seconds"
    open(unit=35, file="SegSA.dat", action="write", status="replace")
    call output(35, conf, Lx,Ly)
    close(35)
    
    write(36, *) "J = ", JJ
    call calcError(Lx,Ly, MR_correct, conf, dif)
    write(36,*) dif
    
    ! J=0
    call CPU_TIME(startTime)
    call random_init(conf, Lx,Ly)
    call CPU_TIME(endTime)
    print*,"Initialising done in ", endTime - startTime, " seconds"
    
    call CPU_TIME(startTime)
    
    beta=1/t_i
    do i = 1, warm_up
        call do_sweep(mag, Inn,Jnn, 0d0,NN, Lx, Ly, beta,  acceptRate, midRR, MW, Sigma, conf)       
    end do
    call CPU_TIME(endTime)
    print*,"Warm up done in ", endTime - startTime, " seconds with <r>: ", midRR
    call CPU_TIME(startTime)
    T=t_i
    do while(T > t_f)
        
        beta = 1d0/T
        do i = 1, step
            call do_sweep(mag, Inn,Jnn, 0d0,NN, Lx, Ly, beta,  acceptRate, midRR, MW, Sigma, conf)
            
        end do
        write(34, *) T, acceptRate
        if (acceptRate < 1d-4) then
            exit
        endif
        print*, "T: ", T, "rate: ", acceptRate
        T = T/lambda
    end do
    call CPU_TIME(endTime)
    print*,"Calculation done in ", endTime - startTime, " seconds"
    open(unit=39, file="SegLocal.dat", action="write", status="replace")
    call output(39, conf, Lx,Ly)
    close(39)
    
    write(36,*) ''
    write(36, *) "J = 0"
    call calcError(Lx,Ly, MR_correct, conf, dif)
    write(36,*) dif
    
    !J = 10
    JJ = 0.9
    call CPU_TIME(startTime)
    call random_init(conf, Lx,Ly)
    call CPU_TIME(endTime)
    print*,"Initialising done in ", endTime - startTime, " seconds"
    
    call CPU_TIME(startTime)
    
    beta=1/t_i
    do i = 1, warm_up
        call do_sweep(mag, Inn,Jnn, 0d0,NN, Lx, Ly, beta,  acceptRate, midRR, MW, Sigma, conf)       
    end do
    call CPU_TIME(endTime)
    print*,"Warm up done in ", endTime - startTime, " seconds with <r>: ", midRR
    call CPU_TIME(startTime)
    T=t_i
    do while(T > t_f)
        
        beta = 1d0/T
        do i = 1, step
            call do_sweep(mag, Inn,Jnn, 0d0,NN, Lx, Ly, beta,  acceptRate, midRR, MW, Sigma, conf)
            
        end do
        write(34, *) T, acceptRate
        if (acceptRate < 1d-4) then
            exit
        endif
        print*, "T: ", T, "rate: ", acceptRate
        T = T/lambda
    end do
    call CPU_TIME(endTime)
    print*,"Calculation done in ", endTime - startTime, " seconds"
    
    write(36,*) ''
    write(36, *) "J = ", JJ
    call calcError(Lx,Ly, MR_correct, conf, dif)
    write(36,*) dif
    
    close(34)
    close(36)
    
    
    
    
    
    
    
    
    
    
    
    contains
    
    subroutine random_init(conf, Lx,Ly)
    
    integer, intent(in) :: Lx, Ly
    integer, intent(inout) :: conf(Lx,Ly)

    integer :: i,j
    real*8 :: zufall
    do i=1,Lx
        do j=1,Ly
            call RANDOM_NUMBER(zufall)
            conf(i,j) = INT(zufall*5.d0)+1
        enddo
    enddo

    return
    end subroutine random_init

    subroutine do_sweep(mag, Inn,Jnn, JJ,NN, Lx, Ly, beta,  acceptRate, midRR, MW, Sigma, conf)
        integer, intent(in) :: Lx, Ly
        integer, dimension(Lx, Ly), intent(inout) :: mag, conf
        integer, dimension(NN), intent(in) :: Inn, Jnn
        real*8, intent(in) :: JJ
        integer, intent(in) :: NN
        real*8, intent(in) :: beta
        real*8, intent(out) :: acceptRate
        integer :: MW(5), Sigma(5)
        integer :: oldconf 
        integer :: i, j, k
        real*8 :: deltaE
        real*8 :: midRR
        real*8 :: zufall, RR
        
        
        acceptRate = 0_8
        midRR = 0_8
        
        do k = 1, Lx * Ly
            call RANDOM_NUMBER(zufall) !zufall = rand() !immer gleiche Zufallzahlen
            i = ceiling(size(mag,1) * zufall) 
            call RANDOM_NUMBER(zufall) 
            j = ceiling(size(mag,2) * zufall) !Gitterplatz bestimmt
            
            
            deltaE = CalcSubE(mag, conf,i,j,Lx,Ly,Jnn,Inn,NN,JJ,MW, Sigma)
            
            !random tissue in 1..5:
            call RANDOM_NUMBER(zufall)
            oldconf = conf(i, j)
            conf(i,j)=INT(zufall*5.d0)+1 !configuration flip
             
            deltaE =  CalcSubE(mag, conf,i,j,Lx,Ly,Jnn,Inn,NN,JJ,MW, Sigma) - deltaE
            
            
            
            !check for floating point overflow
            if ((deltaE) * beta .gt. 2d1)  then
                RR=0.d0
            else
                if ((deltaE)*beta .lt. 5d-2)  then
                    RR=1.d0
                else
                    RR=exp(-(deltaE)*beta)
                end if
            end if
            
            RR = min(1d0, RR)
            midRR = midRR + RR
            
            call RANDOM_NUMBER(zufall)
            if (zufall < RR) then !Falls rand < min(1,r) -> neue Konfig, dh wir behalten unsere ände    
                acceptRate = acceptRate + 1
            else
                conf(i,j) = oldconf !rand > ... -> zurück auf alte Konfig
            end if
        end do
        
        acceptRate = acceptRate/dble(Lx*Ly)
        midRR =  midRR / dble(Lx*Ly)
    end subroutine do_sweep
    
    subroutine output(unit, conf, Lx, Ly)
        integer, intent(in) :: Lx,Ly
        integer, intent(in) :: conf(Lx,Ly)
        integer :: unit
        
        integer :: i, j
        !write out image
        do i=1,Lx
            do j=1,Ly
                write(unit,'(3I5)') i,j,conf(i,j)
            enddo
            write (unit,*)
        enddo
    end subroutine output
    
    subroutine input(unit, mag, Lx, Ly)
        integer, intent(in) :: Lx,Ly
        integer, intent(inout) :: mag(Lx,Ly)
        integer :: unit
        
        integer :: i, j
        integer :: dummyi, dummyj
        !read in image
        do i=1,Lx
            do j=1,Ly
                read(unit,'(3I5)') dummyi,dummyj,mag(i,j)
            enddo
            read (unit,*)
        enddo
    
    end subroutine input
    
    real function CalcE(mag, conf,Lx,Ly,Jnn,Inn,NN,JJ,mean, deviation)
        integer mag(Lx,Ly), conf(Lx, Ly), Jnn(NN),Inn(NN),NN
        integer ::  mean(5), deviation(5)
        real*8 :: JJ
        integer :: i,j,k,Inew,Jnew, Lx,Ly
        real*8 :: Energy
        

        !/* Determine the initial energy */
        Energy = 0.0
        do i=1, Lx
            do j=1, Ly

                
                if (JJ .ne. 0) then
                    !/* Loop over nearest neighbors */
                    do k=1, NN
                        Inew = i + Inn(k)
                        Jnew = j + Jnn(k)

                            !new boundary condiition (period possible since BG->BG)
                            !/* Check periodic boundary conditions */
                            if (Inew .le. 0) then
                                Inew = Lx
                            else
                                if(Inew .gt. Lx ) then
                                    Inew = 1
                                endif
                            endif
                            if (Jnew .le. 0) then
                                Jnew = Ly
                            else
                                if(Jnew .gt. Ly) then
                                    Jnew = 1
                                endif
                            endif
                        !/* Update the energy */
                        if (conf(i,j) == conf(Inew, Jnew)) then
                            Energy = Energy-JJ
                        end if
                    enddo
                end if
                !/*Calculate the contribution from the field H */
                Energy = Energy + (mag(i,j) - mean(conf(i,j)))**2_8/(2_8*deviation(conf(i,j))**2) + &
                dlog(dble(deviation(conf(i,j)))) ;
            enddo
        enddo

        !/* Account for double counting */
        Energy = Energy/2.d0;

        CalcE=Energy
        return
    end function CalcE
    
    subroutine calcError(Lx,Ly, correct_conf, conf, dif)
        integer, intent(in) :: Lx, Ly
        integer, intent(in) :: correct_conf(Lx,Ly), conf(Lx,Ly)
        real*8, intent(out) :: dif(5)
        integer :: amount(5)
        
        integer :: i, j
        do i = 1, Lx
            do j= 1, Ly
                amount(correct_conf(i,j)) = amount(correct_conf(i,j)) + 1
                if (correct_conf(i,j) .ne. conf(i,j)) then
                    dif(conf(i,j)) = dif(conf(i,j)) + 1
                end if
            end do
        end do
        
        dif = dif / amount
    end subroutine calcError
    
    real function CalcSubE(mag, conf,i, j,Lx, Ly,Jnn,Inn,NN,JJ,mean, deviation)
        integer mag(Lx,Ly), conf(Lx, Ly), Jnn(NN),Inn(NN),NN
        integer ::  mean(5), deviation(5)
        real*8 :: JJ
        integer :: i,j
        integer :: k,Inew,Jnew, Lx,Ly
        real*8 :: Energy
        

        !/* Determine the initial energy */
        Energy = 0.0
                
        if (JJ .ne. 0) then
            !/* Loop over nearest neighbors */
            do k=1, NN
                Inew = i + Inn(k)
                Jnew = j + Jnn(k)
                !new boundary condiition (period possible since BG->BG)
                !/* Check periodic boundary conditions */
                if (Inew .le. 0) then
                    Inew = Lx
                else
                    if(Inew .gt. Lx ) then
                        Inew = 1
                    endif
                endif
                if (Jnew .le. 0) then
                    Jnew = Ly
                else
                    if(Jnew .gt. Ly) then
                        Jnew = 1
                    endif
                endif
                !/* Update the coupled energy */
                if (conf(i,j) == conf(Inew, Jnew)) then
                    Energy = Energy-JJ
                end if
            enddo
        end if
        !/*Calculate the contribution from the intensities */
        Energy = Energy + (mag(i,j) - mean(conf(i,j)))**2_8/(2_8*deviation(conf(i,j))**2) + &
        dlog(dble(deviation(conf(i,j)))) ;
    

        CalcSubE=Energy
        return
    end function CalcSubE
    
    end program MRTSegmentierung

