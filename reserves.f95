!subroutine energygradient(na,dux,duy,duz)
!    use lattice_mod
!    use interaction_mod
!    use carbon_IIA_mod
!    integer i, k, na
!
!    real(8) delX,delY,delZ,dux,duy,duz,r
!    real(8) dphi_i, phi_i, x_dphi_i, y_dphi_i, z_dphi_i, reciprocal_phi !for summs
!    real(8) cos_alp, cos_bet, cos_gam
!    real(8) pairwise_deriv1,edensity,edensity_deriv1,dpoti,dpotuC
!    phi_i    = 0.0d0
!    x_dphi_i = 0.0d0
!    y_dphi_i = 0.0d0
!    z_dphi_i = 0.0d0
!    dux = 0d0
!    duy = 0d0
!    duz = 0d0
!!R_curr(1,na) = R_curr(1,na) -0.00065_8 !test for forces
!!R_curr(2,na) = R_curr(2,na) -0.00025_8 !
!!R_curr(3,na) = R_curr(3,na) -0.00035_8 !
!!write(*,*) inact(0,na)
!!F_x = -0.5*SUM_of_k`[V_r(r_k)*x_k/r_k] + 0.5*((SUM_of_k`[PHI(r_k)])**-1)*SUM_of_k`[Phi_r(r_k)*x_k/r_k]
!    do k=1,inact(0,na)
!        i=inact(k,na)
!        delX    = R_curr(1,na)-R_curr(1,i)!+1d-30
!        delY    = R_curr(2,na)-R_curr(2,i)!+1d-30
!        delZ    = R_curr(3,na)-R_curr(3,i)!+1d-30
!        r       = dsqrt(delX**2+delY**2+delZ**2 +1d-30)
!        cos_alp = delX/r
!        cos_bet = delY/r
!        cos_gam = delZ/r
!
!        dpoti  = pairwise_deriv1(r) !dpotu(r)
!        dphi_i = edensity_deriv1(r)
!        !phi_i  = edensity(r) + phi_i
!
!        x_dphi_i = x_dphi_i + dphi_i*cos_alp
!        y_dphi_i = y_dphi_i + dphi_i*cos_bet
!        z_dphi_i = z_dphi_i + dphi_i*cos_gam
!        dux = dux - dpoti*cos_alp
!        duy = duy - dpoti*cos_bet
!        duz = duz - dpoti*cos_gam
!    enddo
!    reciprocal_phi = 0d0*dsqrt(1d0/(phi_i+1d-30))  !that is not fair embedding function
!
!    dux = 5d-1*(dux + x_dphi_i*reciprocal_phi ) !first 0.5 is from formula
!    duy = 5d-1*(duy + y_dphi_i*reciprocal_phi ) !see the article
!    duz = 5d-1*(duz + z_dphi_i*reciprocal_phi ) !second is from derivative of sqrt
!        delX    = R_curr(1,na)-X_carb(1)!carbon
!        delY    = R_curr(2,na)-Y_carb(1)
!        delZ    = R_curr(3,na)-Z_carb(1)
!        !Write(*,*) Z_carb(1), "no", delz
!        r       = dsqrt(delX**2+delY**2+delZ**2 +1d-30)
!        dpoti   = dpotuC(r)
!        r       = 1d0/r
!        dux = dux - dpoti*delX*r !r is 1/r
!        duy = duy - dpoti*delY*r
!        duz = duz - dpoti*delZ*r
!endsubroutine energygradient
!
!!subroutine      relax_atom(na)
!!    use lattice_mod
!!    use phys_parameters_mod
!!    integer na, decrease, moves, i
!!    integer, parameter :: max_decrease = 14, max_moves = 60, retries=1!00
!!    real(8) ::  step = 1d-1, denominator = 125d-3, small_move = 1d-07!d-7 is ok but a bit rude
!!    real(8)     fx,fy,fz,xsh,ysh,zsh,normal
!!    real(8)     x_prev,y_prev,z_prev
!!    real(8)     energy_current, energy_best, eatom, rndm(retries+max_moves)
!!!    OPEN(11, FILE = 'GRADIATION.TXT')
!!    fx  =0d0; fy  =0d0; fz  =0d0
!!    xsh =0d0; ysh =0d0; zsh =0d0
!!!    Write(*,*) R_curr(1,na),R_curr(2,na),R_curr(3,na)
!!!    Energy_current = eatom(na)
!!!    Write(*,*) "performing a relaxation of atom# ", na, energy_current
!!!R_curr(1,na)=R_curr(1,na)-193227d-6
!!!R_curr(2,na)=R_curr(2,na)-325123d-6
!!!R_curr(3,na)=R_curr(3,na)-264123d-6
!!!    Write(*,*) R_curr(1,na),R_curr(2,na),R_curr(3,na)
!!!    energy_current = eatom(na)
!!!    write(*,*) "performing a relaxation of atom# ", na, energy_current
!!    call random_number(rndm)
!!!    Write(*,*) R_curr(1,na),R_curr(2,na),R_curr(3,na)
!!    Do i=1,retries
!!        step = 2d-1*rndm(i)
!!        do decrease = 1,max_decrease
!!            step = step * denominator
!!            do moves = 1,max_moves
!!                energy_current = eatom(na)
!!                call energygradient(na,fx,fy,fz)
!!                normal = dsqrt(fx*fx+fy*fy+fz*fz)
!!                    if (normal .lt. 1d-14) then
!!                        !write(*,*) fx,fy,fz,"why so small normalizator. can not deal with this"
!!                        cycle
!!                    endif
!!                normal=1d0/normal
!!                fx=fx*normal
!!                fy=fy*normal
!!                fz=fz*normal
!!                xsh=fx*step
!!                ysh=fy*step
!!                zsh=fz*step
!!                R_curr(1,na)=R_curr(1,na)+xsh
!!                R_curr(2,na)=R_curr(2,na)+ysh
!!                R_curr(3,na)=R_curr(3,na)+zsh
!!                if (eatom(na) .gt. energy_current) then
!!                    R_curr(1,na)=R_curr(1,na)-xsh*rndm(moves)!*6d-1
!!                    R_curr(2,na)=R_curr(2,na)-ysh*rndm(moves)!*6d-1
!!                    R_curr(3,na)=R_curr(3,na)-zsh*rndm(moves)!*6d-1
!!!                        Write(*,*) "wow",xsh,ysh,zsh
!!!                        ELSE
!!!                        Write(*,*) "wwwwqeqrsw",xsh,ysh,zsh
!!                endif
!!                !WRITE(11,*) R_curr(1,na),R_curr(2,na),R_curr(3,na)
!!            enddo
!!            if((dabs(xsh)+dabs(ysh)+dabs(zsh)) .lt. small_move) exit!e-7 is ok
!!        enddo
!!    Enddo
!!!    Write(*,*) R_curr(1,na),R_curr(2,na),R_curr(3,na)
!!!    Pause
!!    !CLOSE(11)
!!!    Energy_current = eatom(na)
!!!    Write(*,*) "performed  a relaxation of atom# ", na, energy_current
!!!    Write(*,*) R_curr(1,na),R_curr(2,na),R_curr(3,na)
!!endsubroutine   relax_atom

!subroutine wo_matrixes
!!
!!    use latticemod
!!    use slicemod
!!
!!
!!    integer i
!!    character (LEN=28) matrixfilenameUX1, matrixfilenameUX0
!!    character (LEN=28) matrixfilenameUY1, matrixfilenameUY0
!!
!!    matrixWOcounter = matrixWOcounter + 1
!!write(*,*)          "=== writing the matrixes out"
!!write(*,'(A,I2.2)')," === we already did it times: ", matrixWOcounter-1
!!
!!    write (matrixfilenameUX1, 14) ,"Matrix-stage" ,matrixWOcounter,("-z1-Ux-" &
!!    // dislocationtype // ".txt")
!!    open (21, file = matrixfilenameUX1)
!!
!!    write (matrixfilenameUY1, 14) ,"Matrix-stage" ,matrixWOcounter,("-z1-Uy-" &
!!    // dislocationtype // ".txt")
!!    open (22, file = matrixfilenameUY1)
!!
!!    write (matrixfilenameUX0, 14) ,"Matrix-stage" ,matrixWOcounter,("-z0-Ux-" &
!!    // dislocationtype // ".txt")
!!    open (23, file = matrixfilenameUX0)
!!
!!    write (matrixfilenameUY0, 14) ,"Matrix-stage" ,matrixWOcounter,("-z0-Uy-" &
!!    // dislocationtype // ".txt")
!!    open (24, file = matrixfilenameUY0)
!!
!!if (dislocationtype .eq. '001') then
!!    do i=1,sq
!!                !write out matrix of evasions for z=0.5*a0 plane
!!        write(21, 10), X( m05((1+(i-1)*sq ):(i*sq )) ) - Xi( m05((1+(i-1)*sq ):(i*sq )) )
!!
!!        write(22, 10), Y( m05((1+(i-1)*sq ):(i*sq )) ) - Yi( m05((1+(i-1)*sq ):(i*sq )) )
!!    enddo
!!    do i=1,sqq
!!                !write out matrix of evasions for z=0     plane
!!        write(23, 10), X( m00((1+(i-1)*sqq):(i*sqq)) ) - Xi( m00((1+(i-1)*sqq):(i*sqq)) )
!!
!!        write(24, 10), Y( m00((1+(i-1)*sqq):(i*sqq)) ) - Yi( m00((1+(i-1)*sqq):(i*sqq)) )
!!    enddo
!!endif
!!
!!if (dislocationtype .eq. '110') then
!!    do i=1,sq
!!                !write out matrix of evasions for z=0.5*a0 plane
!!        write(21, 10), &
!!        (X( m05((1+(i-1)*sq ):(i*sq )) ) - Xi( m05((1+(i-1)*sq ):(i*sq )) ) )&
!!        *sqrt(5d-1) - &
!!        (Y( m05((1+(i-1)*sq ):(i*sq )) ) - Yi( m05((1+(i-1)*sq ):(i*sq )) ) )&
!!        *sqrt(5d-1)
!!
!!        write(22, 10), &
!!        (Y( m05((1+(i-1)*sq ):(i*sq )) ) - Yi( m05((1+(i-1)*sq ):(i*sq )) ) )&
!!        *sqrt(5d-1) + &
!!        (X( m05((1+(i-1)*sq ):(i*sq )) ) - Xi( m05((1+(i-1)*sq ):(i*sq )) ) )&
!!        *sqrt(5d-1)
!!    enddo
!!    do i=1,sqq
!!                !write out matrix of evasions for z=0     plane
!!        write(23, 10), &
!!        (X( m00((1+(i-1)*sqq):(i*sqq)) ) - Xi( m00((1+(i-1)*sqq):(i*sqq)) ) )&
!!        *sqrt(5d-1) - &
!!        (Y( m00((1+(i-1)*sqq):(i*sqq)) ) - Yi( m00((1+(i-1)*sqq):(i*sqq)) ) )&
!!        *sqrt(5d-1)
!!
!!        write(24, 10), &
!!        (Y( m00((1+(i-1)*sqq):(i*sqq)) ) - Yi( m00((1+(i-1)*sqq):(i*sqq)) ) )&
!!        *sqrt(5d-1) + &
!!        (X( m00((1+(i-1)*sqq):(i*sqq)) ) - Xi( m00((1+(i-1)*sqq):(i*sqq)) ) )&
!!        *sqrt(5d-1)
!!    enddo
!!endif
!!
!!close(21)
!!close(22)
!!close(23)
!!close(24)
!!
!!10 format (1x,62F9.4)
!!11 format (1x,153F11.7)
!!12 format (1x,153F11.7)
!!13 format (1x,153F11.7)
!!14 format (A,I2.2,A)
!
!endsubroutine
!
!
!subroutine wo_grid_highlighted
!!    use latticemod
!!    use shift_valuesmod
!!
!!
!!!variables local
!!    integer i,s, color, phi_index, radrad_index, j
!!    real phi_atom, radrad_atom, shift
!!    character (LEN = 18) filename_grided
!!    logical is_atom_to_write
!!
!!    write(filename_grided,'(A,I2.2,A)'),"VESTA_grided",vestaWOcounter,".mld"
!!    open (13, file = filename_grided)
!!s = 0
!!shift = 1.1*a0
!!
!!    do i=1,qnt
!!        is_atom_to_write = (abs(Z(i)) .lt. a0*c+shift)
!!        if (is_atom_to_write) s=s+1
!!    enddo
!!
!!    write(13, *), "grid_highlighted"
!!    write(13, *), s
!!
!!    do i = 1,bq
!!        is_atom_to_write = (abs(Z(i)) .lt. a0*c+shift)
!!        if (.not. is_atom_to_write) cycle
!!        !now lets color atoms
!!        phi_atom    = atan2(yi(i),abs(xi(i)))
!!        radrad_atom = xi(i)*xi(i) + yi(i)*yi(i)
!!
!!        do j=1,phipoints
!!            if(ux(0,j) .lt. phi_atom ) cycle
!!            phi_index=j-1
!!            exit
!!        enddo
!!        do j=1,radradpoints
!!            if(ux(j,0) .lt. radrad_atom ) cycle
!!            radrad_index=j-1
!!            exit
!!        enddo
!!        !color=5+mod(radrad_index+mod(phi_index,2)*2,7)
!!        color=5+mod(phi_index,7) !smart color(r,phi) does not represent grid good, it is too small
!!        write(13, "(F7.2, 2x, F7.2, 2x, F7.2, 2x, I2)") X(i)*5, Y(i)*5, Z(i)*5, color
!!    enddo
!!
!!    do i = bq+1,qnt
!!        is_atom_to_write = (abs(Z(i)) .lt. a0*c+shift)
!!        !if (.not. is_atom_to_write) cycle
!!
!!        if (is_atom_to_write) &
!!        write(13, "(F7.2, 2x, F7.2, 2x, F7.2, 2x, I2)") X(i)*5, Y(i)*5, Z(i)*5, 3
!!    enddo
!!
!!close(13)
!!    vestaWOcounter=vestaWOcounter+1
!!endsubroutine wo_grid_highlighted
!!
!!subroutine wo_shifts
!!    use shift_valuesmod
!!
!!
!!    integer i,j
!!
!!    open (17, file = 'ux_phi_rsquared.txt')
!!    open (18, file = 'uy_phi_rsquared.txt')
!!    open (19, file = 'uz_phi_rsquared.txt')
!!
!!    write (17,17) ux(0,1:phipoints)
!!    write (18,17) uy(0,1:phipoints)
!!    write (19,17) uz(0,1:phipoints)
!!
!!    write (17,17) ux(1:radradpoints,0)
!!    write (18,17) uy(1:radradpoints,0)
!!    write (19,17) uz(1:radradpoints,0)
!!
!!    do i=1,radradpoints
!!
!!        write(17,17) ux(i,1:phipoints)
!!        write(18,17) uy(i,1:phipoints)
!!        write(19,17) uz(i,1:phipoints)
!!
!!    enddo
!!
!!    close(17)
!!    close(18)
!!    close(19)
!!
!!17  format(100E12.4)
!!endsubroutine
!!
!!subroutine wo_deformations
!!    use deformation_valuesmod
!!
!!
!!    integer i
!!
!!    open (24, file = 'dux_dphi_phi_rsquared.txt')
!!    open (25, file = 'duy_dphi_phi_rsquared.txt')
!!    open (26, file = 'duz_dphi_phi_rsquared.txt')
!!    open (27, file = 'dux_drsquared_phi_rsquared.txt')
!!    open (28, file = 'duy_drsquared_phi_rsquared.txt')
!!    open (29, file = 'duz_drsquared_phi_rsquared.txt')
!!
!!    write (24,17) dux_dphi(   0,1:phipoints)
!!    write (25,17) duy_dphi(   0,1:phipoints)
!!    write (26,17) duz_dphi(   0,1:phipoints)
!!    write (27,17) dux_dradrad(0,1:phipoints)
!!    write (28,17) duy_dradrad(0,1:phipoints)
!!    write (29,17) duz_dradrad(0,1:phipoints)
!!
!!    write (24,17) dux_dphi(   1:radradpoints,0)
!!    write (25,17) duy_dphi(   1:radradpoints,0)
!!    write (26,17) duz_dphi(   1:radradpoints,0)
!!    write (27,17) dux_dradrad(1:radradpoints,0)
!!    write (28,17) duy_dradrad(1:radradpoints,0)
!!    write (29,17) duz_dradrad(1:radradpoints,0)
!!
!!    do i=1,radradpoints
!!
!!        write (24,17) dux_dphi(   i,1:phipoints)
!!        write (25,17) duy_dphi(   i,1:phipoints)
!!        write (26,17) duz_dphi(   i,1:phipoints)
!!        write (27,17) dux_dradrad(i,1:phipoints)
!!        write (28,17) duy_dradrad(i,1:phipoints)
!!        write (29,17) duz_dradrad(i,1:phipoints)
!!
!!    enddo
!!
!!    close(24)
!!    close(25)
!!    close(26)
!!    close(27)
!!    close(28)
!!    close(29)
!!
!!17  format(100E12.4)
!!endsubroutine
!!
!!subroutine wo_summary_and_inputs
!!
!!    use latticemod
!!    use addition_function_mod
!!
!!    open(33,file = 'input_information.txt' )
!!
!!    write(33,*) " "
!!
!!    write(33,30)" a = ",a," ; b = ",b," ; c = ",c
!!    write(33,*) " atoms oberved    :      ", observedatoms(0)
!!    write(33,*) " atoms in main cell      ", bq
!!    write(33,*) " dislocation type is     ", dislocationtype
!!
!!    write(33,*) " different values of r^2 ", radradpoints
!!    write(33,*) " different values of phi ",    phipoints
!!    write(33,*) " multiplicative constant ",  denominator
!!
!!    close(33)
!!
!!
!!30  format(1X,A,I3,1X,A,I3,1X,A,I3)
!endsubroutine

!subroutine      conjugated_relaxation
!    use comp_parameters_mod
!    use positions_mod
!    use cgm_storage_mod
!
!    integer i_passages,i_atoms,i_directions
!    real(8) gamma_curr,gamma_prev,gamma_sum
!    real(8) scalar_curr,scalar_prev,scalar_summ
!    real(8) e_curr,e_prev,D_denominator,shift_magnitude
!    real(8) fx,fy,fz,grad_max,grad_curr
!    logical need_to_rescale_total_gradient
!!real(8), dimension(1:3,atoms_max_array) ::  F_prev,F_curr,D_prev,D_curr gamma_
!
!    print*,"there will be ",total_relax_passages,&
!    "conjugated gradient descent relaxation stages"
!    maxgrad=0d0
!    grad_curr=0d0
!    scalar_summ =   0d0
!    scalar_prev =   1d0
!        !@ obtain current force
!        do  i_atoms=first_relaxable,last__relaxable
!            call energygradient_pw(i_atoms,fx,fy,fz)!gradient of a single atom
!            F_curr(1,i_atoms)=fx
!            F_curr(2,i_atoms)=fy
!            F_curr(3,i_atoms)=fz
!            grad_curr=(fx*fx+fy*fy+fz*fz)
!            scalar_summ=scalar_summ+grad_curr
!            if(maxgrad.lt.grad_curr)then
!                maxgrad = grad_curr
!                print 211,maxgrad , " is maximum atomic gradient"
!            endif
!        enddo
!        scalar_curr=scalar_summ
!        !@ current force is obtained
!        !@ check current force if it is small enough to finalize procedure or rescale
!        if  (scalar_curr.lt.system_gradient_small_enough**2) then
!            print*," total gradient have been weighed ... and found wanting. more precisely, less than "&
!            ,system_gradient_small_enough
!            return
!        endif
!        need_to_rescale_total_gradient=.false.
!        !D_denominator=(5d-1/gamma_curr)
!        if  (scalar_curr.gt.sqrt(dfloat(last__relaxable-first_relaxable+1))) then
!            print*," total gradient have been weighed ... and found big (see it further). rescaling ensues."&
!            ,scalar_curr
!            need_to_rescale_total_gradient=.true.
!            D_denominator=(7d-1/scalar_curr)
!        endif
!        !gamma_curr=scalar_curr/scalar_prev
!
!        D_curr(1:3,first_relaxable:last__relaxable)=&
!        F_curr(1:3,first_relaxable:last__relaxable)!+&
!        !D_prev(1:3,first_relaxable:last__relaxable)*&
!        !gamma_curr !not needed at that step
!        !shift_magnitude=norm2(D_curr(1:3,first_relaxable:last__relaxable))
!        shift_magnitude=sqrt(abs(scalar_prev*D_denominator))
!        if(need_to_rescale_total_gradient)then
!            D_curr(1:3,first_relaxable:last__relaxable)=&
!            D_curr(1:3,first_relaxable:last__relaxable)*D_denominator
!            shift_magnitude=sqrt(scalar_prev*D_denominator)
!        endif
!        scalar_prev=scalar_curr
!        D_prev(1:3,first_relaxable:last__relaxable)=D_curr(1:3,first_relaxable:last__relaxable)
!        !F_prev(1:3,first_relaxable:last__relaxable)=F_curr(1:3,first_relaxable:last__relaxable) !do i really need prev F ?
!        !minimize along d_curr
!        call    energy_of_system(e_curr)
!        do i_directions=1,cgm_direction_switch
!            e_prev  =e_curr
!            R_curr(1:3,first_relaxable:last__relaxable)=&
!            R_curr(1:3,first_relaxable:last__relaxable)+&
!            D_curr(1:3,first_relaxable:last__relaxable)
!            call    set_positions_of_periodic_atoms
!            call    energy_of_system(e_curr)
!            if(e_curr.gt.e_prev) then
!                D_curr(1:3,first_relaxable:last__relaxable)=&
!                D_curr(1:3,first_relaxable:last__relaxable)*(-cgm_denominator)
!                !print"(A,$)","|reversing direction and dividing it by two | "
!                shift_magnitude=shift_magnitude*abs(cgm_denominator)
!                if(shift_magnitude.lt.system_gradient_small_enough)then !maybe another criteria
!                    call    wo_xyz_snapshot!@   write the image of trajectory
!                    exit
!                endif
!            endif
!            call    wo_xyz_snapshot!@   write the image of trajectory
!        enddo!@ minimization along D_curr is done
!
!
!
!    main_cgd_passages_loop : &
!    do i_passages=2,total_relax_passages !@ now do it in a cycle
!        maxgrad=0d0
!        grad_curr=0d0
!
!        scalar_summ =   0d0
!        !@ obtain current force
!        do  i_atoms=first_relaxable,last__relaxable
!            call energygradient_pw(i_atoms,fx,fy,fz)!gradient of a single atom
!            F_curr(1,i_atoms)=fx
!            F_curr(2,i_atoms)=fy
!            F_curr(3,i_atoms)=fz
!            grad_curr=(fx*fx+fy*fy+fz*fz)
!            scalar_summ=scalar_summ+grad_curr
!            if(maxgrad.lt.grad_curr) then
!                maxgrad = grad_curr
!                print 211,maxgrad , " is maximum atomic gradient"
!            endif
!        enddo
!        scalar_curr=scalar_summ
!        !@ current force is obtained
!        !@ check current force if it is small enough to finalize procedure or rescale
!        if  (scalar_curr.lt.system_gradient_small_enough**2) then
!            print*," total gradient have been weighed ... and found wanting. more precisely, less than "&
!            ,system_gradient_small_enough
!            exit main_cgd_passages_loop
!        endif
!        need_to_rescale_total_gradient=.false.
!        if  (scalar_curr.gt.sqrt(dfloat(last__relaxable-first_relaxable+1))) then
!            print*," total gradient have been weighed ... and found big (see it further). rescaling ensues."&
!            ,scalar_curr
!            need_to_rescale_total_gradient=.true.
!            D_denominator=(2d0/scalar_curr)
!        endif
!        gamma_curr=scalar_curr/scalar_prev
!        !@ check is done, now find the direction of descending
!        D_curr(1:3,first_relaxable:last__relaxable)=&
!        F_curr(1:3,first_relaxable:last__relaxable)+&
!        D_prev(1:3,first_relaxable:last__relaxable)*&
!        gamma_curr
!
!        shift_magnitude=sqrt(scalar_prev)
!
!        if(need_to_rescale_total_gradient)then
!            D_curr(1:3,first_relaxable:last__relaxable)=&
!            D_curr(1:3,first_relaxable:last__relaxable)*D_denominator
!            shift_magnitude=sqrt(scalar_prev*D_denominator)
!        endif
!        scalar_prev=scalar_curr
!
!        D_prev(1:3,first_relaxable:last__relaxable)=D_curr(1:3,first_relaxable:last__relaxable) !@ save the directions
!        !@ minimize along D_curr
!        call    energy_of_system(e_curr)
!        do i_directions=1,cgm_direction_switch
!
!            e_prev  =e_curr
!            R_curr(1:3,first_relaxable:last__relaxable)=&
!            R_curr(1:3,first_relaxable:last__relaxable)+&
!            D_curr(1:3,first_relaxable:last__relaxable)
!            call    set_positions_of_periodic_atoms
!            call    energy_of_system(e_curr)
!            !PRINT*, I_PASSAGES, I_DIRECTIONS
!            if(e_curr.gt.e_prev) then
!                D_curr(1:3,first_relaxable:last__relaxable)=&
!                D_curr(1:3,first_relaxable:last__relaxable)*(-cgm_denominator)
!                !print"(A,$)","|reversing direction and dividing it by two | "
!                shift_magnitude=shift_magnitude*abs(cgm_denominator)
!                if(shift_magnitude.lt.system_gradient_small_enough)then !maybe another criteria
!                    call    wo_xyz_snapshot!@   write the image of trajectory
!                    exit
!                endif
!            endif
!            call    wo_xyz_snapshot!@   write the image of trajectory
!        enddo!@ minimization along D_curr is done
!
!    enddo   main_cgd_passages_loop
!211 format(1x,E11.4,A,$)
!    print*," "
!    print*,"we are done with conjugated descending"
!endsubroutine   conjugated_relaxation
