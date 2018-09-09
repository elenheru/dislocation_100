!MS$DECLARE
! UPPERCASE - FOR DEBUGGING(except for main program description)
! to create a directory !: call execute_command_line('mkdir directoryname',WAIT=.TRUE.)
PROGRAM    DISLOCATION_CORE_100_SPLIT

    REAL(8) ANC1,ANC2,ANC3,EATOM
    INTEGER CENTERAL_ATOM_NUMBER
    !;
    call    set_seed_of_random_numbers  !@idnk whether it is needed
    call    spawn_bcc_rectangular_100   !@place the rectangular formed crystal of bcc structure
    !CALL    PART_THREE_ZONES            !@MAIN CELL, PERIODIC BC, FIXED BC
    call    dislocate_100_split         !@dislocate before parting, to erase garbage atoms
    !call    wo_direct
    call    part_three_zones_cylinder            !@main cell, periodic bc, fixed bc
    !call    wo_zoned ; stop "check partition"

    !call    sort_atoms_by_distance      !@relaxation of centeral atoms is more prior
    call    cat_symmetric_pairs         !@along Z there are periodic boundary conditions
    !!call    cast_hirth_shifts
    CALL    ENERGY_OF_SYSTEM(ANC3)
    call    cat_matrixes
    call    cast_gusev_shifts
    call    wo_matrixes_100
    !call    cat_verlet_list
    call    wo_zoned
    !stop"layer check"

!!    PRINT*, "1", EATOM(CENTERAL_ATOM_NUMBER())
!!    CALL    RENEW_VERLET_LIST_ATOM(CENTERAL_ATOM_NUMBER())
!!    PRINT*, "2", EATOM(CENTERAL_ATOM_NUMBER())
!!    CALL    RENEW_VERLET_LIST_NEIGHBORS(CENTERAL_ATOM_NUMBER())
!!    PRINT*, "3", EATOM(CENTERAL_ATOM_NUMBER())
!!    CALL    ENERGYGRADIENT_PW(CENTERAL_ATOM_NUMBER(),ANC1,ANC2,ANC3)
    !CALL    RELAX_ATOM(CENTERAL_ATOM_NUMBER())
    !call    cat_verlet_list             !@neibors lists for fast energy calculation
    call    energy_of_system(anc3)
    !CALL    WO_XYZ_SNAPSHOT
    STOP "matrices"
    !call    proceed_core_relaxation
    !call    conjugated_relaxation
    call    proceed_convergation_burgers_poisson
    !call    proceed_core_relaxation     !@relax every atom considering symmetry
    !call    conjugated_relaxation
    !CALL    ENERGY_OF_SYSTEM(ANC3)
    call    wo_summary_and_results(anc3)      !@results and summary to console and files
    !CALL    ENERGY_OF_SYSTEM(ANC3)
    !call    wo_direct                      !@write out direct VESTA-formatted image of system
    call    wo_matrixes_100
    call    wo_zoned
ENDPROGRAM      DISLOCATION_CORE_100_SPLIT

subroutine      spawn_bcc_rectangular_100
    use positions_mod
    use comp_parameters_mod
    use phys_parameters_mod
    real(8) a_half
    integer cnt_x, cnt_y, cnt_z, cnt_q

    a_half  = a0*5d-1
    cnt_q = 0
    do cnt_x=-x_layers,x_layers
        if (cnt_q + 3 .ge. atoms_max_array) stop "array overflow : too much of atoms"
        do cnt_y=-y_layers,y_layers
            do cnt_z=-z_layers,z_layers
                cnt_q = cnt_q + 1
                R_curr(1,cnt_q) = cnt_x * a0
                R_curr(2,cnt_q) = cnt_y * a0
                R_curr(3,cnt_q) = cnt_z * a0
                if(cnt_x.eq.x_layers) cycle
                if(cnt_y.eq.y_layers) cycle
                if(cnt_z.eq.z_layers) cycle
                cnt_q = cnt_q + 1
                R_curr(1,cnt_q) = cnt_x * a0 + a_half
                R_curr(2,cnt_q) = cnt_y * a0 + a_half
                R_curr(3,cnt_q) = cnt_z * a0 + a_half
            enddo
        enddo
    enddo
    atoms__in_total = cnt_q
    R_perf=R_curr
    write(*,*) "BCC rectangular cell had been placed. Total atoms: ", cnt_q
    write(*,106) " cell size:_",(x_layers*2+1)," x",(y_layers*2+1)," x",(z_layers*2+1),&
    " _;  a_0 lattice is ", a0
106 format  (3(A,I3.1),A,F8.4)
endsubroutine   spawn_bcc_rectangular_100

subroutine      part_three_zones_cylinder !find indexes for main comp cell and
    use positions_mod
    use comp_parameters_mod
    use phys_parameters_mod
    use anisotropy_mod
    real(8)             :: R_swap(1:3,atoms_max_array)!swap massive
    real(8)             :: xedge, yedge, zedge!explicit values for geometry conditions check
    real(8),parameter   :: okr=1d-3           !small epsilon
    integer i, s
    !variables end
    !digital_wth=(0.5*wth+okr)*a0    !value of thickness of second zone [А]
    xedge=a0*(x_layers-x_edges)-okr !explicit value for zones parting at large  X
    yedge=a0*(y_layers-y_edges)-okr !explicit value for zones parting at large  X
    zedge=a0*(z_layers-z_edges)-okr
    R_swap=R_curr

    s=1;first_relaxable=s

    do i=1,atoms__in_total!body atoms
        !if (abs(R_swap(1,i)) .gt. xedge) cycle
        !if (abs(R_swap(2,i)) .gt. yedge) cycle

        if (norm2(R_swap(1:2,i)) .gt. xedge) cycle
        if (norm2(R_swap(1:2,i)) .gt. yedge) cycle

        if (abs(R_swap(3,i)) .gt. zedge) cycle

        R_curr(1:3,s)=R_swap(1:3,i)
        s=s+1
    enddo
    last__relaxable=s-1;first_z__period=s
    if (z_edges.gt.0) then
        do i=1,atoms__in_total!periodic atoms
            if (norm2(R_swap(1:2,i)) .gt. xedge) cycle
            if (norm2(R_swap(1:2,i)) .gt. yedge) cycle
            if (abs(R_swap(3,i)) .lt. zedge) cycle

            R_curr(1:3,s)=R_swap(1:3,i)
            s=s+1
        enddo
        last__z__period=s-1
    else
        last__z__period=s
        print*,"periodic atoms are absent, no symmetric pairs are present "
    endif
    first______wall=s

    do i=1,atoms__in_total!wall atoms
        if ((norm2(R_swap(1:2,i)) .lt. yedge).and.&
            (norm2(R_swap(1:2,i)) .lt. xedge))cycle
        if (R_swap(1,i)       .gt. a0*(1+(x_layers+y_layers+z_layers)*3)) cycle !dislocated
        R_curr(1:3,s)=R_swap(1:3,i)
        s=s+1
    enddo
    last_______wall=s-1;atoms__in_total=s-1

    if  (first_relaxable.ge.last__relaxable) then
        print*,first_relaxable,last__relaxable,"relaxable"
        print*,first_z__period,last__z__period,"periodic"
        print*,first______wall,last_______wall,"wall"
        stop "no atoms in body part"
    endif
    if  (first_z__period.ge.last__z__period) then
        print*,first_relaxable,last__relaxable,"relaxable"
        print*,first_z__period,last__z__period,"periodic"
        print*,first______wall,last_______wall,"wall"
        stop "no atoms in periodic part"
    endif
    if  (first______wall.ge.last_______wall) then

        print*,first_relaxable,last__relaxable,"relaxable"
        print*,first_z__period,last__z__period,"periodic"
        print*,first______wall,last_______wall,"wall"
        if (z_edges.gt.0) stop "no atoms in wall part"
    endif

        print*," mind that 3 layers of 100 type planes are not enough for Ackland."
    R_perf=R_curr
    s=0
    do i=first_relaxable,last__relaxable !find out anisotropy atoms

!        if ((abs(R_curr(2,i)) .lt. yedge-a0*26d-1).and. &
!            (abs(R_curr(1,i)) .lt. xedge-a0*26d-1)) cycle
!        if (abs(R_curr(1,i))  .gt. xedge-a0*19d-1) cycle
!        if (abs(R_curr(2,i))  .gt. yedge-a0*19d-1) cycle
!        if ((abs(R_curr(1,i)) .lt. a0*1d-3 )) cycle!extraplane
!        if ((   (R_curr(3,i)) .lt.-a0*1d-3 )) cycle!not all body, but z=0    area
!        if ((   (R_curr(3,i)) .gt. a0*6d-1 )) cycle!not all body, but z=a0/2 area
        if ((norm2(R_curr(1:2,i)) .lt. yedge-a0*26d-1).and. &
            (norm2(R_curr(1:2,i)) .lt. xedge-a0*26d-1)) cycle

        if (norm2(R_curr(1:2,i))  .gt. xedge-a0*19d-1) cycle
        if (norm2(R_curr(1:2,i))  .gt. yedge-a0*19d-1) cycle

        if ((abs(R_curr(1,i)) .lt. a0*1d-3 )) cycle!extraplane
        if ((   (R_curr(3,i)) .lt.-a0*1d-3 )) cycle!not all body, but z=0    area
        if ((   (R_curr(3,i)) .gt. a0*6d-1 )) cycle!not all body, but z=a0/2 area
        s=s+1
        layer_list(s)=i
    enddo
    max_layer_list_index = s
        print 112," from  ",first_relaxable," to ",last__relaxable," are relaxable"
        print 112," from  ",first_z__period," to ",last__z__period," are periodic"
        print 112," from  ",first______wall," to ",last_______wall," are wall"
        print *  ," total ",atoms__in_total
    112 format  (A,I6.1,A,I6.1,A)
        print *  ," anisotropic sampled atoms ",max_layer_list_index
        print *  ," zones are __cylindric "

endsubroutine   part_three_zones_cylinder

subroutine      part_three_zones_rectangle !find indexes for main comp cell and
    use positions_mod
    use comp_parameters_mod
    use phys_parameters_mod
    use anisotropy_mod
    real(8)             :: R_swap(1:3,atoms_max_array)!swap massive
    real(8)             :: xedge, yedge, zedge!explicit values for geometry conditions check
    real(8),parameter   :: okr=1d-3           !small epsilon
    integer i, s
    !variables end

    !digital_wth=(0.5*wth+okr)*a0    !value of thickness of second zone [А]
    xedge=a0*(x_layers-x_edges)-okr !explicit value for zones parting at large  X
    yedge=a0*(y_layers-y_edges)-okr !explicit value for zones parting at large  X
    zedge=a0*(z_layers-z_edges)-okr
    R_swap=R_curr

    s=1;first_relaxable=s

    do i=1,atoms__in_total!body atoms
        if (abs(R_swap(1,i)) .gt. xedge) cycle
        if (abs(R_swap(2,i)) .gt. yedge) cycle
        if (abs(R_swap(3,i)) .gt. zedge) cycle

        R_curr(1:3,s)=R_swap(1:3,i)
        s=s+1
    enddo
    last__relaxable=s-1;first_z__period=s
    if (z_edges.gt.0) then
        do i=1,atoms__in_total!periodic atoms
            if (abs(R_swap(1,i)) .gt. xedge) cycle
            if (abs(R_swap(2,i)) .gt. yedge) cycle
            if (abs(R_swap(3,i)) .lt. zedge) cycle

            R_curr(1:3,s)=R_swap(1:3,i)
            s=s+1
        enddo
        last__z__period=s-1
    else
        last__z__period=s
        print*,"periodic atoms are absent, no symmetric pairs are present "
    endif
    first______wall=s

    do i=1,atoms__in_total!wall atoms
        if ((abs(R_swap(2,i)) .lt. yedge).and.(abs(R_swap(1,i)) .lt. xedge))cycle
        if (R_swap(1,i)       .gt. a0*(1+(x_layers+y_layers+z_layers)*3)) cycle
        R_curr(1:3,s)=R_swap(1:3,i)
        s=s+1
    enddo
    last_______wall=s-1;atoms__in_total=s-1

    if  (first_relaxable.ge.last__relaxable) then
        print*,first_relaxable,last__relaxable,"relaxable"
        print*,first_z__period,last__z__period,"periodic"
        print*,first______wall,last_______wall,"wall"
        stop "no atoms in body part"
    endif
    if  (first_z__period.ge.last__z__period) then
        print*,first_relaxable,last__relaxable,"relaxable"
        print*,first_z__period,last__z__period,"periodic"
        print*,first______wall,last_______wall,"wall"
        stop "no atoms in periodic part"
    endif
    if  (first______wall.ge.last_______wall) then

        print*,first_relaxable,last__relaxable,"relaxable"
        print*,first_z__period,last__z__period,"periodic"
        print*,first______wall,last_______wall,"wall"
        if (z_edges.gt.0) stop "no atoms in wall part"
    endif

        print*," mind that 3 layers of 100 type planes are not enough for Ackland."
    R_perf=R_curr
    s=0
    do i=first_relaxable,last__relaxable !find out anisotropy atoms
        if ((abs(R_curr(2,i)) .lt. yedge-a0*26d-1).and. &
            (abs(R_curr(1,i)) .lt. xedge-a0*26d-1)) cycle

        if (abs(R_curr(1,i))  .gt. xedge-a0*19d-1) cycle
        if (abs(R_curr(2,i))  .gt. yedge-a0*19d-1) cycle

        if ((abs(R_curr(1,i)) .lt. a0*1d-3 )) cycle!extraplane
        if ((   (R_curr(3,i)) .lt.-a0*1d-3 )) cycle!not all body, but z=0    area
        if ((   (R_curr(3,i)) .gt. a0*6d-1 )) cycle!not all body, but z=a0/2 area
        s=s+1
        layer_list(s)=i
    enddo
    max_layer_list_index = s
        print 112," from  ",first_relaxable," to ",last__relaxable," are relaxable"
        print 112," from  ",first_z__period," to ",last__z__period," are periodic"
        print 112," from  ",first______wall," to ",last_______wall," are wall"
        print *  ," total ",atoms__in_total
    112 format  (A,I6.1,A,I6.1,A)
        print *  ," anisotropic sampled atoms ",max_layer_list_index
        print *  ," zones are rectangular "
endsubroutine part_three_zones_rectangle

subroutine      dislocate_100_split
    use phys_parameters_mod
    use positions_mod
    use comp_parameters_mod
    integer i,s
    real(8),parameter   :: okr=1d-3           !small epsilon
    s=0
    do i=1,atoms__in_total
        if ((R_curr(2,i)) .lt. 1d-1) cycle
        if (abs(abs(R_curr(1,i))-a0*5d-1) .gt. okr) cycle
        R_curr(1:3,i)=a0*((x_layers+y_layers+z_layers)*3 +5)+R_curr(1:3,i) !move out atoms
        s=s+1
    enddo
    print*," Dislocation 100 split is made, ",s," atoms removed."
endsubroutine   dislocate_100_split

subroutine      cat_symmetric_pairs
    use positions_mod
    use comp_parameters_mod
    use phys_parameters_mod
    use periodic_conditions_mod
    integer i,j,s
    real(8),parameter   :: okr=1d-3           !small epsilon
    max_periodic_pair=last__z__period-first_z__period+1
    s=0
    if (z_edges.lt.1) print*,"periodic atoms are absent, no symmetric pairs are present "
    do j=first_z__period,last__z__period
        do i=first_relaxable,last__relaxable
            if(abs(R_curr(1,i)-R_curr(1,j)).gt.okr)cycle
            if(abs(R_curr(2,i)-R_curr(2,j)).gt.okr)cycle
            if(abs(abs(R_curr(3,i)-R_curr(3,j)) -(2*(z_layers-z_edges)-1)*a0 ).gt.okr) cycle
            s=s+1
            periodic_and_body(1,s)=j
            periodic_and_body(2,s)=i
            !PRINT*,"PAIRED ATOMS",I,J,R_CURR(3,I),R_CURR(3,J)
            exit
        enddo
    enddo
    print*,"periodic atoms are caught:",s,"atoms."
    IF(S.LT.MAX_PERIODIC_PAIR) PRINT*,"CANT FIND PERIODIC ATOMS",S,MAX_PERIODIC_PAIR
endsubroutine   cat_symmetric_pairs

subroutine      cat_verlet_list
    use positions_mod
    use interaction_mod
    use comp_parameters_mod
    integer i,j
!    INTEGER :: I_ATOMS
    real(8) dist,vector_dist(3)
    verlet_list = 0
    distan_list = 0d0
    DIST=1D5!initialised
    vector_dist=0d0
!        DO I_ATOMS = 1,ATOMS__IN_TOTAL
!            IF(ISNAN(R_CURR(1,I_ATOMS)))THEN
!            PRINT*,"2 CAME ACROSS WITH NAN COORDINATES WHILE {get_and_apply_coordinates}; ATOM#",I_ATOMS
!            ENDIF
!            IF(ISNAN(R_CURR(2,I_ATOMS)))THEN
!            PRINT*,"2 CAME ACROSS WITH NAN COORDINATES WHILE {get_and_apply_coordinates}; ATOM#",I_ATOMS
!            ENDIF
!            IF(ISNAN(R_CURR(3,I_ATOMS)))THEN
!            PRINT*,"2 CAME ACROSS WITH NAN COORDINATES WHILE {get_and_apply_coordinates}; ATOM#",I_ATOMS
!            ENDIF
!        ENDDO
    do  i=first_relaxable,last__relaxable
    !DO  I=1,ATOMS__IN_TOTAL
!            IF (ISNAN(SUM(R_CURR(1:3,I)))) THEN
!                PRINT*,"POSITION IS NONUMBER,WHILE CAT VERLET LIST"
!                PRINT*,"RADIUS VECTOR IS"
!                PRINT*,R_CURR(1:3,I)
!                PRINT*,"ATOM  IS  #I    ; I= "
!                PRINT*,I
!                STOP "CAT VERLET LIST ERROR REPORT 0 IS OVER."
!            ENDIF
        do  j=1,atoms__in_total
            if(i.eq.j)cycle
            vector_dist(1)=R_curr(1,i)-R_curr(1,j)
            if (vector_dist(1) .gt. cutoff) cycle
            if (vector_dist(1) .lt.-cutoff) cycle
            vector_dist(2)=R_curr(2,i)-R_curr(2,j)
            if (vector_dist(2) .gt. cutoff) cycle
            if (vector_dist(2) .lt.-cutoff) cycle
            vector_dist(3)=R_curr(3,i)-R_curr(3,j)
            if (vector_dist(3) .gt. cutoff) cycle
            if (vector_dist(3) .lt.-cutoff) cycle
            dist = norm2(vector_dist)!
            if (dist .gt. cutoff) cycle
!            IF (ISNAN(DIST)) THEN
!                PRINT*,"POSITION ar noNUMBER,WHILE CAT VERLET LIST"
!                PRINT*,"DISTANCE IS NAN, WHILE CAT VERLET LIST"
!                PRINT*,"DISTANCE VECTOR IS"
!                PRINT*,VECTOR_DIST(1:3)
!                PRINT*,"ATOMS ARE #I #J ; I= J="
!                PRINT*,I,J
!                STOP "CAT VERLET LIST ERROR REPORT 1 IS OVER."
!            ENDIF
!            IF ((DIST).le. 1d-4) THEN
!                PRINT*,"DISTANCE IS NEGATIVE, WHILE CAT VERLET LIST"
!                PRINT*,"DISTANCE VECTOR IS"
!                PRINT*,VECTOR_DIST(1:3)
!                PRINT*,"ATOMS ARE #I #J ; I= J="
!                PRINT*,I,J
!                PRINT*,R_CURR(1:3,I)
!                PRINT*,R_CURR(1:3,J)
!                PRINT*,"_____________"
!                !STOP "CAT VERLET LIST ERROR REPORT 2 IS OVER."
!            ENDIF
            verlet_list(0,i) = verlet_list(0,i) + 1
            verlet_list(verlet_list(0,i),i)=j
            distan_list(verlet_list(0,i),i)=dist
        enddo
    enddo
    !print*, "Verlet's list and distances list are obtained. Fast."
endsubroutine   cat_verlet_list

subroutine      renew_verlet_list_atom(na) !particular
    use positions_mod
    use interaction_mod
    use comp_parameters_mod

    integer i,j,na
    real(8) dist,vector_dist(3)

    verlet_list(:,na) = 0
    distan_list(:,na) = 0d0
    DIST=1D5
    vector_dist=0d0
    i=na
    do  j=1,atoms__in_total
        if(i.eq.j)cycle
        vector_dist(1)=R_curr(1,i)-R_curr(1,j)
        if (vector_dist(1) .gt. cutoff) cycle
        if (vector_dist(1) .lt.-cutoff) cycle
        vector_dist(2)=R_curr(2,i)-R_curr(2,j)
        if (vector_dist(2) .gt. cutoff) cycle
        if (vector_dist(2) .lt.-cutoff) cycle
        vector_dist(3)=R_curr(3,i)-R_curr(3,j)
        if (vector_dist(3) .gt. cutoff) cycle
        if (vector_dist(3) .lt.-cutoff) cycle
        dist = norm2(vector_dist)
        if (dist .gt. cutoff) cycle
        verlet_list(0,i) = verlet_list(0,i) + 1
        verlet_list(verlet_list(0,i),i)=j
        distan_list(verlet_list(0,i),i)=dist
    enddo
    !print*, "Verlet's list and distances list are obtained. Fast."
endsubroutine   renew_verlet_list_atom

subroutine      renew_verlet_list_neighbors(na) !particular
    use positions_mod
    use interaction_mod
    use comp_parameters_mod

    integer i,j,k,na
    real(8) dist,vector_dist(3)

    DIST=1D5
    vector_dist=0d0
    do  k=1,verlet_list(0,na)
        i=verlet_list(k,na)
        verlet_list(0,i) = 0
        distan_list(:,i) = 0d0
        do  j=1,atoms__in_total
            if(i.eq.j)cycle
            vector_dist(1)=R_curr(1,i)-R_curr(1,j)
            if (vector_dist(1) .gt. cutoff) cycle
            if (vector_dist(1) .lt.-cutoff) cycle
            vector_dist(2)=R_curr(2,i)-R_curr(2,j)
            if (vector_dist(2) .gt. cutoff) cycle
            if (vector_dist(2) .lt.-cutoff) cycle
            vector_dist(3)=R_curr(3,i)-R_curr(3,j)
            if (vector_dist(3) .gt. cutoff) cycle
            if (vector_dist(3) .lt.-cutoff) cycle
            dist = norm2(vector_dist)
            if (dist .gt. cutoff) cycle
            verlet_list(0,i) = verlet_list(0,i) + 1
            verlet_list(verlet_list(0,i),i)=j
            distan_list(verlet_list(0,i),i)=dist
        enddo
    enddo
    !print*, "Verlet's list and distances list are obtained. Fast.N"
endsubroutine   renew_verlet_list_neighbors

subroutine      cast_hirth_shifts
    use phys_parameters_mod
    use positions_mod
    use comp_parameters_mod
    integer i
    real(8),parameter   :: okr=1d-3           !small epsilon
    real(8) x,y,z, ux,uy,uz
    print*, "casting Hirth shifts: b::",burgers," nu: ",poisson
    print*, "but Hirth & co. field is not working properly"
    print*, "but Hirth & co. field is not working properly"
    print*, "but Hirth & co. field is not working properly"
!    ux=0d0;uy=0d0;uz=0d0
    do i=1,atoms__in_total!first______wall,last_______wall
        ux=0d0;uy=0d0;uz=0d0
        x=R_perf(1,i); y=R_perf(2,i); z=R_perf(3,i)
        if(abs(x).lt.okr)cycle !atoms on extraplane

        if((x .gt.0d0).and.(y .gt.0d0))then
!            ux= burgers/(pi*2d0)*(&
!            (atan2(y,x) + pi*( 05d-1) + x*y*5d-1/(1-poisson)/(x*x+y*y)) )
            ux= burgers/(pi*2d0)*(&
            (atan(y/x) + pi*( 05d-1) + x*y*5d-1/(1-poisson)/(x*x+y*y)) )
        else
!            ux=burgers/(pi*2d0)*& !
!            (atan2(y,x) + pi*(-15d-1) + x*y*5d-1/(1-poisson)/(x*x+y*y))
            ux=burgers/(pi*2d0)*& !
            (atan(y/x) + pi*(-15d-1) + x*y*5d-1/(1-poisson)/(x*x+y*y))
        endif
        uy=-burgers/(pi*4d0)*(&
        (1d0-2d0*poisson)/(2d0-2d0*poisson)*log(x*x+y*y)-&
        y*y/((1-poisson)*(x*x+y*y) ) )

        R_curr(1,i)=R_perf(1,i)+ux
        R_curr(2,i)=R_perf(2,i)+uy
    enddo

endsubroutine   cast_hirth_shifts

subroutine      cast_gusev_shifts
    use phys_parameters_mod
    use positions_mod
    use comp_parameters_mod
    integer i
    real(8),parameter   :: okr=1d-3           !small epsilon
    real(8) x,y,z, ux,uy,uz
    if(core_sign.gt.0)then
        print*, "casting Gusev shifts: b::",burgers," nu: ",poisson,"__compressed core ---"
    else
        print*, "casting Gusev shifts: b::",burgers," nu: ",poisson,"decompressed core +++"
    endif

!    ux=0d0;uy=0d0;uz=0d0
    do i=1,atoms__in_total
    !do i=first______wall,last_______wall
        ux=0d0;uy=0d0;uz=0d0
        x=R_perf(1,i); y=R_perf(2,i); z=R_perf(3,i)
        if(abs(x).gt.okr) then

            if((x .gt.0d0))then
                ux= burgers/(pi*2d0)*(&
                (atan(y/x) + pi*( 05d-1) + x*y*5d-1/(1-poisson)/(x*x+y*y)) )
            else
                ux=burgers/(pi*2d0)*& !
                (atan(y/x) + pi*(-05d-1) + x*y*5d-1/(1-poisson)/(x*x+y*y))
            endif
            uy=-burgers/(pi*4d0)*(&
            (1d0-2d0*poisson)/(2d0-2d0*poisson)*log(x*x+y*y)-&
            y*y/((1-poisson)*(x*x+y*y) ) )
        else!extraplane
            if(y.gt.0d0)then
                uy=-a0*5d-1*core_sign
                uz= a0*5d-1
                if(abs(y).gt.okr)then !the y shift must contain logarithmic part anyway
                    uy=uy+&
                        (&
                        -burgers/(pi*4d0)*(&
                        (1d0-2d0*poisson)/(2d0-2d0*poisson)*log(x*x+y*y)-&
                        y*y/((1-poisson)*(x*x+y*y) ) )&
                        )
                endif
            else
                if(abs(y).gt.okr)then !the y shift must contain logarithmic part anyway
                    uy= (&
                        -burgers/(pi*4d0)*(&
                        (1d0-2d0*poisson)/(2d0-2d0*poisson)*log(x*x+y*y)-&
                        y*y/((1-poisson)*(x*x+y*y) ) )&
                        )
                endif
            endif
        endif
        R_curr(1,i)=R_perf(1,i)+ux
        R_curr(2,i)=R_perf(2,i)+uy
        R_curr(3,i)=R_perf(3,i)+uz
    enddo

endsubroutine   cast_gusev_shifts

subroutine      proceed_core_relaxation
    !;  prorelaksirovatq atomy centra
    !;  vystavitq pozicii periodicxeskih atomov
    use comp_parameters_mod
    integer i_passages,i_atoms
    !print*,"there will be ",total_relax_passages,"relaxation passages (stages)"
    do i_passages=1,4!total_relax_passages
        print*,"relaxaton stage #",i_passages
        if  (mod(i_passages,2).eq.0) then
            do  i_atoms=first_relaxable,last__relaxable
                if  (mod(i_atoms,last__relaxable/15).eq.0) then
                    print"(I6.4,2x,$)",i_atoms
                    call    wo_xyz_snapshot!@   write the image of trajectory
                endif
                call    relax_atom(i_atoms)

            enddo
        else
            do  i_atoms=last__relaxable,first_relaxable,-1
                if  (mod(i_atoms,last__relaxable/15).eq.0) then
                    print"(I6.4,2x,$)",i_atoms
                    call    wo_xyz_snapshot!@   write the image of trajectory
                endif
                call    relax_atom(i_atoms)
                !call    wo_xyz_snapshot!@   write the image of trajectory
            enddo
        endif
        call    set_positions_of_periodic_atoms
    enddo
endsubroutine   proceed_core_relaxation

subroutine      set_positions_of_periodic_atoms
    use periodic_conditions_mod
    use positions_mod
    integer i_couple,periodic_atom,body_atom

    do i_couple=1,max_periodic_pair
        periodic_atom   =periodic_and_body(1,i_couple)
        body_atom       =periodic_and_body(2,i_couple)
        R_curr(1:2,periodic_atom)=R_curr(1:2,body_atom)
    enddo
endsubroutine   set_positions_of_periodic_atoms

subroutine      energygradient_pw(na,dux,duy,duz) !this is antigradient
    use positions_mod
    use interaction_mod
    integer i, k, na

    real(8) delX,delY,delZ,dux,duy,duz,r!,SCALING
    real(8) cos_alp, cos_bet, cos_gam
    real(8) JPotOrig_d1,dpoti

    !SCALING=25D-1
    dux = 0d0; duy = 0d0; duz = 0d0
!R_CURR(1,NA) = R_CURR(1,NA) -65D-4*scaling !TEST FOR FORCES
!R_CURR(2,NA) = R_CURR(2,NA) -65D-4*scaling !
!R_CURR(3,NA) = R_CURR(3,NA) -35D-4*scaling !
!R_CURR(1,NA) = R_CURR(1,NA) -65D-4*scaling !TEST FOR FORCES
!R_CURR(2,NA) = R_CURR(2,NA) -65D-4*scaling !
!R_CURR(3,NA) = R_CURR(3,NA) -286D-4*scaling !
!R_CURR(1,NA) = R_CURR(1,NA) -286D-2*scaling !
!CALL CAT_VERLET_LIST
!WRITE(*,*) verlet_list(0,NA),NA
!F_x = -0.5*SUM_of_k`[V_r(r_k)*x_k/r_k]
    !renew_verlet_list_atom(na)!gradiens is calcuated right after energy, so it is not needed
    do k=1,verlet_list(0,na)
        i       = verlet_list(k,na)
        delX    = R_curr(1,na)-R_curr(1,i)!+1d-30
        delY    = R_curr(2,na)-R_curr(2,i)!+1d-30
        delZ    = R_curr(3,na)-R_curr(3,i)!+1d-30
        r       = distan_list(k,na)
        cos_alp = delX/r
        cos_bet = delY/r
        cos_gam = delZ/r

        dpoti   = JPotOrig_d1(r) !dpotu(r)

        dux = dux - dpoti*cos_alp
        duy = duy - dpoti*cos_bet
        duz = duz - dpoti*cos_gam
    enddo
    dux = 5d-1*(dux) !first 0.5 is from formula
    duy = 5d-1*(duy) !see the article
    duz = 5d-1*(duz) !second is from derivative of sqrt
    !WRITE(*,*) DUX,DUY,DUZ,"GRADIENT"
    !WRITE(*,*) R_CURR(1:3,NA),"POSITION"
endsubroutine energygradient_pw

subroutine      relax_atom(na)
    use positions_mod
    use phys_parameters_mod
    integer na, i_denominate, i_step, i_retries
    integer, parameter :: max_denominate = 11, max_moves = 30, max_retries=4
    real(8), parameter :: small_move = 8d-06,denominator =1d0/o"10"! 125d-3 != 8d-07
    !small_move = 5d-07 gives about 13 correct signs
    real(8)     fx,fy,fz,xsh,ysh,zsh,normal,step
!    real(8)     x_prev,y_prev,z_prev,energy_best
    real(8)     energy_current,enew, eatom, rndm(max_retries+max_moves)
!    OPEN(11, FILE = 'GRADIATION.TXT')
    fx  =0d0; fy  =0d0; fz  =0d0
    xsh =0d0; ysh =0d0; zsh =0d0
!    ENERGY_CURRENT = EATOM(NA)
!    WRITE(*,*) "PERforming a relaxation of atom# ", NA, ENERGY_CURRENT
!    WRITE(*,*) R_CURR(1,NA),R_CURR(2,NA),R_CURR(3,NA)
!    R_CURR(1,NA)=R_CURR(1,NA)-193227D-6
!    R_CURR(2,NA)=R_CURR(2,NA)-125123D-6
!    R_CURR(3,NA)=R_CURR(3,NA)-264123D-6
!    WRITE(*,*) R_CURR(1,NA),R_CURR(2,NA),R_CURR(3,NA)
!    ENERGY_CURRENT = EATOM(NA)
!    WRITE(*,*) "PERFORming a relaxation of atom# ", NA, ENERGY_CURRENT
    call random_number(rndm)
    do i_retries=1,max_retries
        step = rndm(i_retries)*denominator !it is suitable to use denominatior there
        denomination_loop : do i_denominate = 1,max_denominate
            step = step * denominator
            do i_step = 1,max_moves
                energy_current = eatom(na)
                call energygradient_pw(na,fx,fy,fz)!gradient izmenilsia.
                normal = (fx*fx+fy*fy+fz*fz)
                !if (normal .lt. 1d-14) cycle denomination_loop
                if (normal .lt. 1d-09) return
                normal = sqrt(fx*fx+fy*fy+fz*fz)
                normal=1d0/normal

                fx=fx*normal
                fy=fy*normal
                fz=fz*normal
                xsh=fx*step
                ysh=fy*step
                zsh=fz*step
                R_curr(1,na)=R_curr(1,na)+xsh
                R_curr(2,na)=R_curr(2,na)+ysh
                R_curr(3,na)=R_curr(3,na)+zsh
                enew=eatom(na)
                if (enew .gt. energy_current) then
                    R_curr(1,na)=R_curr(1,na)-xsh*rndm(i_step)!*6d-1
                    R_curr(2,na)=R_curr(2,na)-ysh*rndm(i_step)!*6d-1
                    R_curr(3,na)=R_curr(3,na)-zsh*rndm(i_step)!*6d-1
                    call renew_verlet_list_atom(na)
!                        WRITE(*,*) "WOW",XSH,YSH,ZSH
!                        ELSE
!                        WRITE(*,*) "WWWWQEQRSW",XSH,YSH,ZSH
                endif
                !WRITE(11,*) R_curr(1,na),R_curr(2,na),R_curr(3,na)
            enddo
            if((abs(xsh)+abs(ysh)+abs(zsh)) .lt. small_move) exit!e-7 is ok
        enddo denomination_loop
    enddo
!    WRITE(*,*) R_CURR(1,NA),R_CURR(2,NA),R_CURR(3,NA)
    !CLOSE(11)
!    ENERGY_CURRENT = EATOM(NA)
!    Write(*,*) "performed  a relaxation of atom# ", na, energy_current
!    WRITE(*,*) R_CURR(1,NA),R_CURR(2,NA),R_CURR(3,NA)
    call renew_verlet_list_neighbors(na)!not sure if that is needed
!    energy_current = eatom(na)
!    write(*,*) "performed_ a relaxation of atom# ", na, energy_current
endsubroutine   relax_atom

subroutine      conjugated_relaxation
    use comp_parameters_mod
    use positions_mod
    use cgm_storage_mod

    integer i_passages,i_atoms,i_directions,grad_max_index,cgm_direction_multiplicator,thro
    real(8) backstep_value(3000)
    real(8) gamma_curr,gamma_prev,gamma_sum
    real(8) scalar_curr,scalar_prev,scalar_summ
    real(8) e_curr,e_prev,D_denominator,shift_magnitude
    real(8) fx,fy,fz,grad_curr,maxgrad
    logical need_to_rescale_total_gradient,relax_the_max
!real(8), dimension(1:3,atoms_max_array) ::  F_prev,F_curr,D_prev,D_curr gamma_
    print*,"there will be ",total_relax_passages,&
    "conjugated gradient descent relaxation stages"
    relax_the_max=.true.
    relax_the_max=.false.
    maxgrad=0d0
    grad_curr=0d0
    scalar_summ =   0d0
    scalar_prev =   1d0
    grad_max_index=0
    thro=0;call    random_number(backstep_value) !thro is to know random number
    backstep_value=backstep_value*25d-2+75d-2
        !@ obtain current force
        do  i_atoms=first_relaxable,last__relaxable
            call energygradient_pw(i_atoms,fx,fy,fz)!gradient of a single atom
            F_curr(1,i_atoms)=fx
            F_curr(2,i_atoms)=fy
            F_curr(3,i_atoms)=fz
            grad_curr=(fx*fx+fy*fy+fz*fz)
            scalar_summ=scalar_summ+grad_curr
            if(maxgrad.lt.grad_curr)then
                maxgrad = grad_curr
                grad_max_index=i_atoms
            endif
        enddo
        print 211,maxgrad , " is maximum grad,atom# ",grad_max_index," "
        scalar_curr=scalar_summ
        !@ current force is obtained
        !@ check current force if it is small enough to finalize procedure or rescale
        if  (scalar_curr.lt.system_gradient_small_enough**2) then
            print*," total gradient have been weighed ... and found wanting. more precisely, less than "&
            ,system_gradient_small_enough
            return
        endif
        need_to_rescale_total_gradient=.false.
        !D_denominator=(5d-1/gamma_curr)
        if  (scalar_curr.gt.sqrt(1d-1*dfloat(last__relaxable-first_relaxable+1))) then
            print*," total gradient have been weighed ... and found big (see it further). rescaling ensues."&
            ,scalar_curr
            need_to_rescale_total_gradient=.true.
            D_denominator=(7d-1/scalar_curr)
        endif
        !gamma_curr=scalar_curr/scalar_prev

        D_curr(1:3,first_relaxable:last__relaxable)=&
        F_curr(1:3,first_relaxable:last__relaxable)!+&
        !D_prev(1:3,first_relaxable:last__relaxable)*&
        !gamma_curr !not needed at that step
        !shift_magnitude=norm2(D_curr(1:3,first_relaxable:last__relaxable))
        shift_magnitude=sqrt(abs(scalar_prev*D_denominator))
        if(need_to_rescale_total_gradient)then
            D_curr(1:3,first_relaxable:last__relaxable)=&
            D_curr(1:3,first_relaxable:last__relaxable)*D_denominator
            shift_magnitude=sqrt(scalar_prev*D_denominator)
        endif
        scalar_prev=scalar_curr
        D_prev(1:3,first_relaxable:last__relaxable)=D_curr(1:3,first_relaxable:last__relaxable)
        !F_prev(1:3,first_relaxable:last__relaxable)=F_curr(1:3,first_relaxable:last__relaxable) !do i really need prev F ?
        !minimize along d_curr
        call    energy_of_system(e_curr)
        cgm_direction_multiplicator=1
        if(e_curr.gt.0d0)cgm_direction_multiplicator=3
        do i_directions=1,cgm_direction_switch*cgm_direction_multiplicator+(i_passages/2)
            e_prev  =e_curr
            R_curr(1:3,first_relaxable:last__relaxable)=&
            R_curr(1:3,first_relaxable:last__relaxable)+&
            D_curr(1:3,first_relaxable:last__relaxable)
            call    set_positions_of_periodic_atoms
            call    energy_of_system(e_curr)

            if((e_curr-e_prev).gt.abs(1d-3*e_prev)) then
            !this is dramatic increase of energy. something went very wrong
                call    wo_xyz_snapshot!@   write the image of trajectory
                R_curr(1:3,first_relaxable:last__relaxable)=&
                R_curr(1:3,first_relaxable:last__relaxable)-&
                D_curr(1:3,first_relaxable:last__relaxable)
            endif
            if(e_curr.gt.e_prev) then
                thro=thro+1
                D_curr(1:3,first_relaxable:last__relaxable)=&
                D_curr(1:3,first_relaxable:last__relaxable)*&
                (-cgm_denominator*backstep_value(mod(thro,1000)))
                !print"(A,$)","|reversing direction and dividing it by two | "
                shift_magnitude=shift_magnitude*abs(cgm_denominator*backstep_value(mod(thro,1000)))
                if(shift_magnitude.lt.system_gradient_small_enough)then !maybe another criteria
                    call    wo_xyz_snapshot!@   write the image of trajectory
                    print*,"Current energy is ",e_curr
                    exit
                endif
            endif
            call    wo_xyz_snapshot!@   write the image of trajectory
        enddo!@ minimization along D_curr is done


    if(relax_the_max)call    relax_atom(grad_max_index)

    main_cgd_passages_loop : &
    do i_passages=2,total_relax_passages !@ now do it in a cycle
        maxgrad=0d0
        grad_curr=0d0

        scalar_summ =   0d0
        !@ obtain current force
        do  i_atoms=first_relaxable,last__relaxable
            call energygradient_pw(i_atoms,fx,fy,fz)!gradient of a single atom
            F_curr(1,i_atoms)=fx
            F_curr(2,i_atoms)=fy
            F_curr(3,i_atoms)=fz
            grad_curr=(fx*fx+fy*fy+fz*fz)
            scalar_summ=scalar_summ+grad_curr
            if(maxgrad.lt.grad_curr) then
                maxgrad = grad_curr
                grad_max_index=i_atoms
            endif
        enddo
        print 211,maxgrad , " is maximum grad,atom# ",grad_max_index," "
        scalar_curr=scalar_summ
        !@ current force is obtained
        !@ check current force if it is small enough to finalize procedure or rescale
        if  (scalar_curr.lt.system_gradient_small_enough**2) then
            print*," total gradient have been weighed ... and found wanting. more precisely, less than "&
            ,system_gradient_small_enough
            exit main_cgd_passages_loop
        endif
        need_to_rescale_total_gradient=.false.
        if  (scalar_curr.gt.sqrt(dfloat(last__relaxable-first_relaxable+1))) then
            print*," total gradient have been weighed ... and found big (see it further). rescaling ensues."&
            ,scalar_curr
            need_to_rescale_total_gradient=.true.
            D_denominator=(2d0/scalar_curr)
        endif
        gamma_curr=scalar_curr/scalar_prev
        !@ check is done, now find the direction of descending
        D_curr(1:3,first_relaxable:last__relaxable)=&
        F_curr(1:3,first_relaxable:last__relaxable)+&
        D_prev(1:3,first_relaxable:last__relaxable)*&
        gamma_curr

        shift_magnitude=sqrt(scalar_prev)

        if(need_to_rescale_total_gradient)then
            D_curr(1:3,first_relaxable:last__relaxable)=&
            D_curr(1:3,first_relaxable:last__relaxable)*D_denominator
            shift_magnitude=sqrt(scalar_prev*D_denominator)
        endif
        scalar_prev=scalar_curr

        D_prev(1:3,first_relaxable:last__relaxable)=D_curr(1:3,first_relaxable:last__relaxable) !@ save the directions
        !@ minimize along D_curr
        call    energy_of_system(e_curr)
        cgm_direction_multiplicator=1
        if(e_curr.gt.0d0)cgm_direction_multiplicator=3
        do i_directions=1,cgm_direction_switch*cgm_direction_multiplicator+(i_passages/2)

            e_prev  =e_curr
            R_curr(1:3,first_relaxable:last__relaxable)=&
            R_curr(1:3,first_relaxable:last__relaxable)+&
            D_curr(1:3,first_relaxable:last__relaxable)
            call    set_positions_of_periodic_atoms
            call    energy_of_system(e_curr)
            !PRINT*, I_PASSAGES, I_DIRECTIONS

            if((e_curr-e_prev).gt.abs(1d-3*e_prev)) then
            !this is dramatic increase of energy. something went very wrong
                call    wo_xyz_snapshot!@   write the image of trajectory
                R_curr(1:3,first_relaxable:last__relaxable)=&
                R_curr(1:3,first_relaxable:last__relaxable)-&
                D_curr(1:3,first_relaxable:last__relaxable)
            endif

            if(e_curr.gt.e_prev) then
                thro=thro+1
                D_curr(1:3,first_relaxable:last__relaxable)=&
                D_curr(1:3,first_relaxable:last__relaxable)*&
                (-cgm_denominator*backstep_value(mod(thro,1000)))
                !print"(A,$)","|reversing direction and dividing it by two | "
                shift_magnitude=shift_magnitude*abs(cgm_denominator*backstep_value(mod(thro,1000)))
                if(shift_magnitude.lt.system_gradient_small_enough)then !maybe another criteria
                    call    wo_xyz_snapshot!@   write the image of trajectory
                    print*,"Current energy is ",e_curr
                    exit
                endif
            endif
            call    wo_xyz_snapshot!@   write the image of trajectory
        enddo!@ minimization along D_curr is done
        print*,"Current energy is ",e_curr
        if(relax_the_max)call    relax_atom(grad_max_index)
    enddo   main_cgd_passages_loop
211 format(1x,E11.4,A,I6.2,A,$)
    print*," "
    print*,"we are done with conjugated descending"
endsubroutine   conjugated_relaxation

subroutine      proceed_convergation_burgers_poisson
    use comp_parameters_mod
    !; otrelaksirovatq sistemy
    !call    conjugated_relaxation
    !; najti znacxenija burgers i poisson v zavisimosti ot ugla
    !; primenitq k uprugoj zone smehxenija sootvetstvujuhxije uglu
    !; cycle
    integer i_anisotropy
        CALL    CONJUGATED_RELAXATION2
        call    conjugated_relaxation
        CALL    CONJUGATED_RELAXATION2
        call    conjugated_relaxation
    do i_anisotropy=1,anisotropy_passages
        call    wo_matrixes_100
        print*,"anisotropy passage #",i_anisotropy," of ",anisotropy_passages
        !CALL    CONJUGATED_RELAXATION2
        call    conjugated_relaxation
        call    get_burgers_poisson_modified
        call    cast_gusev_shifts_anisotropic
    enddo
        print*, "relaxation for the final writeout"
        call    conjugated_relaxation
        call    wo_matrixes_100
endsubroutine   proceed_convergation_burgers_poisson

subroutine      energy_of_system(e_sys)
    use comp_parameters_mod
    use positions_mod
    use interaction_mod
    real(8),intent(out)::e_sys
    real(8) esum,JPotOrig,dist
    integer i, j , k

    !print"(A, $ )", "_calculating energy of system ... inside the energy_of_system subroutine ... "
    call    cat_verlet_list
    esum=0d0
    do  i=first_relaxable,last__relaxable
        do  k=1,verlet_list(0,i)
            j   =verlet_list(k,i)
            if (j.le.i) cycle
            dist=distan_list(k,i)
            esum=esum+JPotOrig(dist)
        enddo
    enddo
    e_sys=esum
    !print*, "energy of system is ",e_sys,"eV"
endsubroutine   energy_of_system

subroutine      cat_matrixes
    use positions_mod
    use matrixes_mod
    use comp_parameters_mod
    use phys_parameters_mod, only:a0

    integer our_plane((2*x_layers+1)*(2*y_layers+1)*(2)),atoms_in_plane
    integer i,s,i_x,i_y
    real(8) z_const,y_const,x_const
    logical found_a_node

    matrix_z00=0;matrix_z05=0
    !plane z=00
    z_const=a0*0d0  ; s=0
    do i=1,atoms__in_total
        if(abs(R_curr(3,i)-z_const).gt.1d-3) cycle
        s=s+1
        our_plane(s)=i
    enddo
    atoms_in_plane=s; s=0
    do i_y=y_layers,-y_layers,-1
    y_const=a0*i_y
    do i_x=-x_layers,x_layers
        x_const=a0*i_x
        found_a_node=.false.
        do i=1,atoms_in_plane
            if(((R_curr(2,our_plane(i))-y_const)*&
                (R_curr(2,our_plane(i))-y_const)+&
                (R_curr(1,our_plane(i))-x_const)*&
                (R_curr(1,our_plane(i))-x_const)).gt.1d-3) cycle
            matrix_z00(i_x,i_y)=our_plane(i)
            found_a_node=.true.
            s=s+1
            exit
        enddo
        if(.not. found_a_node) matrix_z00(i_x,i_y)=-1 !@ a special code, for approximated nodes.
    enddo!print*,"cant find node ",i_x,i_y,"z=00"
        !PRINT'(100(I5.1,1X))',MATRIX_Z00(-X_LAYERS:X_LAYERS,I_Y)
    enddo
    print*,"z=00; found nodes: ",s,". exepected to find: ",atoms_in_plane
    !plane z=05
    z_const=a0*5d-1 ; s=0
    do i=1,atoms__in_total
        if(abs(R_curr(3,i)-z_const).gt.1d-3) cycle
        s=s+1
        our_plane(s)=i
    enddo
    atoms_in_plane=s; s=0
    do i_y=y_layers-1,-y_layers,-1
    y_const=a0*i_y+a0*5d-1
    do i_x=-x_layers,x_layers-1
        x_const=a0*i_x+a0*5d-1
        found_a_node=.false.
        do i=1,atoms_in_plane
            if(((R_curr(2,our_plane(i))-y_const)*&
                (R_curr(2,our_plane(i))-y_const)+&
                (R_curr(1,our_plane(i))-x_const)*&
                (R_curr(1,our_plane(i))-x_const)).gt.1d-3) cycle
            matrix_z05(i_x,i_y)=our_plane(i)
            found_a_node=.true.
            s=s+1
            exit
        enddo
        if(.not. found_a_node) matrix_z05(i_x,i_y)=-1 !@ a special code, for approximated nodes.
    enddo!print*,"cant find node ",i_x,i_y,"z=05"
        !PRINT'(100(I5.1,1X))',MATRIX_Z05(-X_LAYERS:X_LAYERS,I_Y)
    enddo
    print*,"z=05; found nodes: ",s,". exepected to find: ",atoms_in_plane
endsubroutine   cat_matrixes

subroutine      get_burgers_poisson_modified !# 110
    use anisotropy_mod
    use positions_mod
    use comp_parameters_mod
    character (LEN = 28) filename_relaxed
    integer i,j,na,nb
    real(8) ux, uy, uz, x, y, z, phi_1,phi_2, burgers_temp, poisson_temp
    real(8) alpha,beta,gamma_,c_,d_
    !integer ::  layer_list(atoms_max_array/4)
    !integer ::  max_layer_list_index
    burgers_poisson_WOcounter=burgers_poisson_WOcounter+1
    write(filename_relaxed,'(A,I2.2,A)'),&
    "angle_burgers_poisson_",burgers_poisson_WOcounter,".txt"
    open (110, file = filename_relaxed)
    if(burgers_poisson_WOcounter.eq.1) then
    print*,"sorting layer of sample atoms by phi from -pi to pi"
        do j=1,max_layer_list_index !sorry, but this is bubblesort
        do i=1,max_layer_list_index-1
            na=layer_list(i)
            x=R_perf(1,na);     y=R_perf(2,na);   phi_1=atan2(y,x)
            nb=layer_list(i+1)
            x=R_perf(1,nb);     y=R_perf(2,nb);   phi_2=atan2(y,x)
            if(phi_1.le.phi_2) cycle
            layer_list(i)   =nb
            layer_list(i+1) =na
        enddo
        enddo
    endif
    print*,"Writing out phi Burgers Poisson dependancy #",burgers_poisson_WOcounter
    do i=1,max_layer_list_index
        na=layer_list(i)
         x=R_perf(1,na)  ;   y=R_perf(2,na)  ;   z=R_perf(3,na)  ;phi_1=atan2(y,x)
        ux=R_curr(1,na)-x;  uy=R_curr(2,na)-y;  uz=R_curr(3,na)-z
        alpha=atan(y/(x-1d-6))
        beta=pi*5d-1 ; if(x.lt.0d0)beta=-beta
        gamma_=1d0*x*y*5d-1/(x*x+y*y) !- or +
        c_=log(x*x+y*y)
        d_=1d0*y*y/(x*x+y*y)!*5d-1!- or +
        !calculate ux and uy with default poisson burgers analyticaly
        !check whats with angle,

!        poisson_temp=1d0-&
!        (2d0*uy*gamma_+ux*(d_+5d-1*c_))   /(ux*c_-2d0*uy*(alpha+beta))
!        poisson_temp=1d0-&
!        (-uy*gamma_+ux*(2*d_+c_))   /(ux*c_*2d0+uy*(alpha+beta))
        poisson_temp=1d0-&
        (-2d0*uy*gamma_+ux*(d_+5d-1*c_))   /(ux*c_+2d0*uy*(alpha+beta))

        burgers_temp=2d0*pi*(ux-2d0*uy)/&
        (alpha+beta+c_+(gamma_-d_-5d-1*c_)/(1d0-poisson_temp))

        write(110,'(12(F10.4, 1x))'),phi_1*180/pi,burgers_temp,poisson_temp!&
        !,x,y,ux,uy,alpha,beta,gamma_,c_,d_

        burgers_list(i)=burgers_temp
        poisson_list(i)=poisson_temp
        phi_____list(i)=phi_1
!        PRINT*,"% ",PHI_____LIST(I),"%"
!        PRINT*,"% ",BURGERS_LIST(I),"%"
!        PRINT*,"% ",POISSON_LIST(I),"%"
    enddo
    PHI_____LIST(MAX_LAYER_LIST_INDEX+1)=-1234567891011D-10
    close(110)
    print*,"Modified Burgers Poisson are obtained #",burgers_poisson_WOcounter
endsubroutine   get_burgers_poisson_modified

subroutine      wo_direct  !#102
    use positions_mod
    use comp_parameters_mod
    integer i
    character (LEN = 18) filename_direct

    write(filename_direct,'(A,I2.2,A)'),"VESTA_direct",vestaWOcounter,".mld"
    open (102, file = filename_direct)
    write(102, *), "direct"
    write(102, *), atoms__in_total

    do i = 1,atoms__in_total
        write(102, 152) R_curr(1:3,i)*5 , 1
    enddo
152 format ( 3(E10.3, 1x), I3.1)
    close(102)
    write(*,*)  "Writing out_direct_VESTA. Total atoms: ",atoms__in_total
    vestaWOcounter=vestaWOcounter+1
endsubroutine   wo_direct

subroutine      wo_zoned !#103
    use positions_mod
    use comp_parameters_mod
    use anisotropy_mod
    USE PERIODIC_CONDITIONS_MOD
!variables local
!    REAL(8) R_TEST(1:3)
!    INTEGER I_TEST
    integer i,j
    character (LEN = 18) filename_relaxed
    logical test_list

    write(filename_relaxed,'(A,I2.2,A)'),"VESTA__zoned",vestaWOcounter,".mld"
    open (103, file = filename_relaxed)

    write(103, *), "zoned"
    write(103, *), atoms__in_total+1!+2*3

    do i=first_relaxable,last__relaxable
        test_list=.false.
        do j=1,max_layer_list_index
            if(layer_list(j).eq.i)then
                test_list=.true.
                exit
            endif
        enddo
        if(test_list)then
            write(103, 153) R_curr(1:3,i)*5, 16 !Sulfur S is Sample
        else
            write(103, 153) R_curr(1:3,i)*5, 5  !Boron B is Body
        endif
    enddo
    do i=first_z__period,last__z__period
        write(103, 153) R_curr(1:3,i)*5, 15     !Phosphorus P is periodic
    enddo
    do i=first______wall,last_______wall
        write(103, 153) R_curr(1:3,i)*5, 74     ! Tungsten W is Wall
    enddo
    write(103, 153) 0d0,0d0,0d0, 58!@ Cerium is Center of XYZ

!        I_TEST=110
!        R_TEST(1:3)=R_CURR(1:3,PERIODIC_AND_BODY(1,I_TEST))
!        WRITE(103, 153) R_TEST(1:3)*5D0, 27!@ COBALT IS COPY
!        R_TEST(1:3)=R_CURR(1:3,PERIODIC_AND_BODY(2,I_TEST))
!        WRITE(103, 153) R_TEST(1:3)*5D0, 38!@ STRONCIUM IS SOURCE(BODY ATOM)
!
!        I_TEST=140
!        R_TEST(1:3)=R_CURR(1:3,PERIODIC_AND_BODY(1,I_TEST))
!        WRITE(103, 153) R_TEST(1:3)*5D0, 27!@ COBALT IS COPY
!        R_TEST(1:3)=R_CURR(1:3,PERIODIC_AND_BODY(2,I_TEST))
!        WRITE(103, 153) R_TEST(1:3)*5D0, 38!@ STRONCIUM IS SOURCE(BODY ATOM)
!
!        I_TEST=70
!        R_TEST(1:3)=R_CURR(1:3,PERIODIC_AND_BODY(1,I_TEST))
!        WRITE(103, 153) R_TEST(1:3)*5D0, 27!@ COBALT IS COPY
!        R_TEST(1:3)=R_CURR(1:3,PERIODIC_AND_BODY(2,I_TEST))
!        WRITE(103, 153) R_TEST(1:3)*5D0, 38!@ STRONCIUM IS SOURCE(BODY ATOM)

153 format ( 3(E10.3, 1x), I3.1)

    close(103)
    write(*,*)  "Writing out_zoned__VESTA. Total atoms: ",atoms__in_total
    vestaWOcounter=vestaWOcounter+1
endsubroutine   wo_zoned

subroutine      wo_xyz_snapshot !#104
    use positions_mod
    use comp_parameters_mod
    USE phys_parameters_mod, ONLY:A0
    integer i_atoms,S,w
    character (LEN = 24) filename_trajectory
    !write(filename_relaxed,'(A,I2.2,A)'),  "VESTA__zoned",xyz__WOcounter,".mld"
    xyz__WOcounter=xyz__WOcounter+1

    write(filename_trajectory,'(A,I5.5,A)'),"XYZ_trajectory_",xyz__WOcounter,".xyz"
    open (104, file = filename_trajectory)
!    write(104,*),last__relaxable-first_relaxable+1
!    write(104,*),"image of relaxation # ",xyz__WOcounter
!    do i_atoms=first_relaxable,last__relaxable
!        write(104,194)  "Fe ", R_curr(1:3,i_atoms)
!    enddo
        S=0
        DO I_ATOMS=FIRST_RELAXABLE,LAST__RELAXABLE
            IF (R_CURR(3,I_ATOMS).LT.      0-1d-1) CYCLE
            IF (R_CURR(3,I_ATOMS).GT.5D-1*A0+1d-1) CYCLE
            S=S+1
        ENDDO
        W=0
        DO I_ATOMS=FIRST______WALL,LAST_______WALL
            IF (R_CURR(3,I_ATOMS).LT.      0-1d-1) CYCLE
            IF (R_CURR(3,I_ATOMS).GT.5D-1*A0+1d-1) CYCLE
            W=W+1
        ENDDO

        WRITE(104,*),S+W
        WRITE(104,*),"IMAGE OF RELAXATION # ",XYZ__WOCOUNTER
        DO I_ATOMS=FIRST_RELAXABLE,LAST__RELAXABLE
            IF (R_CURR(3,I_ATOMS).LT.      0-1d-1) CYCLE
            IF (R_CURR(3,I_ATOMS).GT.5D-1*A0+1d-1) CYCLE
            WRITE(104,194)  "Fe ", R_CURR(1:3,I_ATOMS)
        ENDDO
        DO I_ATOMS=FIRST______WALL,LAST_______WALL
            IF (R_CURR(3,I_ATOMS).LT.      0-1d-1) CYCLE
            IF (R_CURR(3,I_ATOMS).GT.5D-1*A0+1d-1) CYCLE
            WRITE(104,194)  " W ", R_CURR(1:3,I_ATOMS)
        ENDDO
    if(mod(xyz__WOcounter,25).eq.0)&
    print'(A,I5.5,A)'," Writing out XYZ picture # ",xyz__WOcounter," . "
    !print'(A,I5.5,A,$)'," Writing out XYZ picture # ",xyz__WOcounter," . "
    close(104)
194 format(A2,3(3x,ES11.4))
endsubroutine   wo_xyz_snapshot

subroutine      wo_summary_and_results(e_sys_init)!#133
    use comp_parameters_mod
    use phys_parameters_mod, only:a0

    real(8) e_sys,e_sys_init
    open( 133,file = 'input_information.txt' )

    call    energy_of_system(e_sys)
    write(133,*) "Total atoms: ", atoms__in_total
    write(133,'(3(A,I3.1),A,F8.4)') &
    " cell size:_",(x_layers*2+1)," x",(y_layers*2+1)," x",(z_layers*2+1),&
    " _;  a_0 lattice is ", a0
    write(133,*) "Energy of system befor relaxation was ",e_sys_init," eV"
    write(133,*) "Energy of system after ",total_relax_passages," relaxations is ",e_sys," eV"
    close(133)
!130  format(1X,A,I3,1X,A,I3,1X,A,I3)
endsubroutine   wo_summary_and_results

subroutine      wo_matrixes_100 !## 221:229
    use positions_mod
    use matrixes_mod
    use comp_parameters_mod
    integer i_x,i_y,na,n_left,n_rigt
    character (LEN=28) matrixfilenameUX1, matrixfilenameUX0, matrixfilenameUXU
    character (LEN=28) matrixfilenameUY1, matrixfilenameUY0, matrixfilenameUYU
    character (LEN=28) matrixfilenameUZ1, matrixfilenameUZ0, matrixfilenameUZU
    real(8) ux,uy,uz
    matrixWOcounter = matrixWOcounter + 1
    print*, "=== writing the matrixes out"
!write(*,'(A,I2.2)')," === we already did it times: ", matrixWOcounter-1

    write (matrixfilenameUX1, 14) ,"Matrix-stage" ,matrixWOcounter,("-z1-Ux-" &
    // "100" // ".txt");    open (221, file = matrixfilenameUX1)
    write (matrixfilenameUY1, 14) ,"Matrix-stage" ,matrixWOcounter,("-z1-Uy-" &
    // "100" // ".txt");    open (222, file = matrixfilenameUY1)
    write (matrixfilenameUZ1, 14) ,"Matrix-stage" ,matrixWOcounter,("-z1-Uz-" &
    // "100" // ".txt");    open (223, file = matrixfilenameUZ1)
    write (matrixfilenameUX0, 14) ,"Matrix-stage" ,matrixWOcounter,("-z0-Ux-" &
    // "100" // ".txt");    open (224, file = matrixfilenameUX0)
    write (matrixfilenameUY0, 14) ,"Matrix-stage" ,matrixWOcounter,("-z0-Uy-" &
    // "100" // ".txt");    open (225, file = matrixfilenameUY0)
    write (matrixfilenameUZ0, 14) ,"Matrix-stage" ,matrixWOcounter,("-z0-Uz-" &
    // "100" // ".txt");    open (226, file = matrixfilenameUZ0)
    !united matrices
    write (matrixfilenameUXu, 14) ,"Matrix-stage" ,matrixWOcounter,("-un-Ux-" &
    // "100" // ".txt");    open (227, file = matrixfilenameUXu)
    write (matrixfilenameUYu, 14) ,"Matrix-stage" ,matrixWOcounter,("-un-Uy-" &
    // "100" // ".txt");    open (228, file = matrixfilenameUYu)
    write (matrixfilenameUZu, 14) ,"Matrix-stage" ,matrixWOcounter,("-un-Uz-" &
    // "100" // ".txt");    open (229, file = matrixfilenameUZu)

    matrix_united_ux=0d0; matrix_united_uy=0d0; matrix_united_uz=0d0

    do i_y= y_layers-1,-y_layers,-1
    do i_x=-x_layers, x_layers-1
        !write out matrix of evasions for z=0.5*a0 plane
        na=matrix_z05(i_x,i_y)
        !print'(1x,I3.3)'
        if(na.eq. 0)cycle
        if(na.eq.-1)then
            !print*,matrix_z05(i_x+1,i_y)
            !IF(I_X.LE.-X_LAYERS)    STOP"WHY1"!TO AVOID CALL TO INCORRECT ADRESS BELOW
            !IF(I_X.GE. X_LAYERS-1)  STOP"WHY2"!TO AVOID CALL TO INCORRECT ADRESS BELOW
            if(i_x.le.-x_layers)    cycle!to avoid call to incorrect adress below
            if(i_x.ge. x_layers-1)  cycle!to avoid call to incorrect adress below
            n_left=matrix_z05(i_x-1,i_y)
            n_rigt=matrix_z05(i_x+1,i_y)
            if(n_left.le.0)then
                write(221,210), (R_curr(1,n_rigt) - R_perf(1,n_rigt))
                write(222,210), (R_curr(2,n_rigt) - R_perf(2,n_rigt))
                write(223,210), (R_curr(3,n_rigt) - R_perf(3,n_rigt))
                matrix_united_ux(2*i_x+1,2*i_y+1) = (R_curr(1,n_rigt) - R_perf(1,n_rigt))
                matrix_united_uy(2*i_x+1,2*i_y+1) = (R_curr(2,n_rigt) - R_perf(2,n_rigt))
                matrix_united_uz(2*i_x+1,2*i_y+1) = (R_curr(3,n_rigt) - R_perf(3,n_rigt))
            endif
            if(n_rigt.le.0)then
                write(221,210), (R_curr(1,n_left) - R_perf(1,n_left))
                write(222,210), (R_curr(2,n_left) - R_perf(2,n_left))
                write(223,210), (R_curr(3,n_left) - R_perf(3,n_left))
                matrix_united_ux(2*i_x+1,2*i_y+1) = (R_curr(1,n_left) - R_perf(1,n_left))
                matrix_united_uy(2*i_x+1,2*i_y+1) = (R_curr(2,n_left) - R_perf(2,n_left))
                matrix_united_uz(2*i_x+1,2*i_y+1) = (R_curr(3,n_left) - R_perf(3,n_left))
            endif
            cycle
        endif
        write(221,210), R_curr(1,na) - R_perf(1,na)
        write(222,210), R_curr(2,na) - R_perf(2,na)
        write(223,210), R_curr(3,na) - R_perf(3,na)
        matrix_united_ux(2*i_x+1,2*i_y+1) = R_curr(1,na) - R_perf(1,na)
        matrix_united_uy(2*i_x+1,2*i_y+1) = R_curr(2,na) - R_perf(2,na)
        matrix_united_uz(2*i_x+1,2*i_y+1) = R_curr(3,na) - R_perf(3,na)
    enddo
        if (i_y.ne.-y_layers) then
            write(221,*)," "!finalize the line
            write(222,*)," "
            write(223,*)," "
        endif
    enddo

    do i_y= y_layers,-y_layers,-1
    do i_x=-x_layers, x_layers
        !write out matrix of evasions for z=0.0*a0 plane
        na=matrix_z00(i_x,i_y)
        if(na.eq. 0)cycle
        if(na.eq.-1)then
            !if(i_x.le.-x_layers)!stop"why3"!to avoid call to incorrect adress below
            !if(i_x.ge. x_layers)stop"why4"!to avoid call to incorrect adress below
            n_left=matrix_z00(i_x-1,i_y)
            n_rigt=matrix_z00(i_x+1,i_y)
            if(.true.) then
                write(224,210), (R_curr(1,n_left) - R_perf(1,n_left) &
                                +R_curr(1,n_rigt) - R_perf(1,n_rigt))*5d-1
                write(225,210), (R_curr(2,n_left) - R_perf(2,n_left) &
                                +R_curr(2,n_rigt) - R_perf(2,n_rigt))*5d-1
                write(226,210), (R_curr(3,n_left) - R_perf(3,n_left) &
                                +R_curr(3,n_rigt) - R_perf(3,n_rigt))*5d-1

                matrix_united_ux(2*i_x,2*i_y) = (R_curr(1,n_left) - R_perf(1,n_left) &
                                                +R_curr(1,n_rigt) - R_perf(1,n_rigt))*5d-1
                matrix_united_uy(2*i_x,2*i_y) = (R_curr(2,n_left) - R_perf(2,n_left) &
                                                +R_curr(2,n_rigt) - R_perf(2,n_rigt))*5d-1
                matrix_united_uz(2*i_x,2*i_y) = (R_curr(3,n_left) - R_perf(3,n_left) &
                                                +R_curr(3,n_rigt) - R_perf(3,n_rigt))*5d-1
            endif
!            if(b) then
!                write(224,210), (R_curr(1,n_left) - R_perf(1,n_left) &
!                                +R_curr(1,n_rigt) - R_perf(1,n_rigt))*5d-1
!                write(225,210), (R_curr(2,n_left) - R_perf(2,n_left) &
!                                +R_curr(2,n_rigt) - R_perf(2,n_rigt))*5d-1
!                write(226,210), (R_curr(3,n_left) - R_perf(3,n_left) &
!                                +R_curr(3,n_rigt) - R_perf(3,n_rigt))*5d-1
!
!                matrix_united_ux(2*i_x,2*i_y) = (R_curr(1,n_left) - R_perf(1,n_left) &
!                                                +R_curr(1,n_rigt) - R_perf(1,n_rigt))*5d-1
!                matrix_united_uy(2*i_x,2*i_y) = (R_curr(2,n_left) - R_perf(2,n_left) &
!                                                +R_curr(2,n_rigt) - R_perf(2,n_rigt))*5d-1
!                matrix_united_uz(2*i_x,2*i_y) = (R_curr(3,n_left) - R_perf(3,n_left) &
!                                                +R_curr(3,n_rigt) - R_perf(3,n_rigt))*5d-1
!            endif
            cycle
        endif
        write(224,210), R_curr(1,na) - R_perf(1,na)
        write(225,210), R_curr(2,na) - R_perf(2,na)
        write(226,210), R_curr(3,na) - R_perf(3,na)
        matrix_united_ux(2*i_x,2*i_y) = R_curr(1,na) - R_perf(1,na)
        matrix_united_uy(2*i_x,2*i_y) = R_curr(2,na) - R_perf(2,na)
        matrix_united_uz(2*i_x,2*i_y) = R_curr(3,na) - R_perf(3,na)
    enddo
        if (i_y.ne.-y_layers) then
            write(224,*)," "!finalize the line
            write(225,*)," "
            write(226,*)," "
        endif
    enddo

    do i_y=-2*y_layers, 2*y_layers
        do i_x=-2*x_layers, 2*x_layers
            if((i_x.eq.-2*x_layers).and.(i_y.eq.-2*y_layers))cycle
            if((i_x.eq.-2*x_layers).and.(i_y.eq. 2*y_layers))cycle
            if((i_x.eq. 2*x_layers).and.(i_y.eq.-2*y_layers))cycle
            if((i_x.eq. 2*x_layers).and.(i_y.eq. 2*y_layers))cycle
            if(mod(i_y+i_x,2).eq.0)then
                if(i_x*i_x.ne.1)cycle
!                if(i_x.ne.0)cycle
!                    if(i_y.eq.-2*y_layers) then
!                        matrix_united_uy(i_x  ,i_y  )=&
!                        matrix_united_uy(i_x  ,i_y-1)*15d-1-&
!                        matrix_united_uy(i_x  ,i_y-3)*05d-1
!                        matrix_united_uz(i_x  ,i_y  )=&
!                        matrix_united_uz(i_x  ,i_y-1)*15d-1-&
!                        matrix_united_uz(i_x  ,i_y-3)*05d-1
!                        cycle
!                    endif
                if(i_y.le.0)cycle
                if(i_y.ge.2*y_layers-1) then
                    matrix_united_uy(i_x  ,i_y  )=&
                    matrix_united_uy(i_x  ,i_y-1)*15d-1-&
                    matrix_united_uy(i_x  ,i_y-3)*05d-1
                    matrix_united_uz(i_x  ,i_y  )=&
                    matrix_united_uz(i_x  ,i_y-1)*15d-1-&
                    matrix_united_uz(i_x  ,i_y-3)*05d-1
                    cycle
                endif
                if(i_y.eq.0) then
                    matrix_united_uy(i_x  ,i_y  )=&
                    matrix_united_uy(i_x  ,i_y-1)
                    matrix_united_uz(i_x  ,i_y  )=&
                    matrix_united_uz(i_x  ,i_y-1)
                    cycle
                endif
                    matrix_united_uy(i_x  ,i_y  )=&
                    matrix_united_uy(i_x  ,i_y-1)*15d-1-&
                    matrix_united_uy(i_x  ,i_y-3)*05d-1
                    matrix_united_uz(i_x  ,i_y  )=&
                    matrix_united_uz(i_x  ,i_y-1)*15d-1-&
                    matrix_united_uz(i_x  ,i_y-3)*05d-1
!                matrix_united_uy(i_x  ,i_y  )=&
!                matrix_united_uy(i_x  ,i_y+1)*5d-1+&
!                matrix_united_uy(i_x  ,i_y-1)*5d-1
!                matrix_united_uz(i_x  ,i_y  )=&
!                matrix_united_uz(i_x  ,i_y+1)*5d-1+&
!                matrix_united_uz(i_x  ,i_y-1)*5d-1
                cycle
            endif
            if ((i_x.eq.-2*x_layers) .or. (i_x.eq. 2*x_layers)) then
                matrix_united_ux(i_x  ,i_y  )=&
                matrix_united_ux(i_x  ,i_y+1)*5d-1+&
                matrix_united_ux(i_x  ,i_y-1)*5d-1
                matrix_united_uy(i_x  ,i_y  )=&
                matrix_united_uy(i_x  ,i_y+1)*5d-1+&
                matrix_united_uy(i_x  ,i_y-1)*5d-1
                matrix_united_uz(i_x  ,i_y  )=&
                matrix_united_uz(i_x  ,i_y+1)*5d-1+&
                matrix_united_uz(i_x  ,i_y-1)*5d-1
                cycle
            endif
            if (i_y.eq.-2*y_layers) then
                matrix_united_ux(i_x  ,i_y  )=&
                matrix_united_ux(i_x+1,i_y  )*5d-1+&
                matrix_united_ux(i_x-1,i_y  )*5d-1
                matrix_united_uy(i_x  ,i_y  )=&
                matrix_united_uy(i_x+1,i_y  )*5d-1+&
                matrix_united_uy(i_x-1,i_y  )*5d-1
                matrix_united_uz(i_x  ,i_y  )=&
                matrix_united_uz(i_x+1,i_y  )*5d-1+&
                matrix_united_uz(i_x-1,i_y  )*5d-1
                cycle
            endif
            if (i_y.eq. 2*y_layers)  then
                if (i_x.eq.-1)  then !do not make average with plane of symmetry
                    matrix_united_ux(i_x  ,i_y  )=&
                    matrix_united_ux(i_x-1,i_y  )*15d-1-&
                    matrix_united_ux(i_x-2,i_y  )*05d-1!delta x is a/2 not a
                    matrix_united_uy(i_x  ,i_y  )=&
                    matrix_united_uy(i_x-1,i_y  )*15d-1-&
                    matrix_united_uy(i_x-2,i_y  )*05d-1
                    matrix_united_uz(i_x  ,i_y  )=&
                    matrix_united_uz(i_x-1,i_y  )*15d-1-&
                    matrix_united_uz(i_x-2,i_y  )*05d-1
                    cycle
                endif
                if (i_x.eq. 1)  then
                    matrix_united_ux(i_x  ,i_y  )=&
                    matrix_united_ux(i_x+1,i_y  )*15d-1-&
                    matrix_united_ux(i_x+3,i_y  )*05d-1 ! i_x+2 is not ready
                    matrix_united_uy(i_x  ,i_y  )=&
                    matrix_united_uy(i_x+1,i_y  )*15d-1-&
                    matrix_united_uy(i_x+3,i_y  )*05d-1
                    matrix_united_uz(i_x  ,i_y  )=&
                    matrix_united_uz(i_x+1,i_y  )*15d-1-&
                    matrix_united_uz(i_x+3,i_y  )*05d-1
                    cycle
                endif
                matrix_united_ux(i_x  ,i_y  )=&
                matrix_united_ux(i_x+1,i_y  )*5d-1+&
                matrix_united_ux(i_x-1,i_y  )*5d-1
                matrix_united_uy(i_x  ,i_y  )=&
                matrix_united_uy(i_x+1,i_y  )*5d-1+&
                matrix_united_uy(i_x-1,i_y  )*5d-1
                matrix_united_uz(i_x  ,i_y  )=&
                matrix_united_uz(i_x+1,i_y  )*5d-1+&
                matrix_united_uz(i_x-1,i_y  )*5d-1
                cycle
            endif
            if (i_x.eq.-1)  then
                matrix_united_ux(i_x  ,i_y  )= (       &
                matrix_united_ux(i_x-1,i_y  )*2d0-     &
                matrix_united_ux(i_x-2,i_y  )  )*20d-2+&
                matrix_united_ux(i_x  ,i_y+1)*40d-2+   &
                matrix_united_ux(i_x  ,i_y-1)*40d-2

                matrix_united_uy(i_x  ,i_y  )= (       &
                matrix_united_uy(i_x-1,i_y  )*2d0-     &
                matrix_united_uy(i_x-2,i_y  )  )*20d-2+&
                matrix_united_uy(i_x  ,i_y+1)*40d-2+   &
                matrix_united_uy(i_x  ,i_y-1)*40d-2

                matrix_united_uz(i_x  ,i_y  )= (       &
                matrix_united_uz(i_x-1,i_y  )*2d0-     &
                matrix_united_uz(i_x-2,i_y  )  )*20d-2+&
                matrix_united_uz(i_x  ,i_y+1)*40d-2+   &
                matrix_united_uz(i_x  ,i_y-1)*40d-2
                cycle
            endif
            if (i_x.eq. 1)  then ! i_x=3 is not ready
                matrix_united_ux(i_x  ,i_y  )= (       &
                matrix_united_ux(i_x+1,i_y  )*175d-2+     & !surely plus
                matrix_united_ux(i_x-3,i_y  )*075d-2  )*20d-2+&
                matrix_united_ux(i_x  ,i_y+1)*40d-2+   &
                matrix_united_ux(i_x  ,i_y-1)*40d-2

                matrix_united_uy(i_x  ,i_y  )= (       &
                matrix_united_uy(i_x+1,i_y  )*175d-2-     & !surely minus
                matrix_united_uy(i_x-3,i_y  )*075d-2  )*20d-2+&
                matrix_united_uy(i_x  ,i_y+1)*40d-2+   &
                matrix_united_uy(i_x  ,i_y-1)*40d-2

                matrix_united_uz(i_x  ,i_y  )= (       &
                matrix_united_uz(i_x+1,i_y  )*175d-2-     & !surely minus
                matrix_united_uz(i_x-3,i_y  )*075d-2  )*20d-2+&
                matrix_united_uz(i_x  ,i_y+1)*40d-2+   &
                matrix_united_uz(i_x  ,i_y-1)*40d-2
                cycle !1/3 is not good multiplicator
            endif
            if(i_x.eq.0) then
                if(i_y.ge.0) then
                matrix_united_ux(i_x  ,i_y  )=0d0!&
!                matrix_united_ux(i_x+1,i_y  )*25d-2+&
!                matrix_united_ux(i_x-1,i_y  )*25d-2+&
!                matrix_united_ux(i_x  ,i_y+1)*25d-2+&
!                matrix_united_ux(i_x  ,i_y-1)*25d-2
                    matrix_united_uy(i_x  ,i_y  )=&
                    matrix_united_uy(i_x  ,i_y+1)*50d-2+&
                    matrix_united_uy(i_x  ,i_y-1)*50d-2

                    matrix_united_uz(i_x  ,i_y  )=&
                    matrix_united_uz(i_x  ,i_y+1)*50d-2+&
                    matrix_united_uz(i_x  ,i_y-1)*50d-2
                    cycle
                endif
            endif
            matrix_united_ux(i_x  ,i_y  )=&
            matrix_united_ux(i_x+1,i_y  )*25d-2+&
            matrix_united_ux(i_x-1,i_y  )*25d-2+&
            matrix_united_ux(i_x  ,i_y+1)*25d-2+&
            matrix_united_ux(i_x  ,i_y-1)*25d-2

            matrix_united_uy(i_x  ,i_y  )=&
            matrix_united_uy(i_x+1,i_y  )*25d-2+&
            matrix_united_uy(i_x-1,i_y  )*25d-2+&
            matrix_united_uy(i_x  ,i_y+1)*25d-2+&
            matrix_united_uy(i_x  ,i_y-1)*25d-2

            matrix_united_uz(i_x  ,i_y  )=&
            matrix_united_uz(i_x+1,i_y  )*25d-2+&
            matrix_united_uz(i_x-1,i_y  )*25d-2+&
            matrix_united_uz(i_x  ,i_y+1)*25d-2+&
            matrix_united_uz(i_x  ,i_y-1)*25d-2
            !cycle
        enddo
    enddo

    do i_y=-2*y_layers, 2*y_layers
        do i_x=-2*x_layers, 2*x_layers
            write(227,210),matrix_united_ux(i_x  ,i_y  )
            write(228,210),matrix_united_uy(i_x  ,i_y  )
            write(229,210),matrix_united_uz(i_x  ,i_y  )
        enddo
        if (i_y.ne.2*y_layers) then
            write(227,*)," "!finalize the line
            write(228,*)," "
            write(229,*)," "
        endif
    enddo

close(221)
close(222)
close(223)
close(224)
close(225)
close(226)
close(227)
close(228)
close(229)
014 format(A,I2.2,A)
210 format(1x,E11.4,$)
endsubroutine   wo_matrixes_100

subroutine      set_seed_of_random_numbers
    integer date_values(8), feed_input(100)
    call date_and_time(VALUES=date_values)
    feed_input = (                                  &
                 (date_values(5)*60+date_values(6)  &
                 ) *60 + date_values(7)             &
                 ) *1000 + date_values(8)
    call random_seed(put=feed_input)
    write(*,*)"Random number generator got date_and_time as seed. OK"
endsubroutine   set_seed_of_random_numbers

integer function centeral_atom_number()
    !returns number of atom that is most close to orgrin(x=0,y=0,z=0)
    use positions_mod
    use comp_parameters_mod

    integer i, challenger
    real(8) dist, dist_min
    challenger = 1
    dist_min = R_curr(1,challenger)**2 + R_curr(2,challenger)**2 + R_curr(3,challenger)**2
    do i=1,atoms__in_total
        dist = R_curr(1,i)**2 + R_curr(2,i)**2 + R_curr(3,i)**2
        if(dist .lt. dist_min) then
            challenger = i
            dist_min = dist
        endif
    enddo
    centeral_atom_number = challenger
endfunction centeral_atom_number

real(8) function eatom(na)
    use positions_mod
    use interaction_mod

    integer i,j,na
    real(8) sum_pw,   dist,JPotOrig

    eatom  = 0d0
    sum_pw = 0d0
    call renew_verlet_list_atom(na)
    if ( verlet_list(0,na) .eq. 0) return !no atoms in range
    do i=1,verlet_list(0,na)
        j       = verlet_list(i,na)
        dist    = distan_list(i,na)
        sum_pw  = sum_pw + JPotOrig(dist)
    enddo
    !eatom = (5d-1)*sum_pw + embeddin(sum_ed)
    !eatom =       sum_pw + embeddin(sum_ed)
    eatom =       sum_pw !+ embeddin(sum_ed)

endfunction eatom

real(8) pure function JPotOrig(r)   !Johnson's potential,->[A], <-[eV]
    real(8),parameter :: &
    A1=-2.195976d0, A2=3.097910d0, A3=2.704060d0, A4=7.436448d0,&
    B1=-0.639230d0, B2=3.115829d0, B3=0.477871d0, B4=1.581570d0,&
    C1=-1.115035d0, C2=3.066403d0, C3=0.466892d0, C4=1.547967d0
    real(8),intent(in):: r
    real(8)q1
    if (r.gt.3.44d0) then
        JPotOrig=0.0d0
        Return
    endif
    if (r.lt.2.40d0) then
        q1=r-A2
        JPotOrig=A1*q1*q1*q1+A3*r-A4
        Return
    endif
    if ((r.ge.2.40d0).and.(r.lt.3.00d0)) then
        q1=r-B2
        JPotOrig=B1*q1*q1*q1+B3*r-B4
         Return
    endif
    if ((r.ge.3.00d0).and.(r.le.3.44d0)) then
        q1=r-C2
        JPotOrig=C1*q1*q1*q1+C3*r-C4
    end if
endfunction JPotOrig

real(8) pure function JPotOrig_d1(r)   !Johnson's potential deriv1,->[A], <-[eV/A]
    real(8),parameter :: &
    A1=-2.195976d0, A2=3.097910d0, A3=2.704060d0, &!A4=7.436448d0,&
    B1=-0.639230d0, B2=3.115829d0, B3=0.477871d0, &!B4=1.581570d0,&
    C1=-1.115035d0, C2=3.066403d0, C3=0.466892d0   !C4=1.547967d0
    real(8),intent(in):: r
    real(8)q1
    if (r.gt.3.44d0) then
        JPotOrig_d1=0.0d0
        Return
    endif
    if (r.lt.2.40d0) then
        q1=r-A2
        JPotOrig_d1=3*A1*q1*q1+A3
        Return
    endif
    if ((r.ge.2.40d0).and.(r.lt.3.00d0)) then
        q1=r-B2
        JPotOrig_d1=3*B1*q1*q1+B3
         Return
    endif
    if ((r.ge.3.00d0).and.(r.le.3.44d0)) then
        q1=r-C2
        JPotOrig_d1=3*C1*q1*q1+C3
    end if
endfunction JPotOrig_d1

!_________________________________________________________________

subroutine      conjugated_relaxation2
    use comp_parameters_mod
    use positions_mod
    use cgm_storage_mod

    integer i_passages,i_atoms,i_directions,grad_max_index,cgm_direction_multiplicator
    real(8) gamma_curr,gamma_prev,gamma_sum
    real(8) scalar_curr,scalar_prev,scalar_summ
    real(8) e_curr,e_prev,D_denominator,shift_magnitude
    real(8) fx,fy,fz,grad_curr,maxgrad
    logical need_to_rescale_total_gradient,relax_the_max
!real(8), dimension(1:3,atoms_max_array) ::  F_prev,F_curr,D_prev,D_curr gamma_

    print*,"there will be ",total_relax_passages,&
    "conjugated gradient descent relaxation stages, this is 1/3 cgr."
    relax_the_max=.true.
    relax_the_max=.false.
    maxgrad=0d0
    grad_curr=0d0
    scalar_summ =   0d0
    scalar_prev =   1d0
    grad_max_index=0
        !@ obtain current force
        do  i_atoms=first_relaxable,last__relaxable
            if(abs(R_curr(1,i_atoms)).gt.26d-1*286d-2)cycle
            call energygradient_pw(i_atoms,fx,fy,fz)!gradient of a single atom
            F_curr(1,i_atoms)=fx
            F_curr(2,i_atoms)=fy
            F_curr(3,i_atoms)=fz
            grad_curr=(fx*fx+fy*fy+fz*fz)
            scalar_summ=scalar_summ+grad_curr
            if(maxgrad.lt.grad_curr)then
                maxgrad = grad_curr
                grad_max_index=i_atoms
            endif
        enddo
        print 212,maxgrad , " is maximum grad,atom# ",grad_max_index," "
        scalar_curr=scalar_summ
        !@ current force is obtained
        !@ check current force if it is small enough to finalize procedure or rescale
        if  (scalar_curr.lt.system_gradient_small_enough**2) then
            print*," total gradient have been weighed ... and found wanting. more precisely, less than "&
            ,system_gradient_small_enough
            return
        endif
        need_to_rescale_total_gradient=.false.
        NEED_TO_RESCALE_TOTAL_GRADIENT=.TRUE.
        D_DENOMINATOR=5D-3
        !D_denominator=(5d-1/gamma_curr)
        if  (scalar_curr.gt.sqrt(1d-1*dfloat(last__relaxable-first_relaxable+1))) then
            print*," total gradient have been weighed ... and found big (see it further). rescaling ensues."&
            ,scalar_curr
            need_to_rescale_total_gradient=.true.
            D_denominator=(7d-1/scalar_curr)
            D_DENOMINATOR=5D-3
        endif
        !gamma_curr=scalar_curr/scalar_prev

        D_curr(1:3,first_relaxable:last__relaxable)=&
        F_curr(1:3,first_relaxable:last__relaxable)!+&
        !D_prev(1:3,first_relaxable:last__relaxable)*&
        !gamma_curr !not needed at that step
        !shift_magnitude=norm2(D_curr(1:3,first_relaxable:last__relaxable))
        shift_magnitude=sqrt(abs(scalar_prev*D_denominator))
        if(need_to_rescale_total_gradient)then
            D_curr(1:3,first_relaxable:last__relaxable)=&
            D_curr(1:3,first_relaxable:last__relaxable)*D_denominator
            shift_magnitude=sqrt(scalar_prev*D_denominator)
        endif
        scalar_prev=scalar_curr
        D_prev(1:3,first_relaxable:last__relaxable)=D_curr(1:3,first_relaxable:last__relaxable)
        !F_prev(1:3,first_relaxable:last__relaxable)=F_curr(1:3,first_relaxable:last__relaxable) !do i really need prev F ?
        !minimize along d_curr
        call    energy_of_system(e_curr)
        cgm_direction_multiplicator=1
        if(e_curr.gt.0d0)cgm_direction_multiplicator=3
        do i_directions=1,cgm_direction_switch*cgm_direction_multiplicator
            e_prev  =e_curr
            R_curr(1:3,first_relaxable:last__relaxable)=&
            R_curr(1:3,first_relaxable:last__relaxable)+&
            D_curr(1:3,first_relaxable:last__relaxable)
            call    set_positions_of_periodic_atoms
            call    energy_of_system(e_curr)

            if((e_curr-e_prev).gt.abs(1d-4*e_prev)) then
            !this is dramatic increase of energy. something went very wrong
                call    wo_xyz_snapshot!@   write the image of trajectory
                R_curr(1:3,first_relaxable:last__relaxable)=&
                R_curr(1:3,first_relaxable:last__relaxable)-&
                D_curr(1:3,first_relaxable:last__relaxable)
            endif

            if(e_curr.gt.e_prev) then
                D_curr(1:3,first_relaxable:last__relaxable)=&
                D_curr(1:3,first_relaxable:last__relaxable)*(-cgm_denominator)
                !print"(A,$)","|reversing direction and dividing it by two | "
                shift_magnitude=shift_magnitude*abs(cgm_denominator)
                if(shift_magnitude.lt.system_gradient_small_enough)then !maybe another criteria
                    call    wo_xyz_snapshot!@   write the image of trajectory
                    print*,"Current energy is ",e_curr
                    exit
                endif
            endif
            call    wo_xyz_snapshot!@   write the image of trajectory
        enddo!@ minimization along D_curr is done


    if(relax_the_max)call    relax_atom(grad_max_index)

    main_cgd_passages_loop : &
    do i_passages=2,total_relax_passages !@ now do it in a cycle
        maxgrad=0d0
        grad_curr=0d0

        scalar_summ =   0d0
        !@ obtain current force
        do  i_atoms=first_relaxable,last__relaxable
            if(abs(R_curr(1,i_atoms)).gt.26d-1*286d-2)cycle
            call energygradient_pw(i_atoms,fx,fy,fz)!gradient of a single atom
            F_curr(1,i_atoms)=fx
            F_curr(2,i_atoms)=fy
            F_curr(3,i_atoms)=fz
            grad_curr=(fx*fx+fy*fy+fz*fz)
            scalar_summ=scalar_summ+grad_curr
            if(maxgrad.lt.grad_curr) then
                maxgrad = grad_curr
                grad_max_index=i_atoms
            endif
        enddo
        print 212,maxgrad , " is maximum grad,atom# ",grad_max_index," "
        scalar_curr=scalar_summ
        !@ current force is obtained
        !@ check current force if it is small enough to finalize procedure or rescale
        if  (scalar_curr.lt.system_gradient_small_enough**2) then
            print*," total gradient have been weighed ... and found wanting. more precisely, less than "&
            ,system_gradient_small_enough
            exit main_cgd_passages_loop
        endif
        need_to_rescale_total_gradient=.false.
        NEED_TO_RESCALE_TOTAL_GRADIENT=.TRUE.
        D_DENOMINATOR=5D-3
        if  (scalar_curr.gt.sqrt(dfloat(last__relaxable-first_relaxable+1))) then
            print*," total gradient have been weighed ... and found big (see it further). rescaling ensues."&
            ,scalar_curr
            need_to_rescale_total_gradient=.true.
            D_denominator=sqrt(2d0/scalar_curr)
            !D_DENOMINATOR=5D-3
        endif
        gamma_curr=scalar_curr/scalar_prev
        !@ check is done, now find the direction of descending
        D_curr(1:3,first_relaxable:last__relaxable)=&
        F_curr(1:3,first_relaxable:last__relaxable)+&
        D_prev(1:3,first_relaxable:last__relaxable)*&
        gamma_curr

        shift_magnitude=sqrt(scalar_prev)

        if(need_to_rescale_total_gradient)then
            D_curr(1:3,first_relaxable:last__relaxable)=&
            D_curr(1:3,first_relaxable:last__relaxable)*D_denominator
            shift_magnitude=sqrt(scalar_prev*D_denominator)
        endif
        scalar_prev=scalar_curr

        D_prev(1:3,first_relaxable:last__relaxable)=D_curr(1:3,first_relaxable:last__relaxable) !@ save the directions
        !@ minimize along D_curr
        call    energy_of_system(e_curr)
        cgm_direction_multiplicator=1
        if(e_curr.gt.0d0)cgm_direction_multiplicator=3
        do i_directions=1,cgm_direction_switch*cgm_direction_multiplicator+(i_passages/2)


            e_prev  =e_curr
            R_curr(1:3,first_relaxable:last__relaxable)=&
            R_curr(1:3,first_relaxable:last__relaxable)+&
            D_curr(1:3,first_relaxable:last__relaxable)
            call    set_positions_of_periodic_atoms
            call    energy_of_system(e_curr)
            !PRINT*, I_PASSAGES, I_DIRECTIONS

            if((e_curr-e_prev).gt.abs(8d-5*e_prev)) then
            !this is dramatic increase of energy. something went very wrong
                call    wo_xyz_snapshot!@   write the image of trajectory
                R_curr(1:3,first_relaxable:last__relaxable)=&
                R_curr(1:3,first_relaxable:last__relaxable)-&
                D_curr(1:3,first_relaxable:last__relaxable)
            endif

            if(e_curr.gt.e_prev) then
                D_curr(1:3,first_relaxable:last__relaxable)=&
                D_curr(1:3,first_relaxable:last__relaxable)*(-cgm_denominator)
                !print"(A,$)","|reversing direction and dividing it by two | "
                shift_magnitude=shift_magnitude*abs(cgm_denominator)
                if(shift_magnitude.lt.system_gradient_small_enough)then !maybe another criteria
                    call    wo_xyz_snapshot!@   write the image of trajectory
                    exit
                endif
            endif
            call    wo_xyz_snapshot!@   write the image of trajectory
        enddo!@ minimization along D_curr is done
        print*,"Current energy is ",e_curr
        if(relax_the_max)call    relax_atom(grad_max_index)
    enddo   main_cgd_passages_loop
212 format(1x,E11.4,A,I6.2,A,$)
    print*," "
    print*,"we are done with conjugated descending2"
endsubroutine   conjugated_relaxation2

subroutine      cast_gusev_shifts_anisotropic
    use phys_parameters_mod
    use positions_mod
    use comp_parameters_mod
    integer i,i_beauty
    real(8),parameter   :: okr=1d-3           !small epsilon
    real(8) x,y,z, ux,uy,uz,phi_1,bv,nu
    if(core_sign.gt.0)then
        print*, "casting anisotropic Gusev shifts:     __compressed core ---"
    else
        print*, "casting anisotropic Gusev shifts:     decompressed core +++"
    endif
    !ux=0d0;uy=0d0;uz=0d0
    do i_beauty=1,beauty_denominator
        do i=first______wall,last_______wall
            ux=0d0;uy=0d0;uz=0d0
            x=R_perf(1,i); y=R_perf(2,i); z=R_perf(3,i); phi_1=atan2(y,x)

            call    interpolated_burgers_poisson(phi_1,bv,nu)

            if(abs(x).gt.okr) then
                if((x .gt.0d0))then
                    ux= bv/(pi*2d0)*(&
                    (atan(y/x) + pi*( 05d-1) + x*y*5d-1/(1-nu)/(x*x+y*y)) )
                else
                    ux=bv/(pi*2d0)*& !
                    (atan(y/x) + pi*(-05d-1) + x*y*5d-1/(1-nu)/(x*x+y*y))
                endif
                uy=-bv/(pi*4d0)*(&
                (1d0-2d0*nu)/(2d0-2d0*nu)*log(x*x+y*y)-&
                y*y/((1-nu)*(x*x+y*y) ) )
            else!co-extraplane intermediate
                if(y.gt.0d0)then
                    uy=-a0*5d-1*core_sign
                    uz= a0*5d-1
                    uy=uy+  (&
                        -bv/(pi*4d0)*(&
                        (1d0-2d0*nu)/(2d0-2d0*nu)*log(x*x+y*y)-&
                        y*y/((1-nu)*(x*x+y*y) ) )&
                            )
                else
                    if(y.lt.-okr) then
                        uy=-bv/(pi*4d0)*(&
                        (1d0-2d0*nu)/(2d0-2d0*nu)*log(x*x+y*y)-&
                        y*y/((1-nu)*(x*x+y*y) ) )
                    endif
                endif

            endif
            R_curr(1,i)=R_perf(1,i)+ux*float(i_beauty-1)/float(beauty_denominator)
            R_curr(2,i)=R_perf(2,i)+uy*float(i_beauty-1)/float(beauty_denominator)
            R_curr(3,i)=R_perf(3,i)+uz*float(i_beauty-1)/float(beauty_denominator)
        enddo
            call    wo_xyz_snapshot!@   write the image of trajectory
    enddo !beauty

endsubroutine   cast_gusev_shifts_anisotropic

subroutine      interpolated_burgers_poisson(angle,bv,nu)
    use anisotropy_mod
    use comp_parameters_mod
    use phys_parameters_mod
    real(8),intent(in )  :: angle
    real(8),intent(out)  :: bv,nu
    real(8) phi_a,phi_b,burg_a,burg_b,pois_a,pois_b
    integer i,j,na,nb

    BV=BURGERS; NU=POISSON; RETURN

    na=max_layer_list_index;nb=1 !closure for the cirle at +- pi
    phi_a   =phi_____list(na)
    phi_b   =phi_____list(nb)+2*pi
!    BURG_A  =-123D-4
!    BURG_B  =-123D-4
!    POIS_A  =-123D-4
!    POIS_B  =-123D-4
    do j=2,max_layer_list_index
!        PRINT*,"#",PHI_____LIST(J)," RADIANS   # ",j,nb
!        PRINT*,"#",PHI_____LIST(J+1)," RADIANS   # ",j+1,na
!        PRINT*,"#",burgers_list(J)," BURGERSES # ",layer_list(j-1)
!        PRINT*,"#",poisson_list(J)," POISSONS  # ",layer_list(j)
        if(phi_____list(j).lt. angle)cycle
        na=j-1!layer_list(j-1)
        nb=j!layer_list(j  )
        phi_a   =phi_____list(na)
        phi_b   =phi_____list(nb)
        !PRINT'(1x,I5.1,1x,I6.2,1x,$)',J,NB
        if(abs(phi_a-phi_b).lt.1d-6)cycle
!            PRINT*,"SMALL1",NA,NB,PHI_A,PHI_B,PHI_____LIST(NA),PHI_____LIST(NB),ANGLE
        exit
    enddo
    !if(abs(phi_a-phi_b).lt.1d-6)then
        !PRINT*,"SMALL2",NA,NB,PHI_A,PHI_B,PHI_____LIST(NA),PHI_____LIST(NB),ANGLE
    !endif
    burg_a  =burgers_list(na)
    burg_b  =burgers_list(nb)
    pois_a  =poisson_list(na)
    pois_b  =poisson_list(nb)
    !PRINT*,BURG_A,BURG_B,POIS_A,POIS_B,"AAAAAA"
    !STOP"ENOUGH"
    bv=burg_a+((burg_b-burg_a)/(phi_b-phi_a))*(angle-phi_a)
    nu=pois_a+((pois_b-pois_a)/(phi_b-phi_a))*(angle-phi_a)
endsubroutine   interpolated_burgers_poisson
