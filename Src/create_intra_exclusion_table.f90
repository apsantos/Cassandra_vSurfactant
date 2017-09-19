!********************************************************************************
!   Cassandra - An open source atomistic Monte Carlo software package
!   developed at the University of Notre Dame.
!   http://cassandra.nd.edu
!   Prof. Edward Maginn <ed@nd.edu>
!   Copyright (2013) University of Notre Dame du Lac
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*******************************************************************************

SUBROUTINE Create_Intra_Exclusion_Table
  
  !-----------------------------------------------------------------------------
  ! This routine takes the input vdw and charge scaling rules and computes
  ! a arrays of dimensions (natoms,natoms,nspecies) that then get multiplied
  ! by the potential interactions for those atom pairs when performing 
  ! intramolecular potential energy calculations.
  ! Net result is creation of two arrays:
  ! vdw_intra_scale(ii,jj,is)
  ! charge_intra_scale(ii,jj,is)
  ! where ii and jj refer to atom number is species is. To compute 
  ! intramolecular NB energy, we loop over all pairs and multiply 
  ! by these scaling factors. 
  !
  ! Called by
  !
  !  gcmc_control
  !  gemc_control
  !  nptmc_control
  !  nvtmc_control
  !  nvt_mc_fragment_control
  !
  ! Revision history
  !
  !   12/10/13 : Beta Release
  !-----------------------------------------------------------------------------  
  USE Run_Variables
  USE Type_Definitions
  USE File_Names

  IMPLICIT NONE

  INTEGER :: is,ii,jj,kk, max_natoms, itype, jtype
!-----------------------------------------------------------------------------

  max_natoms = MAXVAL(natoms)
  ALLOCATE(vdw_intra_scale(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)
  ALLOCATE(charge_intra_scale(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)
  vdw_intra_scale(:,:,:) = 0.0_DP
  charge_intra_scale(:,:,:) = 0.0_DP

  ALLOCATE(vdw_in_param1_table(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)
  ALLOCATE(vdw_in_param2_table(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)
  ALLOCATE(vdw_in_param3_table(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)
  ALLOCATE(vdw_in_param4_table(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)
  ALLOCATE(vdw_in_param5_table(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)
  ALLOCATE(vdw_in_param6_table(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)
  ALLOCATE(vdw_in_param7_table(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)
  ALLOCATE(vdw_in_param8_table(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)
  ALLOCATE(vdw_in_param9_table(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)
  ALLOCATE(vdw_in_param10_table(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)
  ALLOCATE(vdw_in_param11_table(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)
  ALLOCATE(vdw_in_param12_table(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)
  ALLOCATE(vdw_in_param13_table(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)

  vdw_in_param1_table = 0.0_DP
  vdw_in_param2_table = 0.0_DP
  vdw_in_param3_table = 0.0_DP
  vdw_in_param4_table = 0.0_DP
  vdw_in_param5_table = 0.0_DP
  vdw_in_param6_table = 0.0_DP
  vdw_in_param7_table = 0.0_DP
  vdw_in_param8_table = 0.0_DP
  vdw_in_param9_table = 0.0_DP
  vdw_in_param10_table = 0.0_DP
  vdw_in_param11_table = 0.0_DP
  vdw_in_param12_table = 0.0_DP
  vdw_in_param13_table = 0.0_DP

  ALLOCATE(rcut_in_vdw_mix(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)
  rcut_in_vdw_mix(:,:,:) = rcut_vdw(1)
  ALLOCATE(rcut_in_vdwsq_mix(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)
  rcut_in_vdwsq_mix(:,:,:) = 0.0_DP
  ALLOCATE(int_in_vdw_style_mix(max_natoms,max_natoms,nspecies,15), Stat=AllocateStatus)
  int_in_vdw_style_mix(:,:,:,:) = .false.
  int_in_vdw_style_mix(:,:,:,int_vdw_style(1)) = .true.
  ALLOCATE(int_in_vdw_sum_style_mix(max_natoms,max_natoms,nspecies), Stat=AllocateStatus)
  int_in_vdw_sum_style_mix(:,:,:) = int_vdw_sum_style(1)

  IF (AllocateStatus .NE. 0) THEN
     err_msg = ''
     err_msg(1) = ' ERROR: Not enough memory for scaling tables '
     CALL Clean_Abort(err_msg,'ceate_intra_exclusion_table')
  END IF
  
  DO is=1,nspecies

     DO ii=1,natoms(is)

        DO jj = 1,natoms(is)

           ! Self interactions are automatically excluded. 
           IF (ii == jj) THEN
              vdw_intra_scale(ii,jj,is) = 0.0_DP
              charge_intra_scale(ii,jj,is) = 0.0_DP
           ELSE
              !Set all other interactions to default "1_N" scaling (normally 1.0) 
              ! and turn selected values off using otherinput scaling values
              vdw_intra_scale(ii,jj,is) = scale_1_N_vdw(is)
              charge_intra_scale(ii,jj,is) = scale_1_N_charge(is)
           ENDIF

        ENDDO

     ENDDO

  ENDDO

  ! Now apply specific 1-2, 1-3 and 1-4 scaling rules. 

  SpeciesLoop: DO is = 1,nspecies

     ! 1-4 scaling via dihedrals. Do this first. Later, if atoms
     ! connected by a dihedral and by an angle (for example) are
     ! encountered, then the priority scaling overwrites (i.e. 1-2
     ! has priority over 1-3, which has priority over 1-4).
     ! For example, in the compound below 1-2-3-4 is a 1-4 interaction 
     ! between 1 and 4,but via 1-5-4 it is a 1-3 interaction. We thus
     ! use the 1-3 interaction value. 
     ! 
     !                       1
     !                   /       \           7   10
     !                  /         \          |    |
     !                 2           5----6----9---12
     !                  \         /          |    |
     !                   3------ 4           8   11
     !

     DO kk=1,ndihedrals(is)
        ii = dihedral_list(kk,is)%atom1
        jj = dihedral_list(kk,is)%atom4

        vdw_intra_scale(ii,jj,is) = scale_1_4_vdw(is)
        vdw_intra_scale(jj,ii,is) = scale_1_4_vdw(is)
        charge_intra_scale(ii,jj,is) = scale_1_4_charge(is)
        charge_intra_scale(jj,ii,is) = scale_1_4_charge(is)
        
     ENDDO


     ! 1-3 scaling via angles
     DO kk = 1,nangles(is)

        ii = angle_list(kk,is)%atom1
        jj = angle_list(kk,is)%atom3

        vdw_intra_scale(ii,jj,is) = scale_1_3_vdw(is)
        vdw_intra_scale(jj,ii,is) = scale_1_3_vdw(is)
        charge_intra_scale(ii,jj,is) = scale_1_3_charge(is)
        charge_intra_scale(jj,ii,is) = scale_1_3_charge(is)

     ENDDO

     ! 1-2 scaling via bonds     
     DO kk=1,nbonds(is)

        ii = bond_list(kk,is)%atom1
        jj = bond_list(kk,is)%atom2

        vdw_intra_scale(ii,jj,is) = scale_1_2_vdw(is)
        vdw_intra_scale(jj,ii,is) = scale_1_2_vdw(is)
        charge_intra_scale(ii,jj,is) = scale_1_2_charge(is)
        charge_intra_scale(jj,ii,is) = scale_1_2_charge(is)

     ENDDO
     
  ENDDO SpeciesLoop

  DO is=1,nspecies
     ! set all values based on scaling
     DO ii=1,natoms(is)
        itype = nonbond_list(ii,is)%atom_type_number
        DO jj = 1,natoms(is)
           jtype = nonbond_list(jj,is)%atom_type_number

           vdw_in_param1_table(ii,jj,is) = vdw_param1_table(itype,jtype) * vdw_intra_scale(ii,jj,is)
           vdw_in_param2_table(ii,jj,is) = vdw_param2_table(itype,jtype) * vdw_intra_scale(ii,jj,is)
           vdw_in_param3_table(ii,jj,is) = vdw_param3_table(itype,jtype) * vdw_intra_scale(ii,jj,is)
           vdw_in_param4_table(ii,jj,is) = vdw_param4_table(itype,jtype) * vdw_intra_scale(ii,jj,is)
           vdw_in_param5_table(ii,jj,is) = vdw_param5_table(itype,jtype) * vdw_intra_scale(ii,jj,is)
           vdw_in_param6_table(ii,jj,is) = vdw_param6_table(itype,jtype) * vdw_intra_scale(ii,jj,is)
           vdw_in_param7_table(ii,jj,is) = vdw_param7_table(itype,jtype) * vdw_intra_scale(ii,jj,is)
           vdw_in_param8_table(ii,jj,is) = vdw_param8_table(itype,jtype) * vdw_intra_scale(ii,jj,is)
           vdw_in_param9_table(ii,jj,is) = vdw_param9_table(itype,jtype) * vdw_intra_scale(ii,jj,is)
           vdw_in_param10_table(ii,jj,is) = vdw_param10_table(itype,jtype) * vdw_intra_scale(ii,jj,is)
           vdw_in_param11_table(ii,jj,is) = vdw_param11_table(itype,jtype) * vdw_intra_scale(ii,jj,is)
           vdw_in_param12_table(ii,jj,is) = vdw_param12_table(itype,jtype) * vdw_intra_scale(ii,jj,is)
           vdw_in_param13_table(ii,jj,is) = vdw_param13_table(itype,jtype) * vdw_intra_scale(ii,jj,is)

        ENDDO
     ENDDO

     ! Overwrite default values if defined by the user
     IF (intrafile_name(is) .NE. "") THEN
        CALL Read_Intra_Exclusion_Table(is)
     ENDIF
  ENDDO
     

  ! report info to log
  WRITE(logunit,*)
  WRITE(logunit,*) '*** Creating exclusion table ***'
  WRITE(logunit,*)    

  WRITE(logunit,'(5x,a)') 'species   atom1      atom2  vdw scale q-q scale'

  DO is=1,nspecies
     DO ii=1,natoms(is)
        DO jj = 1,natoms(is)
           WRITE(logunit,'(3(I10),2(F10.3))') is,ii,jj,vdw_intra_scale(ii,jj,is),&
                charge_intra_scale(ii,jj,is)
        ENDDO
     ENDDO
  ENDDO
  

  WRITE(logunit,*)
  WRITE(logunit,*) '*** Completed construction of exclusion table ***'
  WRITE(logunit,*)    


END SUBROUTINE Create_Intra_Exclusion_Table

SUBROUTINE Read_Intra_Exclusion_Table(is)
  !-----------------------------------------------------------------------------
  ! Reads in a file that explicitly defines the intra scaling of user-defined
  ! atoms.
  ! If they are not defined by the user then, the default values are used.
  ! Throws errors if the names in the file do not coincide with MCF files or if
  ! the two atoms are not part of the same molecule.
  !
  ! Called by
  !
  !   Create_Intra_Exclusion_Table
  !
  ! Revision history
  !
  !   8/20/15 
  !
  !   - Andrew P. Santos
  !-----------------------------------------------------------------------------  
  USE Run_Variables
  USE Type_Definitions
  USE IO_Utilities
  USE File_Names

  IMPLICIT NONE

  INTEGER, INTENT(in) :: is
  INTEGER :: ia,ii,jj,js,ja, i
  INTEGER :: ierr,nbr_entries, i_line
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: temp_type

  CHARACTER(charLength) :: line_string, line_array(lineArrayLength)
  CHARACTER(charLength) :: temp_name, pot_type
!-----------------------------------------------------------------------------
  ! Open intra scalingfile and find the line where the data begins
  OPEN(UNIT=intrafile_unit,FILE=intrafile_name(is),STATUS="OLD",IOSTAT=openstatus,ACTION="READ")

  CALL Read_String(intrafile_unit,line_string,ierr)
  IF (ierr .NE. 0) THEN
      err_msg = ""
      err_msg(1) = "Error reading intra scaling values."
      err_msg(2) = intrafile_name(is)
      CALL Clean_Abort(err_msg,'Read_Intra_Exclusion_Table')
  END IF
  i_line = 1
  ALLOCATE(temp_type( SUM(natoms(:)), nspecies ) )
  temp_type = 0
  DO 
      IF (ierr < 0) EXIT

      IF (line_string(1:12) == '# Atom_Types') THEN
         DO ia = 1, natoms(is)
            CALL Parse_String(intrafile_unit, i_line, 1, nbr_entries, line_array, ierr)
            temp_name = line_array(1)
            IF (ANY(temp_name == nonbond_list(:,is)%atom_name ) .eqv. .FALSE.) THEN
               err_msg(1) = "Atom name and type in Intra Scaling table (" &
                             //TRIM(temp_name)//") does not match the MCF file."
               err_msg(2) = "Make sure the atom list in the file starts with "// &
                            "species 1 then 2 and so on."
               CALL Clean_Abort(err_msg,'Read_Intra_Exclusion_Table')
            ENDIF

            temp_type(ia, is) = String_To_Int(line_array(2)) 

            i_line = i_line + 1
         ENDDO
         DO js = is+1, nspecies
            IF (ANY(temp_type(:natoms(is),is) == temp_type(:natoms(js),js))) THEN
               err_msg(1) = "Atom types between among species"
               err_msg(2) = "must be unique in the intra scaling table."
               CALL Clean_Abort(err_msg,'Read_Intra_Exclusion_Table')
            ENDIF
         ENDDO

      ELSEIF (line_string(1:15) == '# Intra_Scaling') THEN
         ! read in values
         DO
            CALL Parse_String(intrafile_unit, i_line, 2, nbr_entries, line_array, ierr)
            IF (TRIM(line_array(2)) == 'Done_Intra_Scaling') THEN
              EXIT
            ENDIF
 
            ii = String_To_Int(line_array(1))
            jj = String_To_Int(line_array(2))

            IF (( ANY(ii == temp_type(:,is)) .eqv. .FALSE.) .or. &
                ( ANY(jj == temp_type(:,is)) .eqv. .FALSE.)) THEN
               err_msg(1) = "Intra Scaling can only be done among atoms in the same molecule."
               CALL Clean_Abort(err_msg,'Read_Intra_Exclusion_Table')
               EXIT
            ENDIF

            DO ia = 1, natoms(is)
               IF ( ii == temp_type(ia,is)) THEN
                  DO ja = 1, natoms(is)
                     IF ( jj == temp_type(ja,is)) THEN
                        IF (int_vdw_style(1) /= vdw_none) THEN
                           vdw_intra_scale(ia,ja,is) = String_To_Double(line_array(3))
                           vdw_intra_scale(ja,ia,is) = String_To_Double(line_array(3))
                        ENDIF
                        IF (int_charge_style(1) /= charge_none) THEN
                           IF (nbr_entries == 4) THEN
                              charge_intra_scale(ia,ja,is) = String_To_Double(line_array(4))
                              charge_intra_scale(ja,ia,is) = String_To_Double(line_array(4))
                           ELSE
                              charge_intra_scale(ia,ja,is) = String_To_Double(line_array(3))
                              charge_intra_scale(ja,ia,is) = String_To_Double(line_array(3))
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            i_line = i_line + 1
         ENDDO

      ELSEIF (line_string(1:15) == '# Mixing_Values') THEN
         ! read in values
         DO
            CALL Parse_String(intrafile_unit, i_line, 2, nbr_entries, line_array, ierr)
            IF (TRIM(line_array(2)) == 'Done_Mixing_Values') THEN
              EXIT
            ENDIF
 
            ii = String_To_Int(line_array(1))
            jj = String_To_Int(line_array(2))

            IF (( ANY(ii == temp_type(:,is)) .eqv. .FALSE.) .or. &
                ( ANY(jj == temp_type(:,is)) .eqv. .FALSE.)) THEN
               err_msg(1) = "Intra Scaling can only be done among atoms in the same molecule."
               CALL Clean_Abort(err_msg,'Read_Intra_Exclusion_Table')
               EXIT
            ENDIF

            DO ia = 1, natoms(is)
               IF ( ii == temp_type(ia,is)) THEN
                  DO ja = 1, natoms(is)
                     IF ( jj == temp_type(ja,is)) THEN
                        int_in_vdw_style_mix(ia,ja,is,:) = .false.
                        int_in_vdw_style_mix(ja,ia,is,:) = .false.
                
                        DO i = 3, SIZE(line_array)
                            pot_type = line_array(i)
                
                            IF (pot_type == 'LJ' .or. pot_type == 'LJ126' .or. &
                                pot_type == 'LJ124' .or. &
                                pot_type == 'LJ96') THEN

                                IF (pot_type == 'LJ' .or. pot_type == 'LJ126') THEN
                                   int_in_vdw_style_mix(ia,ja,is,vdw_lj) = .true.
                                   int_in_vdw_style_mix(ja,ia,is,vdw_lj) = .true.
                                ELSE IF (pot_type == 'LJ124') THEN
                                   int_in_vdw_style_mix(ia,ja,is,vdw_lj124) = .true.
                                   int_in_vdw_style_mix(ja,ia,is,vdw_lj124) = .true.
                                ELSE IF (pot_type == 'LJ96') THEN
                                   int_in_vdw_style_mix(ia,ja,is,vdw_lj96) = .true.
                                   int_in_vdw_style_mix(ja,ia,is,vdw_lj96) = .true.
                                END IF

                                vdw_in_param1_table(ia,ja,is) = String_To_Double(line_array(i+1))
                                vdw_in_param2_table(ia,ja,is) = String_To_Double(line_array(i+2))
                                vdw_in_param1_table(ja,ia,is) = vdw_in_param1_table(ia,ja,is)
                                vdw_in_param2_table(ja,ia,is) = vdw_in_param2_table(ia,ja,is)
                
                            ELSEIF (pot_type == 'HYDR') THEN
                                int_in_vdw_style_mix(ia,ja,is,vdw_hydra) = .true.
                                int_in_vdw_style_mix(ja,ia,is,vdw_hydra) = .true.
                                vdw_in_param3_table(ia,ja,is) = String_To_Double(line_array(i+1))
                                vdw_in_param4_table(ia,ja,is) = String_To_Double(line_array(i+2))
                                vdw_in_param5_table(ia,ja,is) = String_To_Double(line_array(i+3))
                                vdw_in_param3_table(ja,ia,is) = vdw_in_param3_table(ia,ja,is)
                                vdw_in_param4_table(ja,ia,is) = vdw_in_param4_table(ia,ja,is)
                                vdw_in_param5_table(ja,ia,is) = vdw_in_param5_table(ia,ja,is)
                
                            ELSEIF (pot_type == 'CORR') THEN
                                int_in_vdw_style_mix(ia,ja,is,vdw_corr) = .true.
                                int_in_vdw_style_mix(ja,ia,is,vdw_corr) = .true.
                                vdw_in_param6_table(ia,ja,is) = String_To_Double(line_array(i+1))
                                vdw_in_param6_table(ja,ia,is) = vdw_in_param6_table(ia,ja,is)
                                vdw_in_param7_table(ia,ja,is) = String_To_Double(line_array(i+1))
                                vdw_in_param7_table(ja,ia,is) = vdw_in_param7_table(ia,ja,is)
                
                            ELSEIF (pot_type == 'Yukawa') THEN
                                int_in_vdw_style_mix(ia,ja,is,vdw_yukawa) = .true.
                                int_in_vdw_style_mix(ja,ia,is,vdw_yukawa) = .true.
                                vdw_in_param8_table(ia,ja,is) = String_To_Double(line_array(i+1))
                                vdw_in_param8_table(ja,ia,is) = vdw_in_param8_table(ia,ja,is)
                                vdw_in_param9_table(ia,ja,is) = String_To_Double(line_array(i+2))
                                vdw_in_param9_table(ja,ia,is) = vdw_in_param9_table(ia,ja,is)
                
                            ELSEIF (pot_type == 'SCR') THEN
                                int_in_vdw_style_mix(ia,ja,is,vdw_screen) = .true.
                                int_in_vdw_style_mix(ja,ia,is,vdw_screen) = .true.
                                vdw_in_param12_table(ia,ja,is) = &
                                    SQRT( (String_To_Double(line_array(i+1)) * navogadro / m3_to_A3 ) * &
                                          (charge_factor(1) * 8.0 * PI) / beta(1) )
                                vdw_in_param13_table(ia,ja,is) = String_To_Double(line_array(i+2)) 
                                vdw_in_param12_table(ja,ia,is) = vdw_in_param12_table(ia,ja,is)
                                vdw_in_param13_table(ja,ia,is) = vdw_in_param13_table(ia,ja,is)

                            ELSEIF (pot_type == 'SW') THEN
                                int_in_vdw_style_mix(ia,ja,is,vdw_sw) = .true.
                                int_in_vdw_style_mix(ja,ia,is,vdw_sw) = .true.
                                vdw_in_param10_table(ia,ja,is) = String_To_Double(line_array(i+1))
                                vdw_in_param10_table(ja,ia,is) = vdw_in_param10_table(ia,ja,is)
                                vdw_in_param11_table(ia,ja,is) = String_To_Double(line_array(i+2))
                                vdw_in_param11_table(ja,ia,is) = vdw_in_param11_table(ia,ja,is)
                
                            ENDIF
                
                        ENDDO
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            i_line = i_line + 1
         ENDDO

      ELSEIF (line_string(1:11) == '# VDW_Style') THEN
         ! read in values
         DO
            CALL Parse_String(intrafile_unit, i_line, 2, nbr_entries, line_array, ierr)
            IF (TRIM(line_array(2)) == 'Done_VDW_Style') THEN
              EXIT
            ENDIF
 
            ii = String_To_Int(line_array(1))
            jj = String_To_Int(line_array(2))

            IF (( ANY(ii == temp_type(:,is)) .eqv. .FALSE.) .or. &
                ( ANY(jj == temp_type(:,is)) .eqv. .FALSE.)) THEN
               err_msg(1) = "Intra Scaling can only be done among atoms in the same molecule."
               CALL Clean_Abort(err_msg,'Read_Intra_Exclusion_Table')
               EXIT
            ENDIF

            DO ia = 1, natoms(is)
               IF ( ii == temp_type(ia,is)) THEN
                  DO ja = 1, natoms(is)
                     IF ( jj == temp_type(ja,is)) THEN
                        rcut_in_vdw_mix(ia,ja,is) = String_To_Double(line_array(4))
                        rcut_in_vdw_mix(ja,ia,is) = String_To_Double(line_array(4))

                        IF (line_array(3) == 'cut') THEN
                           int_in_vdw_sum_style_mix(ia,ja,is) = vdw_cut
                        ELSE IF (line_array(3) == 'cut_shift') THEN
                           int_in_vdw_sum_style_mix(ia,ja,is) = vdw_cut_shift
                        ELSE IF (line_array(3) == 'cut_switch') THEN
                           int_in_vdw_sum_style_mix(ia,ja,is) = vdw_cut_switch
                        ELSE IF (line_array(3) == 'CHARMM') THEN
                           int_in_vdw_sum_style_mix(ia,ja,is) = vdw_charmm
                        ENDIF
                        int_in_vdw_sum_style_mix(ja,ia,is) = int_in_vdw_sum_style_mix(ia,ja,is)
                
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            i_line = i_line + 1
         ENDDO

      ENDIF

      i_line = i_line + 1
      CALL Read_String(intrafile_unit,line_string,ierr)
  ENDDO
  CLOSE(UNIT=intrafile_unit)

END SUBROUTINE Read_Intra_Exclusion_Table
