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

  INTEGER :: is,ii,jj,kk
!-----------------------------------------------------------------------------

  ALLOCATE(vdw_intra_scale(MAXVAL(natoms),MAXVAL(natoms),nspecies), Stat=AllocateStatus)
  ALLOCATE(charge_intra_scale(MAXVAL(natoms),MAXVAL(natoms),nspecies), Stat=AllocateStatus)

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

  ! Overwrite values if defined by the user
  IF (intrafile_name .NE. "") THEN
     CALL Read_Intra_Exclusion_Table
  ENDIF

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

SUBROUTINE Read_Intra_Exclusion_Table
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

  INTEGER :: is, ia, ja, ii, jj
  INTEGER :: ierr,nbr_entries, i_line
  INTEGER :: t_atoms
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: temp_type

  CHARACTER(240) :: line_string, line_array(80)
  CHARACTER(240) :: temp_name
!-----------------------------------------------------------------------------
    ! Open intra scalingfile and find the line where the data begins
    OPEN(UNIT=intrafile_unit,FILE=intrafile_name,STATUS="OLD",IOSTAT=openstatus,ACTION="READ")
    CALL Read_String(intrafile_unit,line_string,ierr)
    IF (ierr .NE. 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading intra scaling values."
        CALL Clean_Abort(err_msg,'Read_Intra_Exclusion_Table')
    END IF
    i_line = 1
    ALLOCATE(temp_type( SUM(natoms(:)), nspecies ) )
    DO 
        IF (line_string(1:12) == '# Atom_Types') THEN
            DO is = 1, nspecies
               DO ia = 1, natoms(is)
                  CALL Parse_String(intrafile_unit, i_line, 1, nbr_entries, line_array, ierr)
                  temp_name = line_array(1)
                  IF (ANY(temp_name == nonbond_list(:,is)%atom_name ) .eqv. .FALSE.) THEN
                     err_msg(1) = "Atom name and type in Intra Scaling table (" &
                                   //TRIM(temp_name)//") does not match the MCF file."
                     err_msg(2) = "Make sure the atom list in the file begins with "
                     err_msg(3) = "species 1 then 2 and so on."
                     CALL Clean_Abort(err_msg,'Read_Intra_Exclusion_Table')
                  ENDIF

                  temp_type(ia, is) = String_To_Int(line_array(2)) 

                  i_line = i_line + 1
               ENDDO
            ENDDO
            !DO is = 1, nspecies
            !   DO js = is+1, nspecies
            !      IF (ANY(temp_type(:natoms(is),is) == temp_type(:natoms(js),js))) THEN
            !         err_msg(1) = "Atom types between species"
            !         err_msg(2) = "must be unique in the intra scaling table."
            !         CALL Clean_Abort(err_msg,'Read_Intra_Exclusion_Table')
            !      ENDIF
            !   ENDDO
            !ENDDO
        ELSEIF (line_string(1:15) == '# Intra_Scaling') THEN
           ! read in values
           DO
              CALL Parse_String(intrafile_unit, i_line, 2, nbr_entries, line_array, ierr)
              IF (TRIM(line_array(2)) == 'Done_Intra_Scaling') THEN
                EXIT
              ENDIF
   
              ii = String_To_Int(line_array(1))
              jj = String_To_Int(line_array(2))
              t_atoms = 0
              DO is = 1, nspecies
                 t_atoms = t_atoms + temp_type(natoms(is),is)
                 IF ( ii <= t_atoms ) THEN
                    IF ( jj > t_atoms .OR. jj <= (t_atoms - natoms(is)) ) THEN
                       err_msg(1) = "Intra Scaling can only be done among atoms in the same molecule."
                       CALL Clean_Abort(err_msg,'Read_Intra_Exclusion_Table')
                    ENDIF
                    EXIT
                 ENDIF
              ENDDO

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

           EXIT  
        ENDIF

        i_line = i_line + 1
        CALL Read_String(intrafile_unit,line_string,ierr)
    ENDDO
    CLOSE(UNIT=intrafile_unit)

END SUBROUTINE Read_Intra_Exclusion_Table
