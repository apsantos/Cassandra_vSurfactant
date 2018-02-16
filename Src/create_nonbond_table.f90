!********************************************************************************
! CASSANDRA - Computational Atomistic Simulation Software at Notre Dame 
! for Research in Academia.
! http://molsim.wiki.zoho.com/
! Copyright (2007) University of Notre Dame.
! Authors: Ed Maginn (ed@nd.edu); Jindal Shah (jshah@nd.edu)
!********************************************************************************

!********************************************************************************
  SUBROUTINE Create_Nonbond_Table
!********************************************************************************
    ! This routine determines the number and identity of all unique atom types
    ! in the system. It then takes the associated vdw paramters of each atom
    ! and creates an interaction table using the specified mixing rule.
    ! Currently only supports LJ vdw type and LB or geometric mixing rules.

    ! The final product of this code are a bunch of matrices
    ! vdw_param1_table(itype,jtype), vdw_param2_table(itype,jtype), etc. 
    ! where for the LJ potential, vdw_param1_table(itype,jtype) is the epsilon
    ! to use between atoms of itype and jtype, and vdw_param2_table(itype,jtype)
    ! is the sigma to use between these two atom type. 

    ! Example: let's say you need to compute the LJ energy betweem atom 4 of species 2
    ! and atom 6 of species 1. 
    ! You first identify the unique atom type of these two atoms:
    ! itype = nonbond_list(4,2)%atom_type_number
    ! jtype = nonbond_list(6,1)%atom_type_number

    ! Then get the parameters:
    ! epsilon(itype, jtype) = vdw_param1_table(itype,jtype)
    ! sigma(itype,jtype) = vdw_param2_table(itype,jtype)
    ! Note that these matrices are symmetric, so 
    ! vdw_param1_table(itype,jtype) = vdw_param1_table(jtype,itype)

    ! Written: Sat Oct  6 08:58:06 MDT 2007
    ! Author: E. Maginn
    
    ! *** CALLS ***
    ! Clean_Abort
    ! 
    ! *** CALLED BY ***
    ! NVTMC_Control
    !
    ! Revision history:

    ! 04/08/09 (TR) : Removed any character evaluations.
    !
    ! 08/26/11 (JS) : Steele potential constants calculations
    !
    ! 03/27/14 Eliseo Rimoldi : custom (manually input) nonbond potentials added

!********************************************************************************
    USE Run_Variables
    USE Type_Definitions
    USE IO_Utilities
    USE File_Names
    
    IMPLICIT NONE

    INTEGER :: is, ia, itype, jtype, iset,jset,k
    REAL(DP), DIMENSION(max_nonbond_params) :: temp_param_i, temp_param_j

    !custom mixing rules
    INTEGER :: ierr,line_nbr,nbr_entries, is_1, is_2, ia_1, ia_2, itype_custom, jtype_custom
    CHARACTER(charLength) :: line_string, line_array(lineArrayLength)


!********************************************************************************
    ! APS
    ! Read in user-defined table if specified
    IF (mix_rule == 'table') THEN
       CALL Read_Nonbond_Table
       RETURN
    ENDIF

    ! Compute the number of different atom types
    CALL Compute_Atom_Types

    ! Create a character array containing the names of each unique atom type, with the index equal
    ! to the atom type number

    ALLOCATE(atom_type_list(nbr_atomtypes), Stat=AllocateStatus)

    ! allocate arrays containing vdw parameters for all interaction pairs.
    ALLOCATE(vdw_param1_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param2_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param3_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param4_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param5_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param6_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param7_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param8_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param9_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param10_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param11_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param12_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param13_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    vdw_param1_table = 0.0_DP
    vdw_param2_table = 0.0_DP
    vdw_param3_table = 0.0_DP
    vdw_param4_table = 0.0_DP
    vdw_param5_table = 0.0_DP
    vdw_param6_table = 0.0_DP
    vdw_param7_table = 0.0_DP
    vdw_param8_table = 0.0_DP
    vdw_param9_table = 0.0_DP
    vdw_param10_table = 0.0_DP
    vdw_param11_table = 0.0_DP
    vdw_param12_table = 0.0_DP
    vdw_param13_table = 0.0_DP

    ALLOCATE(rcut_vdw_mix(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    rcut_vdw_mix(:,:) = rcut_vdw(1)
    ALLOCATE(rcut_vdwsq_mix(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    rcut_vdwsq_mix(:,:) = 0.0_DP
    ALLOCATE(int_vdw_style_mix(nbr_atomtypes,nbr_atomtypes,15), Stat=AllocateStatus)
    int_vdw_style_mix(:,:,:) = .false.
    int_vdw_style_mix(:,:,int_vdw_style(1)) = .true.
    ALLOCATE(int_vdw_sum_style_mix(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    int_vdw_sum_style_mix(:,:) = int_vdw_sum_style(1)

    ! Allocate memory for total number bead types in each box
    ALLOCATE(nint_beads(nbr_atomtypes,nbr_boxes))
    
    IF (AllocateStatus .NE. 0) THEN
       err_msg = ''
       err_msg(1) = ' ERROR: Not enough memory for vdw interaction table '
       CALL Clean_Abort(err_msg,'ceate_nonbond_table')
    END IF

    atom_type_list = ""

    ! Now determine the set of vdw parameters for each type of interaction and load them into vdw_param_table
    ! This is a brute force search - but it is fast.

    WRITE(logunit,*) '*** Creating VDW interaction table ***'
    WRITE(logunit,'(A,T25,A)') 'Mixing rule used is:', mix_rule
    WRITE(logunit,*)

    ! Write header for logfile output. Specific for the vdw style
    IF (int_vdw_style(1) == vdw_lj) THEN
       WRITE(logunit,'(A6,5x,A6,2x,T20,A,T50,A)') 'Atom 1','Atom 2', 'epsilon (amu A^2/ps^2)', 'sigma (A)'
    ENDIF

    DO itype = 1, nbr_atomtypes

       DO jtype = 1, nbr_atomtypes


        IF (mix_rule /= 'custom') THEN

            ! Flags that are tripped when a particular atomtype parameter set is located
            iset = 0
            jset = 0

            DO is = 1, nspecies

               DO ia = 1, natoms(is)

                  IF (iset == 0) THEN
                     ! Search for atomtype i because the parameters have not been found

                     IF (nonbond_list(ia,is)%atom_type_number == itype) THEN
                        ! This type of atom has been found. Grab its vdw parameters

                        DO k=1, nbr_vdw_params(is)

                           temp_param_i(k) = nonbond_list(ia,is)%vdw_param(k)
                         
                        ENDDO

                        ! This atomtype parameters are now determined. Trip flag
                        iset = 1
                      
                        atom_type_list(itype) = nonbond_list(ia,is)%atom_name
  
                     ENDIF

                  ENDIF

                  IF(jset == 0) THEN

                  ! Now search for parameters of type j

                     IF (nonbond_list(ia,is)%atom_type_number == jtype) THEN

                        DO k=1, nbr_vdw_params(is)
                         
                           temp_param_j(k) = nonbond_list(ia,is)%vdw_param(k)

                        ENDDO
                      
                        ! This atom type has been located
                        jset = 1

                        atom_type_list(jtype) = nonbond_list(ia,is)%atom_name

                     ENDIF

                  ENDIF

               ENDDO

            ENDDO

            ! atom types i and j have been found and vdw parameters loaded into temporary arrays.
            ! Apply mixing rules and load vdw_table


            ! Adapt this to other potential types and mixing rules. Custom mixing rules can 
            ! be created manually and then this routine should be bypassed.


            IF (int_vdw_style(1) == vdw_lj   .or. & 
                int_vdw_style(1) == vdw_lj96 .or. &
                int_vdw_style(1) == vdw_lj124.or. &
                int_vdw_style(1) == vdw_mie) THEN
               ! There are two vdw parameters
               int_vdw_style_mix(itype, jtype,int_vdw_style(1)) = .true.

               ! Set LJ epsilon
               IF ( (temp_param_i(1) <= tiny_number) .OR. (temp_param_j(1) <= tiny_number) ) THEN
               ! for parameters with zero, avoid overflow and set to zero
                  vdw_param1_table(itype,jtype) = 0.0_DP

               ELSE

                  ! Use specified mixing rule
                  ! LB mixing rule: epsij = (epsi * epsj)^(1/2); sigmaij = 1/2 (sigmai + sigmaj)
                  IF (mix_rule == 'LB') &
                       vdw_param1_table(itype,jtype) = dsqrt(temp_param_i(1)*temp_param_j(1))

                  ! geometric mixing rule: epsij = (epsi * epsj)^(1/2); sigmaij = (sigmai * sigmaj)^(1/2)
                  IF (mix_rule == 'geometric')  &
                       vdw_param1_table(itype,jtype) = dsqrt(temp_param_i(1)*temp_param_j(1))
                
               ENDIF

               ! Set LJ sigma
               IF ( (temp_param_i(2) <= 1.0E-05) .OR. (temp_param_j(2) <= 1.0E-05) ) THEN
                  vdw_param2_table(itype,jtype) = 0.0_DP

               ELSE

                  IF (mix_rule == 'LB') &
                       vdw_param2_table(itype,jtype) = (temp_param_i(2) + temp_param_j(2)) * 0.5
                  IF (mix_rule == 'geometric') &
                       vdw_param2_table(itype,jtype) = dsqrt(temp_param_i(2) * temp_param_j(2))

               ENDIF

              ! Report parameters to logfile. Format is specific to vdw type. Add others here if 
              ! other than LJ potential is used.
                WRITE(logunit,'(A6,5x,A6,2x,T20,f10.4,T50,f10.4)') &
                      atom_type_list(itype), atom_type_list(jtype), &
                      vdw_param1_table(itype,jtype), vdw_param2_table(itype,jtype)

            ELSE IF (int_vdw_style(1) == vdw_yukawa) THEN
               ! There are two vdw parameters
               int_vdw_style_mix(itype, jtype, int_vdw_style(1)) = .true.

               ! Set epsilon
               IF ( (temp_param_i(1) <= tiny_number) .OR. (temp_param_j(1) <= tiny_number) ) THEN
               ! for parameters with zero, avoid overflow and set to zero
                  vdw_param8_table(itype,jtype) = 0.0_DP

               ELSE
                  ! Use specified mixing rule
                  ! LB mixing rule: epsij = (epsi * epsj)^(1/2); sigmaij = 1/2 (sigmai + sigmaj)
                  IF (mix_rule == 'LB') &
                       vdw_param8_table(itype,jtype) = dsqrt(temp_param_i(1)*temp_param_j(1))

                  ! geometric mixing rule: epsij = (epsi * epsj)^(1/2); sigmaij = (sigmai * sigmaj)^(1/2)
                  IF (mix_rule == 'geometric')  &
                       vdw_param8_table(itype,jtype) = dsqrt(temp_param_i(1)*temp_param_j(1))
                
               ENDIF

               ! Set kappa
               IF ( (temp_param_i(2) <= 1.0E-05) .OR. (temp_param_j(2) <= 1.0E-05) ) THEN
                  vdw_param9_table(itype,jtype) = 0.0_DP

               ELSE

                  IF (mix_rule == 'LB') &
                       vdw_param9_table(itype,jtype) = (temp_param_i(2) + temp_param_j(2)) * 0.5
                  IF (mix_rule == 'geometric') &
                       vdw_param9_table(itype,jtype) = dsqrt(temp_param_i(2) * temp_param_j(2))

               ENDIF

              ! Report parameters to logfile. Format is specific to vdw type. Add others here if 
              ! other than LJ potential is used.
                WRITE(logunit,'(A6,5x,A6,2x,T20,f10.4,T50,f10.4)') &
                      atom_type_list(itype), atom_type_list(jtype), &
                      vdw_param9_table(itype,jtype), vdw_param10_table(itype,jtype)

            ELSE IF (int_vdw_style(1) == vdw_screen) THEN
               ! There are two vdw parameters
               int_vdw_style_mix(itype, jtype, int_vdw_style(1)) = .true.

               ! Set kappa
               IF ( (temp_param_i(1) <= tiny_number) .OR. (temp_param_j(1) <= tiny_number) ) THEN
               ! for parameters with zero, avoid overflow and set to zero
                  vdw_param12_table(itype,jtype) = 0.0_DP

               ELSE
                  IF (mix_rule == 'LB') &
                       vdw_param12_table(itype,jtype) = dsqrt(temp_param_i(1)*temp_param_j(1))

                  ! geometric mixing rule: epsij = (epsi * epsj)^(1/2); sigmaij = (sigmai * sigmaj)^(1/2)
                  IF (mix_rule == 'geometric')  &
                       vdw_param12_table(itype,jtype) = dsqrt(temp_param_i(1)*temp_param_j(1))
                
               ENDIF

               ! Set LJ sigma
               IF ( (temp_param_i(2) <= 1.0E-05) .OR. (temp_param_j(2) <= 1.0E-05) ) THEN
                  vdw_param13_table(itype,jtype) = 0.0_DP

               ELSE

                  IF (mix_rule == 'LB') &
                       vdw_param13_table(itype,jtype) = (temp_param_i(2) + temp_param_j(2)) * 0.5
                  IF (mix_rule == 'geometric') &
                       vdw_param13_table(itype,jtype) = dsqrt(temp_param_i(2) * temp_param_j(2))

               ENDIF

              ! Report parameters to logfile. Format is specific to vdw type. Add others here if 
              ! other than LJ potential is used.
                WRITE(logunit,'(A6,5x,A6,2x,T20,f10.4,T50,f10.4)') &
                      atom_type_list(itype), atom_type_list(jtype), &
                      vdw_param12_table(itype,jtype), vdw_param13_table(itype,jtype)

            ELSE IF (int_vdw_style(1) == vdw_sw) THEN
               ! There are two vdw parameters
               int_vdw_style_mix(itype, jtype, int_vdw_style(1)) = .true.

               ! Set LJ epsilon
               IF ( (temp_param_i(1) <= tiny_number) .OR. (temp_param_j(1) <= tiny_number) ) THEN
                  vdw_param10_table(itype,jtype) = 0.0_DP
               ! for parameters with zero, avoid overflow and set to zero

               ELSE
                  ! Use specified mixing rule
                  ! LB mixing rule: epsij = (epsi * epsj)^(1/2); sigmaij = 1/2 (sigmai + sigmaj)
                  IF (mix_rule == 'LB') &
                       vdw_param10_table(itype,jtype) = dsqrt(temp_param_i(1)*temp_param_j(1))

                  ! geometric mixing rule: epsij = (epsi * epsj)^(1/2); sigmaij = (sigmai * sigmaj)^(1/2)
                  IF (mix_rule == 'geometric')  &
                       vdw_param10_table(itype,jtype) = dsqrt(temp_param_i(1)*temp_param_j(1))
                
               ENDIF

               ! Set LJ sigma
               IF ( (temp_param_i(2) <= 1.0E-05) .OR. (temp_param_j(2) <= 1.0E-05) ) THEN
                  vdw_param11_table(itype,jtype) = 0.0_DP

               ELSE
                  IF (mix_rule == 'LB') &
                       vdw_param11_table(itype,jtype) = (temp_param_i(2) + temp_param_j(2)) * 0.5
                  IF (mix_rule == 'geometric') &
                       vdw_param11_table(itype,jtype) = dsqrt(temp_param_i(2) * temp_param_j(2))

               ENDIF

              ! Report parameters to logfile. Format is specific to vdw type. Add others here if 
              ! other than LJ potential is used.
                WRITE(logunit,'(A6,5x,A6,2x,T20,f10.4,T50,f10.4)') &
                      atom_type_list(itype), atom_type_list(jtype), &
                      vdw_param10_table(itype,jtype), vdw_param11_table(itype,jtype)

            ENDIF


          ELSE !custom mixing rule if
  
            REWIND(inputunit)

            ierr = 0
            line_nbr = 0

            DO
                line_nbr = line_nbr + 1

                CALL Read_String(inputunit,line_string,ierr)


                IF (ierr .NE. 0) THEN
                    err_msg = ""
                    err_msg(1) = "Error reading mixing rules1."
                    CALL Clean_Abort(err_msg,'Get_Mixing_Rules')
                END IF

                IF (line_string(1:13) == '# Mixing_Rule') THEN
                    READ(inputunit,*)
                    ! Assign the first entry on the line to the mixing rule
                    DO is_1 = 1, nspecies
                        DO ia_1 = 1, natoms(is_1)
                            itype_custom = nonbond_list(ia_1,is_1)%atom_type_number
                            DO is_2 = 1, nspecies
                                DO ia_2 = 1, natoms(is_2)
                                    CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
                                    jtype_custom = nonbond_list(ia_2,is_2)%atom_type_number
                                    !Convert epsilon to atomic units amu A^2/ps^2 
                                    vdw_param1_table(itype_custom,jtype_custom) = &
                                                            kboltz * String_To_Double(line_array(3))
                                    vdw_param2_table(itype_custom,jtype_custom) = &
                                                            String_To_Double(line_array(4))
                                    WRITE(logunit,'(A6,5x,A6,2x,T20,f10.4,T50,f10.4)') &
                                             atom_type_list(itype_custom), atom_type_list(jtype_custom), &
                                             vdw_param1_table(itype_custom,jtype_custom), &
                                             vdw_param2_table(itype_custom,jtype_custom)

                                END DO
                            END DO
                        END DO
                    END DO
                    RETURN
                END IF
            END DO
          END IF
       ENDDO
    
   ENDDO
 
 WRITE(logunit,*)
 WRITE(logunit,*) '*** Completed construction of VDW interaction table ***'
 WRITE(logunit,*)    


END SUBROUTINE Create_Nonbond_Table

SUBROUTINE Read_Nonbond_Table
!********************************************************************************
    ! Written: Jun 25 2015
    ! Author: A. Santos
    
    ! 
    ! *** CALLED BY ***
    ! Create_Nonbond_Table
    !
    ! Revision history:
    !

!********************************************************************************
    USE Run_Variables
    USE Type_Definitions
    USE IO_Utilities
    USE File_Names
    
    IMPLICIT NONE

    INTEGER :: i, is, ia, itype, jtype, tot_natoms, iatom

    !custom mixing rules
    INTEGER :: ierr,nbr_entries
    INTEGER :: i_line, n_params, cur_line
    CHARACTER(charLength) :: line_string, line_array(lineArrayLength)
    CHARACTER(charLength) :: temp_name, pot_type, temp_name_list(30)
    INTEGER :: temp_type_list(30), temp_type, ncheck
  !********************************************************************************

!********************************************************************************
    ! Open mixfile and find the line where the data begins
    OPEN(UNIT=mixfile_unit,FILE=mixfile_name,STATUS="OLD",IOSTAT=openstatus,ACTION="READ")
    CALL Read_String(mixfile_unit,line_string,ierr)
    IF (ierr .NE. 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading mixinf rules."
        CALL Clean_Abort(err_msg,'Get_Mixing_Rules')
    END IF
    i_line = 1
    tot_natoms = SUM(natoms(:))
    temp_type_list = -100
    temp_name_list = ""
    nbr_atomtypes = 0
    DO 
        IF (line_string(1:12) == '# Atom_Types') THEN
            DO iatom = 1, tot_natoms
                CALL Parse_String(mixfile_unit, i_line, 1, nbr_entries, line_array, ierr)
                IF (TRIM(line_array(2)) == 'Done_Atom_Types') THEN
                   EXIT
                ENDIF
                temp_name = line_array(1)
                temp_type = String_To_Int( line_array(2) )
                IF ( ALL(temp_type_list /= temp_type) ) THEN
                    nbr_atomtypes = nbr_atomtypes + 1
                    temp_type_list(nbr_atomtypes) = temp_type
                    temp_name_list(nbr_atomtypes) = temp_name
                END IF
                !Check which speicies and atom corresponds to this name
                ncheck = 0
                DO is = 1, nspecies
                    DO ia = 1, natoms(is)
                        IF (temp_name == nonbond_list(ia,is)%atom_name) THEN
                            nonbond_list(ia,is)%atom_type_number = temp_type
                            ncheck = 1
                            
                        ENDIF
                    ENDDO
                ENDDO
                IF (ncheck == 0) THEN
                    err_msg(1) = "Atom name in mixing table ("//TRIM(temp_name)//") does not match any name in MCF file."
                    CALL Clean_Abort(err_msg,'Get_Mixing_Rules')
                ENDIF
            ENDDO
            CALL Read_String(mixfile_unit,line_string,ierr)
            EXIT
        ENDIF
        i_line = i_line + 1
        CALL Read_String(mixfile_unit,line_string,ierr)
    ENDDO

    ! Create a character array containing the names of each unique atom type, with the index equal
    ! to the atom type number

    ALLOCATE(atom_type_list(nbr_atomtypes), Stat=AllocateStatus)

    DO itype = 1, nbr_atomtypes
       atom_type_list(itype) = temp_name_list(itype)
    END DO

    ! allocate arrays containing vdw parameters for all interaction pairs.
    ALLOCATE(vdw_param1_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param2_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param3_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param4_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param5_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param6_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param7_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param8_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param9_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param10_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param11_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param12_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    ALLOCATE(vdw_param13_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    vdw_param1_table = 0.0_DP
    vdw_param2_table = 0.0_DP
    vdw_param3_table = 0.0_DP
    vdw_param4_table = 0.0_DP
    vdw_param5_table = 0.0_DP
    vdw_param6_table = 0.0_DP
    vdw_param7_table = 0.0_DP
    vdw_param8_table = 0.0_DP
    vdw_param9_table = 0.0_DP
    vdw_param10_table = 0.0_DP
    vdw_param11_table = 0.0_DP
    vdw_param12_table = 0.0_DP
    vdw_param13_table = 0.0_DP

    ALLOCATE(rcut_vdw_mix(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    rcut_vdw_mix(:,:) = rcut_vdw(1)
    ALLOCATE(rcut_vdwsq_mix(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    rcut_vdwsq_mix(:,:) = 0.0_DP
    ALLOCATE(int_vdw_style_mix(nbr_atomtypes,nbr_atomtypes,15), Stat=AllocateStatus)
    int_vdw_style_mix(:,:,:) = .false.
    int_vdw_style_mix(:,:,int_vdw_style(1)) = .true.
    ALLOCATE(int_vdw_sum_style_mix(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    int_vdw_sum_style_mix(:,:) = int_vdw_sum_style(1)

    ! Allocate memory for total number bead types in each box
    ALLOCATE(nint_beads(nbr_atomtypes,nbr_boxes))
    
    IF (AllocateStatus .NE. 0) THEN
       err_msg = ''
       err_msg(1) = ' ERROR: Not enough memory for vdw interaction table '
       CALL Clean_Abort(err_msg,'ceate_nonbond_table')
    END IF

    ! Now determine the set of vdw parameters for each type of interaction and load them into vdw_param_table
    ! This is a brute force search - but it is fast.

    WRITE(logunit,*) '*** Creating VDW interaction table ***'
    WRITE(logunit,'(A,T25,A)') 'Mixing rule used is:', mix_rule
    WRITE(logunit,*)

    REWIND(inputunit)

    DO 
        IF (line_string(1:15) == '# Mixing_Values') THEN
            EXIT
        ENDIF

        i_line = i_line + 1
        CALL Read_String(mixfile_unit,line_string,ierr)
    ENDDO
    n_params = INT(nbr_atomtypes * (nbr_atomtypes + 1 ) / 2)

    ! Read in from the table file
line_loop: DO cur_line = i_line, n_params + i_line - 1
        CALL Parse_String(mixfile_unit, cur_line, 2, nbr_entries, line_array, ierr)
        itype = String_To_Int( line_array(1) )
        jtype = String_To_Int( line_array(2) )
        IF (itype > nbr_atomtypes) THEN
            WRITE(*,'(A,T25,A)') '    WARNING:'
            WRITE(*,'(A,I3,2A)') 'type (', itype, ') in the mix file', mixfile_name,&
                                 'is too large based on # Atom_Types info and number of species, skipping'
            CYCLE line_loop
        END IF
        IF (jtype > nbr_atomtypes) THEN
            WRITE(*,'(A,T25,A)') '    WARNING'
            WRITE(*,'(A,I3,2A)') 'type (', jtype, ') in the mix file', mixfile_name,&
                                 'is too large based on # Atom_Types info and number of species, skipping'
            CYCLE line_loop
        END IF

        int_vdw_style_mix(itype,jtype,:) = .false.
        int_vdw_style_mix(jtype,itype,:) = .false.
                
        DO i = 3, SIZE(line_array)
            pot_type = line_array(i)

            IF (pot_type == 'LJ' .or. pot_type == 'LJ126' .or. &
                pot_type == 'LJ124' .or. &
                pot_type == 'LJ96') THEN

                IF (pot_type == 'LJ' .or. pot_type == 'LJ126') THEN
                    int_vdw_style_mix(itype,jtype,vdw_lj) = .true.
                    int_vdw_style_mix(jtype,itype,vdw_lj) = .true.
                ELSEIF ( pot_type == 'LJ124' ) THEN
                    int_vdw_style_mix(itype,jtype,vdw_lj124) = .true.
                    int_vdw_style_mix(jtype,itype,vdw_lj124) = .true.
                ELSEIF ( pot_type == 'LJ96') THEN
                    int_vdw_style_mix(itype,jtype,vdw_lj96) = .true.
                    int_vdw_style_mix(jtype,itype,vdw_lj96) = .true.
                ENDIF
                vdw_param1_table(itype,jtype) = String_To_Double(line_array(i+1))
                vdw_param2_table(itype,jtype) = String_To_Double(line_array(i+2))
                vdw_param1_table(jtype,itype) = vdw_param1_table(itype,jtype)
                vdw_param2_table(jtype,itype) = vdw_param2_table(itype,jtype)

            ELSEIF (pot_type == 'HYDR') THEN
                int_vdw_style_mix(itype,jtype,vdw_hydra) = .true.
                int_vdw_style_mix(jtype,itype,vdw_hydra) = .true.
                vdw_param3_table(itype,jtype) = String_To_Double(line_array(i+1))
                vdw_param4_table(itype,jtype) = String_To_Double(line_array(i+2))
                vdw_param5_table(itype,jtype) = String_To_Double(line_array(i+3))
                vdw_param3_table(jtype,itype) = vdw_param3_table(itype,jtype)
                vdw_param4_table(jtype,itype) = vdw_param4_table(itype,jtype)
                vdw_param5_table(jtype,itype) = vdw_param5_table(itype,jtype)

            ELSEIF (pot_type == 'CORR') THEN
                int_vdw_style_mix(itype,jtype,vdw_corr) = .true.
                int_vdw_style_mix(jtype,itype,vdw_corr) = .true.
                vdw_param6_table(itype,jtype) = String_To_Double(line_array(i+1))
                vdw_param7_table(itype,jtype) = String_To_Double(line_array(i+2))
                vdw_param6_table(jtype,itype) = vdw_param6_table(itype,jtype)
                vdw_param7_table(jtype,itype) = vdw_param7_table(itype,jtype)

            ELSEIF (pot_type == 'Yukawa') THEN
                int_vdw_style_mix(itype,jtype,vdw_yukawa) = .true.
                int_vdw_style_mix(jtype,itype,vdw_yukawa) = .true.
                vdw_param8_table(itype,jtype) = String_To_Double(line_array(i+1))
                vdw_param9_table(itype,jtype) = String_To_Double(line_array(i+2))
                vdw_param8_table(jtype,itype) = vdw_param8_table(itype,jtype)
                vdw_param9_table(jtype,itype) = vdw_param9_table(itype,jtype)

            ELSEIF (pot_type == 'SCR') THEN
                int_vdw_style_mix(itype,jtype,vdw_screen) = .true.
                int_vdw_style_mix(jtype,itype,vdw_screen) = .true.
                ! line_array(i+1) is the number of added salt in the solution, get from cmc0 + c_salt
                ! convert from mM to num density
                vdw_param12_table(itype,jtype) = SQRT( (String_To_Double(line_array(i+1)) * navogadro / m3_to_A3 ) * &
                                                       (charge_factor(1) * 8.0 * PI) / beta(1) )
                vdw_param13_table(itype,jtype) = String_To_Double(line_array(i+2)) 
                vdw_param12_table(jtype,itype) = vdw_param12_table(itype,jtype)
                vdw_param13_table(jtype,itype) = vdw_param13_table(itype,jtype)

            ELSEIF (pot_type == 'SW') THEN
                int_vdw_style_mix(itype,jtype,vdw_sw) = .true.
                int_vdw_style_mix(jtype,itype,vdw_sw) = .true.
                vdw_param10_table(itype,jtype) = String_To_Double(line_array(i+1))
                vdw_param11_table(itype,jtype) = String_To_Double(line_array(i+2))
                vdw_param10_table(jtype,itype) = vdw_param10_table(itype,jtype)
                vdw_param11_table(jtype,itype) = vdw_param11_table(itype,jtype)

            ENDIF

        ENDDO
    ENDDO line_loop

    ! read the vdw sum style and cutoffs if there
    i_line = 1

    DO 
        CALL Read_String(mixfile_unit,line_string,ierr)
        IF (ierr < 0) EXIT
        IF (line_string(1:11) == '# VDW_Style') THEN
           ! read in values
           DO
              i_line = i_line + 1
              CALL Parse_String(mixfile_unit, i_line, 2, nbr_entries, line_array, ierr)
              IF (TRIM(line_array(2)) == 'Done_VDW_Style') THEN
                EXIT
              ENDIF
   
              itype = String_To_Int( line_array(1) )
              jtype = String_To_Int( line_array(2) )

              IF (itype > nbr_atomtypes) THEN
                  WRITE(*,'(A,T25,A)') '    WARNING'
                  WRITE(*,'(A,I3,2A)') 'type (', itype, ') in the mix file', mixfile_name,&
                                       'is too large based on # Atom_Types info and number of species, skipping'
                  CYCLE 
              END IF
              IF (jtype > nbr_atomtypes) THEN
                  WRITE(*,'(A,T25,A)') '    WARNING'
                  WRITE(*,'(A,I3,2A)') 'type (', jtype, ') in the mix file', mixfile_name,&
                                       'is too large based on # Atom_Types info and number of species, skipping'
                  CYCLE 
              END IF

              rcut_vdw_mix(itype,jtype) = String_To_Double(line_array(4))
              rcut_vdw_mix(jtype,itype) = String_To_Double(line_array(4))

              IF (line_array(3) == 'cut') THEN
                 int_vdw_sum_style_mix(itype,jtype) = vdw_cut
              ELSE IF (line_array(3) == 'cut_shift') THEN
                 int_vdw_sum_style_mix(itype,jtype) = vdw_cut_shift
              ELSE IF (line_array(3) == 'cut_switch') THEN
                 int_vdw_sum_style_mix(itype,jtype) = vdw_cut_switch
              ELSE IF (line_array(3) == 'CHARMM') THEN
                 int_vdw_sum_style_mix(itype,jtype) = vdw_charmm
              ENDIF
              int_vdw_sum_style_mix(jtype,itype) = int_vdw_sum_style_mix(itype,jtype)
                  
           ENDDO
        ENDIF
        i_line = i_line + 1
    ENDDO


    ! Write output
    WRITE(logunit,'(A)') 'itype jtype vdw_param 1 2 3 4 5 6 7 8 9 10 11 12 13'
    DO itype = 1, nbr_atomtypes
        DO jtype = itype, nbr_atomtypes
            WRITE(logunit,'(2I3,13f13.4)') &
                 itype, jtype, &
                 vdw_param1_table(itype,jtype), vdw_param2_table(itype,jtype), &
                 vdw_param3_table(itype,jtype), vdw_param4_table(itype,jtype), &
                 vdw_param5_table(itype,jtype), vdw_param6_table(itype,jtype), &
                 vdw_param7_table(itype,jtype), vdw_param8_table(itype,jtype), &
                 vdw_param9_table(itype,jtype), vdw_param10_table(itype,jtype), &
                 vdw_param11_table(itype,jtype), vdw_param12_table(itype,jtype), &
                 vdw_param13_table(itype,jtype)
        ENDDO
    ENDDO

 WRITE(logunit,*)
 WRITE(logunit,*) '*** Completed construction of VDW interaction table ***'
 WRITE(logunit,*)    


END SUBROUTINE Read_Nonbond_Table
! APS
! Compute the number of different atom types.
! Originally in Create_Nonbond_Table.  
! Moved for clarity and to use in Read_Nonbond_Table.
SUBROUTINE Compute_Atom_Types
!********************************************************************************
    USE Run_Variables
    USE Type_Definitions
    USE IO_Utilities
    USE File_Names
    
    IMPLICIT NONE

    LOGICAL :: repeat_type
    CHARACTER(6), DIMENSION(:), ALLOCATABLE :: temp_atomtypes ! same dimension as atom_name
    INTEGER :: ii, is, ia

!********************************************************************************
    ALLOCATE(temp_atomtypes(1000), Stat=AllocateStatus)
    IF (AllocateStatus .NE. 0) THEN
       err_msg = ''
       err_msg(1) = ' ERROR: Not enough memory for temp_atomtypes '
       CALL Clean_Abort(err_msg,'create_nonbond_table')
    END IF

    ! Initialize the number of atom types and the temp_atomtypes
    nbr_atomtypes = 0
    temp_atomtypes = ''

    DO is = 1, nspecies
       
       DO ia = 1, natoms(is)
          
          repeat_type = .FALSE.

          !----------------------------------------------------------------
          ! Determine whether the atomtype has already been accounted for
          !----------------------------------------------------------------
          IF ((is .EQ. 1).AND.(ia .EQ. 1)) THEN
             ! If this is the first atomtype obtained, store it in temp_atomtypes
             nbr_atomtypes = nbr_atomtypes + 1
             temp_atomtypes(nbr_atomtypes) = nonbond_list(ia,is)%atom_name

             ! Store the unique number identifier
             nonbond_list(ia,is)%atom_type_number = nbr_atomtypes
            
          ELSE
             ! Loop over all current atomtypes to check if the atomtype has been accounted for
             !   If so, turn the fatomtype flag to true
             DO ii = 1, nbr_atomtypes
                IF(nonbond_list(ia,is)%atom_name .EQ. temp_atomtypes(ii)) THEN

                   ! This atom name is already present. Do not advance counter
                   repeat_type = .TRUE.

                   ! Store the unique number identifier
                   nonbond_list(ia,is)%atom_type_number = ii

                ENDIF
             ENDDO

             ! If the atomtype has not been accounted for, add it to the temp_atomtypes list
             IF(.NOT.repeat_type) THEN
                nbr_atomtypes = nbr_atomtypes + 1
                temp_atomtypes(nbr_atomtypes) = nonbond_list(ia,is)%atom_name

                ! Store the unique number identifier
                nonbond_list(ia,is)%atom_type_number = nbr_atomtypes
                
             ENDIF

          ENDIF

       ENDDO

    ENDDO

    ! Write the number of different atom types to the screen and logfile
    WRITE(logunit,'(A)') &
         '  There are '//TRIM(Int_To_String(nbr_atomtypes))//' different atom types in the system '
    DO ii = 1, nbr_atomtypes
       WRITE(logunit,'(3x,I3,2x,A6)') ii, temp_atomtypes(ii)
    ENDDO

    WRITE(logunit,*)

    DO is=1,nspecies
       WRITE(logunit,*)
       WRITE(logunit,'(A,T25,I3,3x,A)') 'species number and name:',is, molfile_name(is)
       WRITE(logunit,*) 'Name      number'
       WRITE(logunit,*) '------    ------'

       DO ia = 1, natoms(is)
          WRITE(logunit,'(A6,T10,I4)') nonbond_list(ia,is)%atom_name, &
               nonbond_list(ia,is)%atom_type_number
       ENDDO

    ENDDO

    WRITE(logunit,*)
    WRITE(logunit,*)

    IF (ALLOCATED(temp_atomtypes)) DEALLOCATE(temp_atomtypes)

END SUBROUTINE Compute_Atom_Types
