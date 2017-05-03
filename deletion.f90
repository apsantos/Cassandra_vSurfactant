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

SUBROUTINE Deletion(this_box)
  
  !*****************************************************************************
  !
  ! PURPOSE: attempt to delete a molecule that was inserted via 
  !          configurational bias monte carlo
  !
  ! Called by
  !
  !  gcmc_driver
  ! 
  ! Revision History
  !
  !   12/10/13 : Beta Release
  !   Version 1.1
  !     04/21/15  Corrected acceptance criteria
  !     05/05/15  Documented this code
  !
  ! DESCRIPTION: This subroutine performs the following steps:
  !
  ! Step 1) Select a species with uniform probability
  ! Step 2) Select a molecule with uniform probability
  ! Step 3) Calculate the bias probability for the reverse insertion move
  ! Step 4) Calculate the change in potential energy if the molecule is deleted
  ! Step 5) Accept or reject the move
  !
  !*****************************************************************************
  
  USE Type_Definitions
  USE Run_Variables
  USE Random_Generators
  USE Simulation_Properties
  USE Energy_Routines
  USE Fragment_Growth
  USE IO_Utilities

  !*****************************************************************************
  ! Declare and Initialize Variables
  !*****************************************************************************

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(INOUT) :: this_box ! attempt to delete a molecule in this_box

  ! Local declarations
  INTEGER :: i, i_type               ! atom indices
  INTEGER :: ifrag                   ! fragment indices
  INTEGER :: im(nspecies), alive(nspecies)               ! molecule indices
  INTEGER :: is, is_rand, is_counter ! species indices
  INTEGER :: kappa_tot
  INTEGER, ALLOCATABLE :: frag_order(:)
  INTEGER :: k, position
  INTEGER :: tn1, tn2, n1, n2, nplocal, npair, dn

  REAL(DP) :: ppt, pp(n_insertable), randnpair, loc_chem_pot
  REAL(DP) :: delta_e, dblocal
  REAL(DP) :: E_bond, E_angle, E_dihedral, E_improper
  REAL(DP) :: f_bond, f_angle, f_dihedral, f_improper
  REAL(DP) :: E_intra_vdw, E_intra_qq
  REAL(DP) :: E_inter_vdw, E_inter_qq
  REAL(DP) :: E_reciprocal_move, E_self_move, E_lrc
  REAL(DP) :: nrg_ring_frag_tot
  REAL(DP) :: ln_pacc, P_seq, P_bias, this_lambda
  REAL(DP) :: E_intra_vdw_igas, E_intra_qq_igas
  REAL(DP) :: fp_bias, fp_seq, f_ring, adden, f_vdw_igas, f_qq_igas
  REAL(DP) :: f_inter_vdw, f_inter_qq, pair_vdw, pair_qq
  REAL(DP) :: f_intra_vdw, f_intra_qq, f_reciprocal, f_self_diff

  LOGICAL :: inter_overlap(nspecies), cbmc_overlap(nspecies), intra_overlap(nspecies)
  LOGICAL :: accept, accept_or_reject, poverlap, isgas, isfrag

  ! Initialize variables
  ln_pacc = 0.0_DP
  P_seq = 1.0_DP
  P_bias = 1.0_DP
  this_lambda = 1.0_DP
  nrg_ring_frag_tot = 0.0_DP
  inter_overlap(:) = .FALSE.
  cbmc_overlap(:) = .FALSE.
  intra_overlap(:) = .FALSE.
  fp_bias = 1.0_DP
  fp_seq = 1.0_DP
  f_ring = 0.0_DP
  f_inter_vdw = 0.0_DP
  f_inter_qq = 0.0_DP
  pair_vdw = 0.0_DP
  pair_qq = 0.0_DP
  poverlap = .FALSE.
  f_intra_vdw = 0.0_DP
  f_intra_qq = 0.0_DP
  f_reciprocal = 0.0_DP
  f_self_diff = 0.0_DP
  isgas = .FALSE.
  isfrag = .FALSE.
  adden = 0.0_DP
  f_bond = 0.0_DP
  f_angle = 0.0_DP
  f_dihedral = 0.0_DP
  f_improper = 0.0_DP

  !*****************************************************************************
  ! Step 1) Select a species with uniform probability
  !*****************************************************************************
  !
  ! All species may not be insertable. For example, in a simulation of dilute
  ! water (species 3) and CO2 (species 4) in an ionic liquid (species 1 and 2), 
  ! the number of ionic liquid molecules may be fixed and only the numbers of
  ! water and CO2 allowed to fluctuate. First, choose a random integer between 1
  ! and the number of insertable species, nspec_insert:

  if(any(species_list(:)%pair_insert) .eqv. .TRUE.) then
  do i = 1, n_insertable
  tn1 = ins_species_index(i,1)
  tn2 = ins_species_index(i,2)
  ppt = prob_species_ins_pair(tn1,tn2)
  if (i == 1) then
     pp(i) = ppt
  else
     pp(i) = ppt + pp(i-1)
  endif
  enddo

  randnpair = rranf()
  do i = n_insertable, 1, -1
  if(randnpair .LE. pp(i)) then
     n1 = ins_species_index(i,1)
     n2 = ins_species_index(i,2)
     npair = i
  endif
  enddo
  else
  is_rand = INT(rranf() * nspec_insert) + 1

  is_counter = 0
  DO is = 1, nspecies
     IF(species_list(is)%int_species_type == int_sorbate) THEN
        is_counter = is_counter + 1
     END IF
     IF(is_counter == is_rand) EXIT ! exit the loop when 'is' has been found
  END DO
  n1 = is
  n2 = is
  endif

  ! Cannot delete a molecule if there aren't any in the box
  IF (nmols(n1,this_box) == 0 .OR. nmols(n2,this_box) == 0) RETURN

  ! Now that a deletion will be attempted, we need to do some bookkeeping:
  !  * Increment the counters to compute success ratios

  dn = n2 - n1
  if (dn .EQ. 0) dn = 1

  do is = n1, n2, dn
  ntrials(is,this_box)%deletion = ntrials(is,this_box)%deletion + 1  
  enddo

  tot_trials(this_box) = tot_trials(this_box) + 1

  !*****************************************************************************
  ! Step 2) Select a molecule with uniform probability
  !*****************************************************************************
  !

  do is = n1, n2, dn
  im(is) = INT(rranf() * nmols(is,this_box)) + 1
  CALL Get_Index_Molecule(this_box,is,im(is),alive(is)) ! sets the value of 'alive(is)'

  ! Save the coordinates of 'alive(is)' because Build_Molecule will erase them if
  ! cbmc_overlap(is) is tripped.

  CALL Save_Old_Cartesian_Coordinates(alive(is),is)
  enddo

  ! Compute the energy of the molecule
  
  !*****************************************************************************
  ! Step 3) Calculate the bias probability for the reverse insertion move
  !*****************************************************************************
  !
  ! The bias probability, P_bias, of the reverse insertion move is required to
  ! calculate the probability of accepting the deletion. P_bias will be 
  ! calculated using the following procedure:
  ! 
  !   3.1) Select the order to insert fragments, with probability P_seq
  !   3.2) Select kappa_ins - 1 trial coordinates, each with uniform probability
  !   3.3) Calculate the probability of the fragment's current COM
  !   3.4) For each additional fragment:
  !          a) Select kappa_dih - 1 trial dihedrals, each with uniform 
  !             probability
  !          b) Calculate the probability of the fragment's current dihedral
  !
  ! These steps are implemented in the subroutine Build_Molecule
  do is = n1, n2, dn
  P_seq = 1.0_DP
  P_bias = 1.0_DP
  nrg_ring_frag_tot = 0.0_DP
  IF(species_list(is)%fragment .AND. &
     (species_list(is)%int_insert .NE. int_igas)) THEN

     ! Build_Molecule places the first fragment, then calls Fragment_Placement 
     ! to place the additional fragments
     del_flag = .TRUE.      ! Don't change the coordinates of 'alive(is)'
     get_fragorder = .TRUE. !
     ALLOCATE(frag_order(nfragments(is)))
     CALL Build_Molecule(alive(is),is,this_box,frag_order,this_lambda, &
             P_seq, P_bias, nrg_ring_frag_tot, cbmc_overlap(is))
     DEALLOCATE(frag_order)
     
     ! cbmc_overlap(is) will only trip if the molecule being deleted had bad
     ! contacts
     IF (cbmc_overlap(is)) THEN
        WRITE(*,*)
        WRITE(*,*) 'Warning....energy overlap(is) detected in old configuration in deletion.f90'
        WRITE(*,*) 'molecule, species', alive(is), is
        WRITE(*,*)

        ! Revert to the COM and Eulerian angles for the molecule
        CALL Revert_Old_Cartesian_Coordinates(alive(is),is)
        atom_list(1:natoms(is),alive(is),is)%exist = .TRUE.
        molecule_list(alive(is),is)%cfc_lambda = this_lambda
        
     END IF

     ! So far P_bias only includes the probability of choosing the 
     ! insertion point from the collection of trial coordinates times the 
     ! probability of choosing each dihedral from the collection of trial 
     ! dihedrals. We need to include the number of trial coordinates, kappa_ins,
     ! and the number of trial dihedrals, kappa_dih, for each dihedral.
     kappa_tot = 1

     IF (nfragments(is) /=0 ) THEN

        kappa_tot = kappa_tot * kappa_ins

        IF (kappa_rot /= 0) THEN
           kappa_tot = kappa_tot * kappa_rot
        END IF
        
        IF (kappa_dih /=0 ) THEN
           DO ifrag = 2, nfragments(is)
              kappa_tot = kappa_tot * kappa_dih
           END DO
        END IF
     
     END IF
     
     P_bias = P_bias * REAL(kappa_tot , DP)

  END IF
  fp_bias = fp_bias * P_bias
  fp_seq = fp_seq * P_seq
  f_ring = f_ring + nrg_ring_frag_tot
  enddo

  !*****************************************************************************
  ! Step 4) Calculate the change in potential energy if the molecule is deleted
  !*****************************************************************************
  !
  ! Whether the deletion will be accepted depends on the change in potential
  ! energy, delta_e. The potential energy will be computed in 5 stages:
  !   4.1) Nonbonded intermolecular energies
  !   4.2) Bonded intramolecular energies
  !   4.3) Nonbonded intramolecular energies
  !   4.4) Ewald energies
  !   4.5) Long-range energy correction
  ! 

  do is = n1, n2, dn
  ! Recompute the COM
  CALL Get_COM(alive(is),is)

  ! Compute the distance of the atom farthest from COM
  CALL Compute_Max_COM_Distance(alive(is),is)

  ! 4.1) Nonbonded intermolecular energies

  IF (l_pair_nrg) THEN
     CALL Store_Molecule_Pair_Interaction_Arrays(alive(is),is,this_box, &
             E_inter_vdw,E_inter_qq)
  ELSE
     CALL Compute_Molecule_Nonbond_Inter_Energy(alive(is),is, &
             E_inter_vdw,E_inter_qq,inter_overlap(is))
  END IF
  f_inter_vdw = f_inter_vdw + E_inter_vdw
  f_inter_qq = f_inter_qq + E_inter_qq
  if(l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)
  enddo

  if (n1 /= n2) then
     CALL Compute_Molecule_Pair_Interaction(alive(n1),n1,alive(n2),n2,this_box,pair_vdw,pair_qq,poverlap)
  endif

  f_inter_vdw = f_inter_vdw - pair_vdw
  f_inter_qq = f_inter_qq - pair_qq

  delta_e = - f_inter_vdw - f_inter_qq

  ! 4.2) Bonded intramolecular energies

  do is = n1, n2, dn
  CALL Compute_Molecule_Bond_Energy(alive(is),is,E_bond)
  CALL Compute_Molecule_Angle_Energy(alive(is),is,E_angle)
  CALL Compute_Molecule_Dihedral_Energy(alive(is),is,E_dihedral)
  CALL Compute_Molecule_Improper_Energy(alive(is),is,E_improper)
  f_bond = f_bond + E_bond
  f_angle = f_angle + E_angle
  f_dihedral = f_dihedral + E_dihedral
  f_improper = f_improper + E_improper
  enddo

  delta_e = delta_e - f_bond - f_angle - f_dihedral - f_improper  
  
  ! 4.3) Nonbonded intramolecular energies

  do is = n1, n2, dn
  CALL Compute_Molecule_Nonbond_Intra_Energy(alive(is),is, &
          E_intra_vdw,E_intra_qq,intra_overlap(is))
  f_intra_vdw = f_intra_vdw + E_intra_vdw
  f_intra_qq = f_intra_qq + E_intra_qq
  enddo

  delta_e = delta_e - f_intra_vdw - f_intra_qq

  ! 4.4) Ewald energies

  do is = n1, n2, dn
  IF ( (int_charge_sum_style(this_box) == charge_ewald) .AND. &
       (has_charge(is)) ) THEN
    if (n1 /= n2) then
       if (is == n1) store_sum = .TRUE.
       if (is == n2) store_sum = .FALSE.
     CALL Ins_Pairs_Ewald_Reciprocal_Energy_Difference(alive(is),alive(is),is,this_box, &
             int_deletion,E_reciprocal_move)
    else
     CALL Compute_Ewald_Reciprocal_Energy_Difference(alive(is),alive(is),is,this_box, &
        int_insertion,E_reciprocal_move)
    endif
     CALL Compute_Ewald_Self_Energy_Difference(is,this_box, &
             int_deletion,E_self_move)

     f_reciprocal = f_reciprocal + E_reciprocal_move
     f_self_diff = f_self_diff + E_self_move

  END IF
  enddo

     delta_e = delta_e + f_self_diff &
                       + (f_reciprocal - energy(this_box)%ewald_reciprocal)
  ! 4.5) Long-range energy correction

  IF (int_vdw_sum_style(this_box) == vdw_cut_tail) THEN

     ! subtract off beads for this species
     nbeads_out(:) = nint_beads(:,this_box)

     do is = n1, n2, dn
     DO i = 1, natoms(is)
        i_type = nonbond_list(i,is)%atom_type_number
        nint_beads(i_type,this_box) = nint_beads(i_type,this_box) - 1
     END DO
     enddo

     CALL Compute_LR_Correction(this_box,e_lrc)
     delta_e = delta_e + ( e_lrc - energy(this_box)%lrc )

  END IF  
  
  !*****************************************************************************
  ! Step 5) Accept or reject the move
  !*****************************************************************************
  !
  ! The following quantity is calculated
  !
  !                  p_m a_mn 
  !    ln_pacc = Log[--------]
  !                  p_n a_nm
  !
  ! and passed to accept_or_reject() which executes the metropolis criterion.
  ! The acceptance criterion to delete a molecule that was inserted via CBMC is
  !
  !                                                          V
  !    ln_pacc = b(dU_mn + U_frag) + b mu' + Log[-----------------------]
  !                                              P_seq P_bias N Lambda^3
  !
  !                                          b f' V
  !            = b(dU_mn + U_frag) + Log[--------------]
  !                                      P_seq P_bias N 
  !
  ! where the primes (') indicate that additional intensive terms have been
  ! absorbed into the chemical potential and fugacity, respectively.

  f_vdw_igas = 0.0_DP
  f_qq_igas = 0.0_DP
  do is = n1, n2, dn
  E_intra_vdw_igas = 0.0_DP
  E_intra_qq_igas = 0.0_DP
  IF(species_list(is)%int_insert == int_igas) THEN
     igas_flag = .TRUE.
     CALL Compute_Molecule_Nonbond_Intra_Energy(alive(is),is, &
             E_intra_vdw_igas,E_intra_qq_igas,intra_overlap(is))
     igas_flag = .FALSE. 
     isgas = .TRUE.
     !ln_pacc = beta(this_box) * (delta_e + E_bond + E_angle &
!                                         + E_dihedral + E_improper &
!                                         + E_intra_vdw_igas + E_intra_qq_igas)
  ELSEIF(species_list(is)%fragment) THEN
     isfrag = .TRUE.
        !ln_pacc = beta(this_box) * (delta_e + E_angle + nrg_ring_frag_tot)
  END IF
  f_vdw_igas = f_vdw_igas + E_intra_vdw_igas
  f_qq_igas = f_qq_igas + E_intra_qq_igas
  enddo

  if(isgas) adden = adden + f_vdw_igas + f_qq_igas
  if(isfrag) adden = adden + f_angle + f_ring

  ln_pacc = beta(this_box) * (delta_e + adden)

  is = n1
  ! P_seq and P_bias equal 1.0 unless changed by Build_Molecule
  if (n1 /= n2) then
     nplocal = MIN(nmols(n1,this_box), nmols(n2,this_box))
     ln_pacc = ln_pacc - DLOG(fp_seq * fp_bias) &
                    - 2.0_DP*DLOG(REAL(nplocal,DP)) &
                    + 2.0_DP*DLOG(box_list(this_box)%volume)
  else
     ln_pacc = ln_pacc - DLOG(fp_seq * fp_bias) &
                    - DLOG(REAL(nmols(is,this_box),DP)) &
                    + DLOG(box_list(this_box)%volume) 
  endif
 
  IF(lchempot) THEN
     ! chemical potential is input
     if (n1 /= n2) then
        dblocal = species_list(n1)%de_broglie(this_box)*&
           species_list(n2)%de_broglie(this_box)
        loc_chem_pot = pair_chem_potential(npair)
     else
        dblocal = species_list(n1)%de_broglie(this_box)
        loc_chem_pot = species_list(n1)%chem_potential
     endif
     ln_pacc = ln_pacc + beta(this_box) * loc_chem_pot &
                       - 3.0_DP*DLOG(dblocal)
  ELSE
     ! fugacity is input
     ln_pacc = ln_pacc + DLOG(species_list(is)%fugacity) &
                       + DLOG(beta(this_box))
  END IF 
 
  accept = accept_or_reject(ln_pacc)

  IF (accept) THEN
     ! Update energies
     energy(this_box)%total = energy(this_box)%total + delta_e
     energy(this_box)%intra = energy(this_box)%intra - f_bond - f_angle &
                            - f_dihedral - f_improper
     energy(this_box)%bond = energy(this_box)%bond - f_bond
     energy(this_box)%angle = energy(this_box)%angle - f_angle
     energy(this_box)%dihedral = energy(this_box)%dihedral - f_dihedral
     energy(this_box)%improper = energy(this_box)%improper - f_improper
     energy(this_box)%intra_vdw = energy(this_box)%intra_vdw - f_intra_vdw
     energy(this_box)%intra_q = energy(this_box)%intra_q - f_intra_qq
     energy(this_box)%inter_vdw = energy(this_box)%inter_vdw - f_inter_vdw
     energy(this_box)%inter_q   = energy(this_box)%inter_q - f_inter_qq

     IF ( int_charge_sum_style(this_box) == charge_ewald .AND. &
          has_charge(is)) THEN
        energy(this_box)%ewald_reciprocal = f_reciprocal
        energy(this_box)%ewald_self = energy(this_box)%ewald_self + f_self_diff
     END IF

     IF ( int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
        energy(this_box)%lrc = e_lrc
     END IF

     ! obtain the original position of the deleted molecule so that the
     ! linked list can be updated

     do is = n1, n2, dn
     CALL Get_Position_Molecule(this_box,is,im(is),position)

     IF (position < SUM(nmols(is,:))) THEN
        DO k = position + 1, SUM(nmols(is,:))
           locate(k-1,is) = locate(k,is)
        END DO
     END IF
     
     ! move the deleted molecule to the end of alive(is) molecules
     locate(nmols(is,this_box),is) = alive(is)
     molecule_list(alive(is),is)%live = .FALSE.
     atom_list(:,alive(is),is)%exist = .FALSE.
     
     ! update the number of molecules
     nmols(is,this_box) = nmols(is,this_box) - 1

     ! Increment counter
     nsuccess(is,this_box)%deletion = nsuccess(is,this_box)%deletion + 1
     enddo

  ELSE

     is = n1

     IF ( (int_charge_sum_style(this_box) == charge_ewald) .AND. &
           (has_charge(is)) ) THEN
        ! Restore cos_sum and sin_sum. Note that these were changed when
        ! difference in reciprocal energies was computed
        cos_sum(:,this_box) = cos_sum_old(:,this_box)
        sin_sum(:,this_box) = sin_sum_old(:,this_box)
     END IF

     IF ( int_vdw_sum_style(this_box) == vdw_cut_tail ) THEN
        ! Restore the total number of bead types
        nint_beads(:,this_box) = nbeads_out(:)
     END IF

  END IF

!  IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)

END SUBROUTINE Deletion
