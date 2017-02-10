!*******************************************************************************
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

SUBROUTINE Insertion(this_box,mcstep,randno)

  !*****************************************************************************
  ! 
  ! PURPOSE: attempt to insert a molecule via configurational bias monte carlo
  !
  ! Called by
  !
  !    gcmc_driver
  !
  ! Revision history
  !
  !   12/10/13  : Beta version created
  !   Version 1.1
  !     04/21/15  Corrected acceptance criteria
  !     05/01/15  Documented this code
  !   
  ! DESCRIPTION: This subroutine performs the following steps:
  !
  ! Step 1) Select a species with uniform probability
  ! Step 2) Choose a position, orientation and conformation for the 
  !         to-be-inserted molecule
  ! Step 3) Calculate the change in potential energy if the molecule is inserted
  ! Step 4) Accept or reject the move
  !
  !*****************************************************************************

  USE Run_Variables
  USE Energy_Routines
  USE IO_Utilities
  USE Random_Generators
  USE Rotation_Routines
  USE Fragment_Growth

  !*****************************************************************************
  ! Declare and Initialize Variables
  !*****************************************************************************

  IMPLICIT NONE

  ! Arguments
  INTEGER :: this_box ! attempt to insert a molecule in this_box
  INTEGER :: mcstep   ! not used
  REAL(DP) :: randno  ! not used

  ! Local declarations
  INTEGER :: i, i_type               ! atom indices
  INTEGER :: ifrag                   ! fragment indices
  INTEGER :: im, alive(2)               ! molecule indices
  INTEGER :: is, is_rand, is_counter ! species indices
  INTEGER :: kappa_tot, which_anchor
  INTEGER, ALLOCATABLE :: frag_order(:)
  INTEGER :: rand_igas, tot_mols
  INTEGER :: tn1, tn2, n1, n2, nplocal, npair, dn

  REAL(DP) :: ppt, pp(n_insertable), randnpair, loc_chem_pot
  REAL(DP) :: dx, dy, dz, delta_e
  REAL(DP) :: E_bond, E_angle, E_dihedral, E_improper
  REAL(DP) :: f_bond, f_angle, f_dihedral, f_improper
  REAL(DP) :: E_intra_vdw, E_intra_qq
  REAL(DP) :: E_inter_vdw, E_inter_qq
  REAL(DP) :: f_inter_vdw, f_inter_qq, f_intra_vdw, f_intra_qq
  REAL(DP) :: pair_vdw, pair_qq, igas_en
  REAL(DP) :: E_reciprocal_move, E_self_move, E_lrc
  REAL(DP) :: f_reciprocal, f_self_diff, suben
  REAL(DP) :: nrg_ring_frag_tot, f_ring, dblocal
  REAL(DP) :: ln_pacc, P_seq, P_bias, this_lambda
  REAL(DP) :: fp_bias, fp_seq

  LOGICAL :: inter_overlap(2), cbmc_overlap(2), intra_overlap(2), poverlap
  LOGICAL :: accept, accept_or_reject, isfrag, isgas, rej_pair, cbmc_rej_pair

  ! Initialize variables
  ln_pacc = 0.0_DP
  P_seq = 1.0_DP
  P_bias = 1.0_DP
  this_lambda = 1.0_DP
  nrg_ring_frag_tot = 0.0_DP
  inter_overlap(:) = .FALSE.
  cbmc_overlap(:) = .FALSE.
  intra_overlap(:) = .FALSE.
  poverlap = .FALSE.
  ppt = 0.0_DP
  pp(:) = 0.0_DP
  tn1 = 1
  n1 = 1
  tn2 = 2
  n2 = 2
  dn = 1
  f_bond = 0.0_DP
  f_angle = 0.0_DP
  f_dihedral = 0.0_DP
  f_improper = 0.0_DP
  f_inter_vdw = 0.0_DP
  f_inter_qq = 0.0_DP
  f_intra_vdw = 0.0_DP
  f_intra_qq = 0.0_DP
  f_reciprocal = 0.0_DP
  f_self_diff = 0.0_DP
  igas_en = 0.0_DP
  pair_qq = 0.0_DP
  pair_vdw = 0.0_DP
  isfrag = .FALSE.
  isgas = .FALSE.
  rej_pair = .FALSE.
  cbmc_rej_pair = .FALSE.
  suben = 0.0_DP
  f_ring = 0.0_DP
  fp_bias = 1.0_DP
  fp_seq = 1.0_DP
  delta_e = 0.0_DP
  dblocal = 0.0_DP
  npair = 1

  !*****************************************************************************
  ! Step 1) Randomly select a species
  !*****************************************************************************
  !
  ! All species may not be insertable. For example, in a simulation of dilute
  ! water (species 3) and CO2 (species 4) in an ionic liquid (species 1 and 2), 
  ! the number of ionic liquid molecules may be fixed and only the numbers of
  ! water and CO2 allowed to fluctuate. First, choose a random integer between 1
  ! and the number of insertable species, nspec_insert:

  if (any(species_list(:)%pair_insert) .eqv. .TRUE.) then
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

  dn = n2 - n1
  if (dn .EQ. 0) dn = 1

  do is = n1, n2, dn
  tot_mols = SUM(nmols(is,:)) ! summed over the number of boxes?

  ! Check that tot_mols is less than the maximum allowable, nmolecules(is)

  IF (tot_mols == nmolecules(is)) THEN
     err_msg = ""
     err_msg(1) = 'Number of molecule exceeds limit of ' // &
                  INT_to_String(tot_mols) 
     err_msg(2) = 'Increase molecule number limit in input file '
     CALL Clean_Abort(err_msg,'Insertion')
     ! exit if we are attempting an insertion above the maximum allowable
  END IF
  enddo

  ! Now that an insertion will be attempted, we need to do some bookkeeping:
  !  * Increment the counters to compute success ratios

  do is = n1, n2, dn
  ntrials(is,this_box)%insertion = ntrials(is,this_box)%insertion + 1

  !  * Assign a locate number for this molecule

  IF ( locate(tot_mols+1,is) == 0 ) THEN
     locate(tot_mols+1,is) = tot_mols + 1
     ! otherwise we will use the locate number of a previously deleted molecule 
     ! that has been moved to the end of the array.
  END IF

  !  * Set properties of the to-be-inserted molecule

  alive(is) = locate(tot_mols+1,is)
  molecule_list(alive(is),is)%which_box = this_box
  molecule_list(alive(is),is)%cfc_lambda = this_lambda
  molecule_list(alive(is),is)%molecule_type = int_normal
  enddo

  tot_trials(this_box) = tot_trials(this_box) + 1

  ! With the bookkeeping completed, we are ready to attempt the insertion
  
  !*****************************************************************************
  ! Step 2) Choose a position, orientation and conformation for the 
  !         to-be-inserted molecule
  !*****************************************************************************
  !
  ! If the molecule:
  !   * has a fragment library and is not an ideal gas, then the conformation
  !     will be grown fragment-by-fragment using CBMC. The resulting 
  !     conformation will be molded to a high probability position.
  !   * has no fragment library, then there is only one conformation. Position 
  !     and orientation are random.
  !   * is an ideal gas, then molecular conformations (not fragment 
  !     conformations?) were sampled according to their Boltzmann weight. One 
  !     is chosen at random. Position and orientation are random.
 
  do is = n1, n2, dn 
  P_seq = 1.0_DP
  P_bias = 1.0_DP
  nrg_ring_frag_tot = 0.0_DP
  IF (species_list(is)%fragment .AND. &
     (species_list(is)%int_insert .NE. int_igas) ) THEN

     ! Build_Molecule places the first fragment, then calls Fragment_Placement
     ! to place the additional fragments 
     del_flag = .FALSE.     ! Change the coordinates of 'alive(is)'
     get_fragorder = .TRUE.
     ALLOCATE(frag_order(nfragments(is)))
     CALL Build_Molecule(alive(is),is,this_box,frag_order,this_lambda, &
             P_seq,P_bias,nrg_ring_frag_tot,cbmc_overlap(is))
     DEALLOCATE(frag_order)

     ! Turn the molecule on
     molecule_list(alive(is),is)%live = .TRUE.
     atom_list(:,alive(is),is)%exist = .TRUE.

     ! So far P_bias only includes the probability of choosing the 
     ! insertion point from the collection of trial coordinates times the 
     ! probability of choosing each dihedral from the collection of trial 
     ! dihedrals. We need to include the number of trial coordinates, kappa_ins,
     ! and the number of trial dihedrals, kappa_dih, for each dihedral.
     kappa_tot = 1

     IF (nfragments(is) /= 0 ) THEN

        kappa_tot = kappa_tot * kappa_ins

        IF (kappa_rot /= 0 ) THEN
           kappa_tot = kappa_tot * kappa_rot
        END IF

        IF (kappa_dih /= 0 ) THEN
           DO ifrag = 1, nfragments(is) - 1
              kappa_tot = kappa_tot * kappa_dih
           END DO
        END IF

     END IF

     P_bias = P_bias * REAL(kappa_tot, DP)

  ELSE
 
     ! The molecule does not have a fragment libary or is an ideal gas. The 
     ! position and orientation will be randomized independent of the molecular 
     ! conformation.

     ! Turn the molecule on
     atom_list(:,alive(is),is)%exist = .TRUE.
     molecule_list(alive(is),is)%live = .TRUE.

     ! Now we need to grab the molecule's conformation

     IF(species_list(is)%int_insert == int_random) THEN
     
        ! A rigid molecule has only one possible conformation

        molecule_list(alive(is),is)%xcom = species_list(is)%xcom
        molecule_list(alive(is),is)%ycom = species_list(is)%ycom
        molecule_list(alive(is),is)%zcom = species_list(is)%zcom
        
        atom_list(:,alive(is),is)%rxp = init_list(:,1,is)%rxp
        atom_list(:,alive(is),is)%ryp = init_list(:,1,is)%ryp
        atom_list(:,alive(is),is)%rzp = init_list(:,1,is)%rzp


     ELSE IF(species_list(is)%int_insert == int_igas) THEN

        ! An ideal gas molecule can adopt any conformation independent of its 
        ! local environment. The conformations are sampled according to their 
        ! Boltzmann weight. Read one in at random:

        rand_igas = (rranf() * n_igas(is)) + 1

        molecule_list(alive(is),is)%xcom = molecule_list_igas(rand_igas,is)%xcom
        molecule_list(alive(is),is)%ycom = molecule_list_igas(rand_igas,is)%ycom
        molecule_list(alive(is),is)%zcom = molecule_list_igas(rand_igas,is)%zcom

        atom_list(:,alive(is),is)%rxp = atom_list_igas(:,rand_igas,is)%rxp
        atom_list(:,alive(is),is)%ryp = atom_list_igas(:,rand_igas,is)%ryp
        atom_list(:,alive(is),is)%rzp = atom_list_igas(:,rand_igas,is)%rzp

     END IF   

     ! Randomize the molecule's orientation.
     
     CALL Rotate_Molecule_Eulerian(alive(is),is)
     
     ! Randomize the molecule's COM position anywhere in the box.

     IF ( box_list(this_box)%int_box_shape == int_cubic ) THEN
        
        molecule_list(alive(is),is)%xcom = &
                             (rranf() - 0.5_DP) * box_list(this_box)%length(1,1)
        molecule_list(alive(is),is)%ycom = &
                             (rranf() - 0.5_DP) * box_list(this_box)%length(2,2)
        molecule_list(alive(is),is)%zcom = &
                             (rranf() - 0.5_DP) * box_list(this_box)%length(3,3)
        
     END IF
     
     ! Move the molecule to the chosen COM position

     IF(species_list(is)%int_insert == int_random) THEN
     
        dx = molecule_list(alive(is),is)%xcom - species_list(is)%xcom
        dy = molecule_list(alive(is),is)%ycom - species_list(is)%ycom
        dz = molecule_list(alive(is),is)%zcom - species_list(is)%zcom

     ELSE

        dx = molecule_list(alive(is),is)%xcom &
           - molecule_list_igas(rand_igas,is)%xcom
        dy = molecule_list(alive(is),is)%ycom &
           - molecule_list_igas(rand_igas,is)%ycom
        dz = molecule_list(alive(is),is)%zcom &
           - molecule_list_igas(rand_igas,is)%zcom

     END IF        
     
     atom_list(:,alive(is),is)%rxp = atom_list(:,alive(is),is)%rxp + dx
     atom_list(:,alive(is),is)%ryp = atom_list(:,alive(is),is)%ryp + dy
     atom_list(:,alive(is),is)%rzp = atom_list(:,alive(is),is)%rzp + dz

  END IF
  fp_bias = fp_bias * P_bias
  fp_seq = fp_seq * P_seq
  f_ring = f_ring + nrg_ring_frag_tot
  enddo

  !*****************************************************************************
  ! Step 3) Calculate the change in potential energy if the molecule is inserted
  !*****************************************************************************
  !
  ! Whether the insertion will be accepted depends on the change in potential
  ! energy, delta_e. The potential energy will be computed in 5 stages:
  !   3.1) Nonbonded energies
  !   3.2) Reject the move if there is any core overlap
  !   3.3) Bonded intramolecular energies
  !   3.4) Ewald energies
  !   3.5) Long-range energy correction
  ! 
  ! 3.1) Nonbonded energies
  ! If the inserted molecule overlaps cores with any other molecule in  
  ! the system, the move will be rejected. If the inserted molecule was grown 
  ! via CBMC and all attempted configurations overlapped cores, the cbmc_overlap
  ! flag equals .TRUE. If the molecule was inserted randomly, we still need to 
  ! detect core overlaps.

!FSL to speed up the calculation, we're skipping the next loop if there's any overlap
  do is = n1, n2, dn
       if (cbmc_overlap(is)) cbmc_rej_pair = .TRUE.
  enddo


  if (.NOT. cbmc_rej_pair) then
  do is = n1, n2, dn
  E_inter_vdw = 0.0_DP
  E_inter_qq = 0.0_DP
  E_intra_vdw = 0.0_DP
  E_intra_qq = 0.0_DP

    ! Molecule COM may be outside the box boundary if grown via CBMC, so wrap
    ! the molecule coordinates back in the box (if needed)
    CALL Fold_Molecule(alive(is),is,this_box)

    ! Recompute the COM in case the molecule was wrapped
    CALL Get_COM(alive(is),is)

    ! Compute the distance of the atom farthest from COM
    CALL Compute_Max_COM_Distance(alive(is),is)

    ! Calculate the potential energy interaction between the inserted molecule
    ! and the rest of the system
    CALL Compute_Molecule_Nonbond_Inter_Energy(alive(is),is, &
            E_inter_vdw,E_inter_qq,inter_overlap(is))

    ! Calculate the nonbonded energy interaction within the inserted molecule
    CALL Compute_Molecule_Nonbond_Intra_Energy(alive(is),is, &
            E_intra_vdw,E_intra_qq,intra_overlap(is))

         f_inter_vdw = f_inter_vdw + E_inter_vdw
         f_inter_qq = f_inter_qq + E_inter_qq
         f_intra_vdw = f_intra_vdw + E_intra_vdw
         f_intra_qq = f_intra_qq + E_intra_qq         

  enddo 
  endif

  ! 3.3) Reject the move if there is any core overlap
  do is = n1, n2, dn
  IF (cbmc_overlap(is) .OR. inter_overlap(is) .OR. intra_overlap(is)) THEN
        rej_pair = .TRUE.
  ENDIF
  enddo

  if (rej_pair) then
  do is = n1, n2, dn
     molecule_list(alive(is),is)%live = .FALSE.
     atom_list(:,alive(is),is)%exist = .FALSE.
     if (is == n2) RETURN
  enddo
  endif

!FSL this calculation should only be executed if n1 /= n2, updating "is" values
  if (n1 /= n2) then
     CALL Compute_Molecule_Pair_Interaction(alive(n1),n1,alive(n2),n2,this_box,pair_vdw,pair_qq,poverlap)
  endif

  f_inter_vdw = f_inter_vdw - pair_vdw
  f_inter_qq = f_inter_qq - pair_qq

  ! There are no overlaps, so we can calculate the change in potential energy.
  !
  ! Already have the change in nonbonded energies

  delta_e = f_inter_vdw + f_inter_qq 
  delta_e = delta_e + f_intra_vdw + f_intra_qq

  ! 3.4) Bonded intramolecular energies
  ! If the molecule was grown via CBMC, we already have the intramolecular 
  ! bond energies? Otherwise we need to compute them.

  do is = n1, n2, dn
  E_bond = 0.0_DP
  E_angle = 0.0_DP
  E_dihedral = 0.0_DP
  E_improper = 0.0_DP
  IF(species_list(is)%int_insert == int_random) THEN

     CALL Compute_Molecule_Bond_Energy(alive(is),is,E_bond)
     CALL Compute_Molecule_Angle_Energy(alive(is),is,E_angle)
     CALL Compute_Molecule_Dihedral_Energy(alive(is),is,E_dihedral)
     CALL Compute_Molecule_Improper_Energy(alive(is),is,E_improper)

  ELSE IF (species_list(is)%int_insert == int_igas) THEN

     E_bond = energy_igas(rand_igas,is)%bond
     E_angle = energy_igas(rand_igas,is)%angle
     E_dihedral = energy_igas(rand_igas,is)%dihedral
     E_improper = energy_igas(rand_igas,is)%improper

  END IF

     f_bond = f_bond + E_bond
     f_angle = f_angle + E_angle
     f_dihedral = f_dihedral + E_dihedral
     f_improper = f_improper + E_improper

  enddo

  delta_e = delta_e + f_bond + f_angle + f_dihedral + f_improper

  ! 3.5) Ewald energies

  do is = n1, n2, dn
  E_reciprocal_move = 0.0_DP
  E_self_move = 0.0_DP
  IF ( (int_charge_sum_style(this_box) == charge_ewald) .AND. &
       (has_charge(is)) ) THEN
    if (n1 /= n2) then
     if (is == n1) store_sum = .TRUE.
     if (is == n2) store_sum = .FALSE.
     CALL Ins_Pairs_Ewald_Reciprocal_Energy_Difference(alive(is),alive(is),is,this_box, &
             int_insertion,E_reciprocal_move)
    else
     CALL Compute_Ewald_Reciprocal_Energy_Difference(alive(is),alive(is),is,this_box, &
             int_insertion,E_reciprocal_move)
    endif

     CALL Compute_Ewald_Self_Energy_Difference(is,this_box, &
             int_insertion,E_self_move)

     f_reciprocal = f_reciprocal + E_reciprocal_move
     f_self_diff = f_self_diff + E_self_move

  END IF
  enddo

     delta_e = delta_e + f_self_diff &
                       + f_reciprocal - energy(this_box)%ewald_reciprocal
  ! 3.6) Long-range energy correction

  IF (int_vdw_sum_style(this_box) == vdw_cut_tail) THEN

     ! increase number of integer beads
     nbeads_in = nint_beads(:,this_box)

     do is = n1, n2, dn
     DO i = 1, natoms(is)
        i_type = nonbond_list(i,is)%atom_type_number
        nint_beads(i_type,this_box) = nint_beads(i_type,this_box) + 1
     END DO
     enddo

     CALL Compute_LR_correction(this_box,e_lrc)
     delta_e = delta_e + e_lrc - energy(this_box)%lrc

  END IF

  !*****************************************************************************
  ! Step 4) Accept or reject the move
  !*****************************************************************************
  !
  ! The following quantity is calculated
  !
  !                  p_m a_mn 
  !    ln_pacc = Log[--------]
  !                  p_n a_nm
  !
  ! and passed to accept_or_reject() which executes the metropolis criterion.
  ! The acceptance criterion to insert a molecule via CBMC is
  !
  !                                            P_seq P_bias (N + 1) Lambda^3
  !    ln_pacc = b(dU_mn-U_frag) - b mu' + Log[-----------------------------]  
  !                                                         V
  !
  !                                    P_seq P_bias (N + 1) 
  !            = b(dU_mn-U_frag) + Log[--------------------]
  !                                           b f' V
  !
  ! where the primes (') indicate that additional intensive terms have been
  ! absorbed into the chemical potential and fugacity, respectively.

  ! Compute the acceptance criterion
  do is = n1, n2, dn
  IF(species_list(is)%int_insert == int_igas) THEN 
     isgas = .TRUE.
     igas_en = igas_en + energy_igas(rand_igas,is)%total
  ELSEIF (species_list(is)%fragment) THEN
     isfrag = .TRUE.
  END IF
  enddo

  if(isgas) suben = suben + igas_en
  if(isfrag) suben = suben + E_angle + f_ring

  ln_pacc = beta(this_box) * (delta_e - suben)

  is = n1
  ! P_seq and P_bias equal 1.0 unless changed by Build_Molecule.
  if (n1 /= n2) then
     nplocal = MIN(nmols(n1,this_box), nmols(n2,this_box))
     ln_pacc = ln_pacc + DLOG(fp_seq * fp_bias) &
                    + 2.0_DP*DLOG(REAL(nplocal+1,DP)) &
                    - 2.0_DP*DLOG(box_list(this_box)%volume) 
  else
     ln_pacc = ln_pacc + DLOG(fp_seq * fp_bias) &
                    + DLOG(REAL(nmols(is,this_box)+1,DP)) &
                    - DLOG(box_list(this_box)%volume) 
  endif

  IF(lchempot) THEN
    if (n1 /= n2) then
       dblocal = species_list(n1)%de_broglie(this_box)*& 
          species_list(n2)%de_broglie(this_box)
       loc_chem_pot = pair_chem_potential(npair)
    else
       dblocal = species_list(n1)%de_broglie(this_box)
       loc_chem_pot = species_list(n1)%chem_potential
    endif
     ! chemical potential is input
     ln_pacc = ln_pacc - loc_chem_pot * beta(this_box) &
                       + 3.0_DP*DLOG(dblocal)
  ELSE
     ! fugacity is input
     ln_pacc = ln_pacc - DLOG(species_list(n1)%fugacity) &
                       - DLOG(beta(this_box))
  END IF
  
  accept = accept_or_reject(ln_pacc)
  
  IF (accept) THEN

     ! update the number of molecules
     do is = n1, n2, dn
     nmols(is,this_box) = nmols(is,this_box) + 1
     enddo
     ! update the energies
     energy(this_box)%total = energy(this_box)%total + delta_e
     energy(this_box)%intra = energy(this_box)%intra + f_bond + f_angle &
                            + f_dihedral + f_improper
     energy(this_box)%bond = energy(this_box)%bond + f_bond
     energy(this_box)%angle = energy(this_box)%angle + E_angle
     energy(this_box)%dihedral = energy(this_box)%dihedral + f_dihedral
     energy(this_box)%improper = energy(this_box)%improper + f_improper
     energy(this_box)%intra_vdw = energy(this_box)%intra_vdw + f_intra_vdw
     energy(this_box)%intra_q = energy(this_box)%intra_q + f_intra_qq
     energy(this_box)%inter_vdw = energy(this_box)%inter_vdw + f_inter_vdw
     energy(this_box)%inter_q = energy(this_box)%inter_q + f_inter_qq

     is = n1
     IF ( int_charge_sum_style(this_box) == charge_ewald .AND. &
          has_charge(is)) THEN
        energy(this_box)%ewald_reciprocal = f_reciprocal
        energy(this_box)%ewald_self = energy(this_box)%ewald_self + f_self_diff
     END IF

     IF (int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
        energy(this_box)%lrc = e_lrc
     END IF

     ! Increment counter
     do is = n1, n2, dn
     nsuccess(is,this_box)%insertion = nsuccess(is,this_box)%insertion + 1
     enddo

  ELSE
  
     do is = n1, n2, dn
     molecule_list(alive(is),is)%live = .FALSE.
     atom_list(:,alive(is),is)%exist = .FALSE.
     molecule_list(alive(is),is)%molecule_type = int_none
     enddo
     
     is = n1
     IF ( int_charge_sum_style(this_box) == charge_ewald .AND. &
          has_charge(is) ) THEN
        ! Restore cos_sum and sin_sum. Note that these were changed when the
        ! difference in reciprocal energies was computed.
        cos_sum(:,this_box) = cos_sum_old(:,this_box)
        sin_sum(:,this_box) = sin_sum_old(:,this_box)
     END IF

     IF ( int_vdw_sum_style(this_box) == vdw_cut_tail ) THEN
        ! Restore the total number of bead types
        nint_beads(:,this_box) = nbeads_in(:)
     END IF
     
  END IF

END SUBROUTINE Insertion
