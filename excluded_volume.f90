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
!********************************************************************************

MODULE Excluded_Volume

  !******************************************************************************
  ! 
  ! - Andrew P. Santos
  !
  !*******************************************************************************

  USE Run_Variables
  USE Random_Generators
  USE IO_Utilities
  USE File_Names
  USE Energy_Routines
  USE Cluster_Routines
  USE Simulation_Properties

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Calculate_Excluded_Volume(this_box)

    !*********************************************************************************
    !
    ! Calculate the excluded volume
    !
    ! Step 1) Temporarily remove molecules that are not in the clusters
    ! Loop over these steps for exvol%n_iter
    !   Step 2) Choose insertion monomer and insertion location
    !   Step 3) See if the test monomer is part of a cluster
    !   Step 4) Calculate the change in potential energy, essentially the same as in insertion
    !   Step 5) Accept or reject the test insertion
    !
    ! 2/19/15  : Andrew P. Santos
    !*********************************************************************************

    INTEGER, INTENT(IN) :: this_box
    INTEGER :: i_mono_conf, i_type, iatom, im, is, i, j
    INTEGER :: ierr, line_nbr, nbr_entries
    INTEGER :: alive(2)               ! molecule indices

    INTEGER :: ispec, imol, is_clus, t_im
    INTEGER :: i_ins

    REAL(DP) :: ln_pacc, P_seq, P_bias, this_lambda, delta_e

    LOGICAL :: accept, accept_or_reject, isfrag, isgas, rej_pair, cbmc_rej_pair
    LOGICAL :: in_cluster, overlap

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: temp_locate
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: temp_live
    LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: temp_exist
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: temp_nmols

    ALLOCATE( temp_locate(MAXVAL(nmolecules),nspecies) )
    ALLOCATE( temp_live(MAXVAL(nmolecules),nspecies) )
    ALLOCATE( temp_exist(MAXVAL(natoms),MAXVAL(nmolecules)+1,nspecies) )
    ALLOCATE( temp_nmols(nspecies,nbr_boxes) )
    ln_pacc = 0.0_DP
    P_seq = 1.0_DP
    P_bias = 1.0_DP

    exvol%excluded = 0
    exvol%criteria = LOG(1E8)
    exvol%trials = exvol%trials + 1

    !*********************************************************************************
    !   Step 1) Temporarily remove molecules that are not in the clusters
    !*********************************************************************************

    ! save the molecule position arrays etc.
    temp_locate = locate
    temp_live = molecule_list%live
    temp_exist = atom_list%exist
    temp_nmols = nmols

    CALL Remove_Small_Clusters(this_box)

    ins_loop: DO i_ins = 1, exvol%n_iter
       !*********************************************************************************
       !   Step 2) Choose insertion monomer and insertion location
       !*********************************************************************************
       is = exvol%species
       im = 1
       locate(im,is) = im
   
       CALL Select_Monomer(im, is, this_box)
   
       !*********************************************************************************
       !   Step 3) See if the test monomer is part of a cluster
       !*********************************************************************************
   
       DO ispec = 1, cluster%n_species_type
           DO imol = 1, nmols(ispec,this_box)
   
               in_cluster = Neighbor(imol, im, cluster%species_type(ispec), is)
               IF (in_cluster) THEN
                   ! remove monomer
                   CALL Select_Monomer(im, is, this_box)
                   exvol%excluded = exvol%excluded + 1
                   CYCLE ins_loop
               END IF
   
           END DO
       END DO
   
       !*********************************************************************************
       !   Step 4) Calculate the change in potential energy, essentially the same as in insertion
       !*********************************************************************************
   
       CALL Get_Monomer_DeltaE(im, is, this_box, overlap, delta_e)
   
       IF (overlap) THEN
           CALL Select_Monomer(im, is, this_box)
           exvol%excluded = exvol%excluded + 1
           CYCLE ins_loop
       END IF
   
   
       !*********************************************************************************
       !   Step 5) Accept or reject the test insertion
       !*********************************************************************************
   
       ln_pacc = beta(this_box) * delta_e 
   
       ! P_seq and P_bias equal 1.0 unless changed by Build_Molecule.
       P_seq = 1.0_DP
       P_bias = 1.0_DP
       ln_pacc = ln_pacc + DLOG(P_seq * P_bias) !&
                         !+ 3.0_DP * DLOG(species_list(is)%de_broglie(this_box))
                         !+ 2.0_DP*DLOG(REAL(nmols(is,this_box)+1,DP)) &
                         !- 2.0_DP*DLOG(box_list(this_box)%volume) 
                         !- species_list(is)%chem_potential * beta(this_box) 
   
       IF (ln_pacc > exvol%criteria) THEN
          exvol%excluded = exvol%excluded + 1
       END IF
       
       ! Remove inserting monomer
       CALL Select_Monomer(im, is, this_box)
   
    END DO ins_loop
    ! restore molecules that were in clusters smaller than the minimum size
    locate = temp_locate 
    molecule_list%live = temp_live 
    atom_list%exist = temp_exist 
    nmols = temp_nmols

  END SUBROUTINE Calculate_Excluded_Volume

  SUBROUTINE Select_Monomer(im, is, this_box)
    INTEGER, INTENT(IN) :: im, is, this_box
    INTEGER :: this_atom, this_fragment, i
    INTEGER :: alive(2)               ! molecule indices

    ! local fragment declarations
    INTEGER :: total_frags  ! total number of conformations for this fragment in the library
    INTEGER :: frag_start   ! random fragment to start growing from
    INTEGER :: frag_total   ! number of non-zero entries in frag_order
    INTEGER :: frag_type
    REAL(DP) :: this_lambda

    alive(is) = locate(im,is)
    IF (molecule_list(im, is)%live == .FALSE.) THEN
    this_lambda = 1.0_DP
       frag_start = 1
       frag_type = frag_list(frag_start,is)%type
       ! Pull from the fragment reservoir with uniform probability
       total_frags = frag_list(frag_start,is)%nconfig
       this_fragment = INT(rranf() * total_frags) + 1
   
       ! Read the coordinates for every atom
       DO i = 1, frag_list(frag_start,is)%natoms 
           this_atom = frag_list(frag_start,is)%atoms(i)
           
           atom_list(this_atom,im,is)%rxp = frag_coords(i,this_fragment,frag_type)%rxp
           atom_list(this_atom,im,is)%ryp = frag_coords(i,this_fragment,frag_type)%ryp
           atom_list(this_atom,im,is)%rzp = frag_coords(i,this_fragment,frag_type)%rzp
   
           atom_list(this_atom,im,is)%exist = .TRUE.
       END DO
   
       ! Turn on the molecule and its individual atoms
       molecule_list(alive(is),is)%which_box = this_box
       molecule_list(alive(is),is)%cfc_lambda = this_lambda
       molecule_list(alive(is),is)%molecule_type = int_normal
     
       molecule_list(alive(is),is)%live = .TRUE.
   
       ! Randomize the molecule's COM position anywhere in the box.
       IF ( box_list(this_box)%int_box_shape == int_cubic ) THEN
           molecule_list(alive(is),is)%xcom = &
                                (rranf() - 0.5_DP) * box_list(this_box)%length(1,1)
           molecule_list(alive(is),is)%ycom = &
                                (rranf() - 0.5_DP) * box_list(this_box)%length(2,2)
           molecule_list(alive(is),is)%zcom = &
                                (rranf() - 0.5_DP) * box_list(this_box)%length(3,3)
       END IF
        
       atom_list(:,alive(is),is)%rxp = atom_list(:,alive(is),is)%rxp + molecule_list(alive(is),is)%xcom
       atom_list(:,alive(is),is)%ryp = atom_list(:,alive(is),is)%ryp + molecule_list(alive(is),is)%ycom
       atom_list(:,alive(is),is)%rzp = atom_list(:,alive(is),is)%rzp + molecule_list(alive(is),is)%zcom
   
       ! Molecule COM may be outside the box boundary if grown via CBMC, so wrap
       ! the molecule coordinates back in the box (if needed)
       CALL Fold_Molecule(alive(is),is,this_box)
   
       ! Recompute the COM in case the molecule was wrapped
       CALL Get_COM(alive(is),is)
   
       ! Compute the distance of the atom farthest from COM
       CALL Compute_Max_COM_Distance(alive(is),is)

    !Remove the Monomer
    ELSE IF (molecule_list(im, is)%live == .TRUE.) THEN
       ! clean up
       molecule_list(alive(is),is)%live = .FALSE.
       atom_list(:,alive(is),is)%exist = .FALSE.
       molecule_list(alive(is),is)%molecule_type = int_none

    END IF

  END SUBROUTINE Select_Monomer


  SUBROUTINE Remove_Small_Clusters(this_box)
    INTEGER, INTENT(IN) :: this_box
    INTEGER :: n_removed, ispec, imol, is_clus, t_im, position
    INTEGER :: k
    n_removed = 0
    DO ispec = 1, cluster%n_species_type
        is_clus = cluster%species_type(ispec)
        MoleculeLoop: DO imol = 1, nmolecules(is_clus)
            ! Make sure that the molecule exists in the simulation
            t_im = locate(imol-n_removed, is_clus)
            IF( .NOT. molecule_list(t_im,is_clus)%live ) CYCLE MoleculeLoop

            ! IF the molecule is not in a cluster of the specified size
            IF (cluster%N( cluster%clabel(t_im, is_clus) ) < cluster%M_olig(is_clus)) THEN

                CALL Get_Position_Molecule(this_box,is_clus,t_im,position)
                DO k = imol + 1 - n_removed, SUM(nmols(is_clus,:))
                    locate(k-1,is_clus) = locate(k,is_clus)
                END DO
     
                ! move the deleted molecule to the end of alive(is_clus) molecules
                locate(nmols(is_clus,this_box),is_clus) = t_im
                molecule_list(t_im,is_clus)%live = .FALSE.
                atom_list(:,t_im,is_clus)%exist = .FALSE.
     
                ! update the number of molecules
                nmols(is_clus,this_box) = nmols(is_clus,this_box) - 1
                n_removed = n_removed + 1
            END IF
        END DO MoleculeLoop
    END DO


  END SUBROUTINE Remove_Small_Clusters

  SUBROUTINE Get_Monomer_DeltaE(im, is, this_box, overlap, delta_e)
    INTEGER, INTENT(IN) :: im, is, this_box
    LOGICAL, INTENT(OUT) :: overlap
    REAL(DP), INTENT(OUT) :: delta_e
    REAL(DP) :: E_bond, E_angle, E_dihedral, E_improper
    REAL(DP) :: E_intra_vdw, E_intra_qq
    REAL(DP) :: E_inter_vdw, E_inter_qq
    REAL(DP) :: pair_vdw, pair_qq, igas_en
    REAL(DP) :: E_reciprocal_move, E_self_move, E_lrc
    LOGICAL :: inter_overlap(2), intra_overlap(2)
    INTEGER :: alive(2), i, i_type               ! molecule indices

    overlap = .FALSE.
    inter_overlap(:) = .FALSE.
    intra_overlap(:) = .FALSE.
    E_inter_vdw = 0.0_DP
    E_inter_qq = 0.0_DP
    E_intra_vdw = 0.0_DP
    E_intra_qq = 0.0_DP

    alive(is) = locate(im,is)

    ! Calculate the potential energy interaction between the inserted molecule
    ! and the rest of the system
    CALL Compute_Molecule_Nonbond_Inter_Energy(alive(is),is, &
            E_inter_vdw,E_inter_qq,inter_overlap(is))

    ! Calculate the nonbonded energy interaction within the inserted molecule
    CALL Compute_Molecule_Nonbond_Intra_Energy(alive(is),is, &
            E_intra_vdw,E_intra_qq,intra_overlap(is))

    IF (inter_overlap(is) .OR. intra_overlap(is)) THEN
        overlap = .TRUE.
        RETURN
    END IF

    ! There are no overlaps, so we can calculate the change in potential energy.
    ! Already have the change in nonbonded energies
    delta_e = E_inter_vdw + E_inter_qq + E_intra_vdw + E_intra_qq

    E_bond = 0.0_DP
    E_angle = 0.0_DP
    E_dihedral = 0.0_DP
    E_improper = 0.0_DP
    CALL Compute_Molecule_Bond_Energy(alive(is),is,E_bond)
    CALL Compute_Molecule_Angle_Energy(alive(is),is,E_angle)
    CALL Compute_Molecule_Dihedral_Energy(alive(is),is,E_dihedral)
    CALL Compute_Molecule_Improper_Energy(alive(is),is,E_improper)

    delta_e = delta_e + E_bond + E_angle + E_dihedral + E_improper

    E_reciprocal_move = 0.0_DP
    E_self_move = 0.0_DP
    IF ( (int_charge_sum_style(this_box) == charge_ewald) .AND. &
         (has_charge(is)) ) THEN
       CALL Compute_Ewald_Reciprocal_Energy_Difference(alive(is),alive(is),is,this_box, &
             int_insertion,E_reciprocal_move)

       CALL Compute_Ewald_Self_Energy_Difference(alive(is),is,this_box, &
             int_insertion,E_self_move)
    END IF

    delta_e = delta_e + E_self_move + E_reciprocal_move - energy(this_box)%ewald_reciprocal

    ! 3.6) Long-range energy correction
    IF (int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
       ! increase number of integer beads
       nbeads_in = nint_beads(:,this_box)

       DO i = 1, natoms(is)
          i_type = nonbond_list(i,is)%atom_type_number
          nint_beads(i_type,this_box) = nint_beads(i_type,this_box) + 1
       END DO

       CALL Compute_LR_correction(this_box,e_lrc)
       delta_e = delta_e + e_lrc - energy(this_box)%lrc
    END IF

    ! Restore the sums made before
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

  END SUBROUTINE Get_Monomer_DeltaE

END MODULE Excluded_Volume
