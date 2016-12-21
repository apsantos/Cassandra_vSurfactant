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

MODULE Cluster_Routines

  !******************************************************************************
  ! 
  ! This module contains routines that identify clusters using a version of the 
  ! Hoshen-Kopelman algorithm.
  !
  ! - Andrew P. Santos
  !
  !*******************************************************************************

  USE Run_Variables
  USE Random_Generators
  USE Energy_Routines
  USE Pair_Nrg_Routines

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Translate_Cluster(this_box)

    !*********************************************************************************
    !
    ! Similar to regular translation
    !
    ! CALLED BY 
    !
    !        *_driver
    !
    ! Step 1) Pick a cluster with uniform probability
    ! Step 2) Choose insertion monomer and insertion location
    ! Step 3) Reject if any monomer has become part of a cluster
    ! Step 4) Calculate the change in potential energy
    ! Step 5) Accept or reject the test insertion
    !   
    ! 2/19/15  : Andrew P. Santos
    !*********************************************************************************

    INTEGER :: this_box
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: neigh_list
    ! Arguments
  
    ! Local declarations
    INTEGER  :: ibox       ! box index
    INTEGER  :: is, iclus, iclus_N, im_clab  ! species index
    INTEGER  :: im, imol, i  ! molecule indices
    INTEGER  :: total_mols ! number of molecules in the system

    INTEGER  :: c, cs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: iclus_mol, iclus_is

    INTEGER, ALLOCATABLE, DIMENSION(:) :: old_N
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: old_clabel
    INTEGER :: old_clusmax
    
    REAL(DP) :: nmols_box(nbr_boxes)
    REAL(DP), ALLOCATABLE :: x_box(:), x_species(:)
    REAL(DP) :: rand_no
    REAL(DP) :: dx, dy, dz
    REAL(DP) :: delta_e, ln_pacc, success_ratio
    REAL(DP) :: E_vdw, E_qq, E_vdw_move, E_qq_move, E_reciprocal_move
    REAL(DP) :: rcut_small
  
    LOGICAL :: inter_overlap, overlap, accept, accept_or_reject
  
    ! Pair_Energy arrays and Ewald implementation
    INTEGER :: position
    REAL(DP), ALLOCATABLE :: cos_mol_old(:,:), sin_mol_old(:,:)

    reject_type = 0
    E_vdw_move = 0.0_DP
    E_qq_move = 0.0_DP
    E_vdw = 0.0_DP
    E_qq = 0.0_DP
    E_reciprocal_move = 0.0_DP
    inter_overlap = .FALSE.
  
    ! Sum the total number of molecules 
    total_mols = 0 ! sum over species, box
    DO ibox = 1, nbr_boxes
      nmols_box(ibox) = 0
      DO is = 1, nspecies
        ! Only count mobile species
        IF ( max_clus_disp(is,ibox) > 0. ) THEN
          total_mols = total_mols + nmols(is,ibox)
          nmols_box(ibox) = nmols_box(ibox) + nmols(is,ibox)
        END IF
      END DO
    END DO
  
    ! If there are no molecules then return
    IF (total_mols == 0) RETURN
  
    ! If needed, choose a box based on its total mol fraction
    IF(nbr_boxes .GT. 1) THEN
  
      ALLOCATE(x_box(nbr_boxes)) 
  
      DO ibox = 1, nbr_boxes
         x_box(ibox) = REAL(nmols_box(ibox),DP)/REAL(total_mols,DP)
         IF ( ibox > 1 ) THEN
            x_box(ibox) = x_box(ibox) + x_box(ibox-1)
         END IF
      END DO
    
      DO ibox = 1, nbr_boxes
         IF ( rranf() <= x_box(ibox)) EXIT
      END DO
  
      this_box = ibox
      DEALLOCATE(x_box)
  
    ELSE
  
      this_box = 1
  
    END IF
  
    ! If there are no molecules in this box then return
    IF( nmols_box(this_box) == 0 ) RETURN
  
    ! Choose species based on the mol fraction, using Golden sampling
    ALLOCATE(x_species(nspecies))
  
    DO is = 1, nspecies
       IF ( max_clus_disp(is,this_box) > 0. ) THEN
         x_species(is) = REAL(nmols(is,this_box), DP)/REAL(nmols_box(this_box),DP)
       ELSE
         x_species(is) = 0.0_DP
       END IF
       IF ( is > 1 ) THEN
          x_species(is) = x_species(is) + x_species(is-1)
       END IF
    END DO
  
    rand_no = rranf()
    DO is = 1, nspecies
       IF( rand_no <= x_species(is)) EXIT
    END DO
  
    DEALLOCATE(x_species)
  
    ! If the molecule can't move then return
    IF ( max_clus_disp(is,this_box) == 0. ) RETURN
  
    !*********************************************************************************
    !   Step 1) Pick a cluster with uniform probability and its location
    !*********************************************************************************
  
    ! Find Clusters
    CALL Find_Clusters(this_box,2)
  
    ! Choose a cluster (bigger than a monomer) at random
    IF ( MAXVAL(cluster%N) < cluster%M_olig(is) ) RETURN
    !IF ( MAXVAL(cluster%N) < 2. ) RETURN

    iclus = INT ( rranf() * cluster%clusmax ) + 1
    DO WHILE (cluster%N(iclus) < cluster%M_olig(is))
    !DO WHILE (cluster%N(iclus) < 2)

        iclus = INT ( rranf() * cluster%clusmax ) + 1
    END DO
    iclus_N = cluster%N(iclus)
  
    ! update the trial counters
    tot_trials(this_box) = tot_trials(this_box) + 1
    ntrials(is,this_box)%cluster_translate = ntrials(is,this_box)%cluster_translate + 1
  
    ! Generate a random displacement vector. Note that the current formalism will
    ! work for cubic shaped boxes. However, it is easy to extend for nonorthorhombic
    ! boxes where displacements along the basis vectors. 
    dx = ( 2.0_DP * rranf() - 1.0_DP) * max_clus_disp(is,this_box)
    dy = ( 2.0_DP * rranf() - 1.0_DP) * max_clus_disp(is,this_box)
    dz = ( 2.0_DP * rranf() - 1.0_DP) * max_clus_disp(is,this_box)
    IF (int_sim_type == sim_test) THEN
    dx = 1
    dy = 0
    dz = 0
    END IF
  
    !*********************************************************************************
    !   Step 2) Save old energy
    !           obtain the energy of the molecules in the cluster before the move.  Note that due to
    !           this move, the interatomic energies such as vdw and electrostatics will
    !           change. Also the ewald_reciprocal energy will change but there will
    !           be no change in intramolecular energies.
    !*********************************************************************************

    ! Loop over molecules in this cluster
    ALLOCATE( iclus_is(iclus_N), iclus_mol(iclus_N) )
    iclus_mol = 0
    iclus_is = is
    i = 1
    DO c = 1, cluster%n_species_type(2)
        cs = cluster%species_type(2, c)

        DO imol = 1, nmolecules(cs)
            im = locate(imol, cs)
            IF( .NOT. molecule_list(im,cs)%live ) CYCLE 
    
            im_clab = cluster%clabel(im, cs)
            ! IF micelle type then some will not be in 'clustered'
            IF (cluster%criteria(2, int_micelle) .eqv. .TRUE.) THEN
                IF (im_clab == 0) CYCLE
            END IF

            IF ( im_clab == iclus ) THEN
                iclus_mol(i) = im
                iclus_is(i) = cs
                i = i + 1
            END IF

        END DO

    END DO

    CALL Compute_Cluster_Nonbond_Inter_Energy(iclus_N, iclus_mol,iclus_is,E_vdw,E_qq,inter_overlap)

    IF (inter_overlap)  THEN
        WRITE(*,*) 'Disaster, overlap in the old configruation'
        WRITE(*,*) 'Cluster_Translate'
        WRITE(*,*) im, is, this_box
    END IF
          
    DO imol = 1, iclus_N
        im = iclus_mol(imol)
        cs = iclus_is(imol)

        ! IF micelle type then some will not be in 'clustered'
        IF (cluster%criteria(2, int_micelle) .eqv. .TRUE.) THEN
            IF (im == 0) CYCLE
        END IF

        ! Store the old positions of the atoms
        CALL Save_Old_Cartesian_Coordinates(im,cs)

        ! Move atoms by the above vector dx,dy,dz and also update the COM
        atom_list(:,im,cs)%rxp = atom_list(:,im,cs)%rxp + dx
        atom_list(:,im,cs)%ryp = atom_list(:,im,cs)%ryp + dy
        atom_list(:,im,cs)%rzp = atom_list(:,im,cs)%rzp + dz

        molecule_list(im,cs)%xcom = molecule_list(im,cs)%xcom + dx
        molecule_list(im,cs)%ycom = molecule_list(im,cs)%ycom + dy
        molecule_list(im,cs)%zcom = molecule_list(im,cs)%zcom + dz

        CALL Fold_Molecule(im,cs,this_box)

    END DO
  
    ALLOCATE(old_N(SIZE(cluster%N)))
    ALLOCATE( old_clabel(MAXVAL(nmolecules(:)), cluster%n_species_type(2)) )
  
    old_clusmax = cluster%clusmax
    old_N = cluster%N
    old_clabel = cluster%clabel

    !*********************************************************************************
    !   Step 3) Reject if any monomer has become part of a cluster
    !*********************************************************************************

    CALL Find_Clusters(this_box,2)

    ! IF cluster translate added a molecule, so move is rejected
    IF ( iclus_N /= cluster%N(iclus) ) THEN 
        accept = .FALSE.
        reject_type = -1

    !*********************************************************************************
    !   Step 4) Calculate the change in potential energy
    !*********************************************************************************
    ELSE

        CALL Compute_Cluster_Nonbond_Inter_Energy(iclus_N, iclus_mol, iclus_is, E_vdw_move, E_qq_move, inter_overlap)

        IF (inter_overlap) THEN ! Move is rejected

            accept = .FALSE.
            reject_type = -2

        ELSE

            IF ((int_charge_sum_style(this_box) == charge_ewald) .AND. (has_charge(is))) THEN
       
               ALLOCATE(cos_mol_old(nvecs(this_box),SUM(nmolecules)),sin_mol_old(nvecs(this_box),SUM(nmolecules)))
               
               !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
               cos_mol_old(:,:) = cos_mol(1:nvecs(this_box),:)
               sin_mol_old(:,:) = sin_mol(1:nvecs(this_box),:)
               !$OMP END PARALLEL WORKSHARE
       
               CALL Compute_Cluster_Ewald_Reciprocal_Energy_Difference(iclus_N, iclus_mol, iclus_is, int_cluster, this_box, E_reciprocal_move)
               
            END IF
            
            ! Compute the difference in old and new energy
            
            delta_e = ( E_vdw_move - E_vdw ) + ( E_qq_move - E_qq ) + E_reciprocal_move
       
            IF (int_sim_type == sim_nvt_min) THEN
               IF (delta_e  <= 0.0_DP) THEN
                  accept = .TRUE.
               ELSE
                  accept = .FALSE.
                  reject_type = -3
               END IF
            ELSE
       
                ln_pacc = beta(this_box) * delta_e
                accept = accept_or_reject(ln_pacc)
       
            END IF

    !*********************************************************************************
    !   Step 5) Accept or reject the test insertion
    !*********************************************************************************
            IF ( accept ) THEN
       
               !print*, 'clus accep', delta_e, (dx**2.0 + dy**2.0 + dz**2.0)**0.5, iclus_N
               ! accept the move and update the global energies
               energy(this_box)%inter_vdw = energy(this_box)%inter_vdw + E_vdw_move - E_vdw
               energy(this_box)%inter_q   = energy(this_box)%inter_q   + E_qq_move - E_qq
               
               IF(int_charge_sum_style(this_box) == charge_ewald .AND. has_charge(is)) THEN
                  energy(this_box)%ewald_reciprocal = energy(this_box)%ewald_reciprocal + E_reciprocal_move
               END IF
               
               energy(this_box)%total = energy(this_box)%total + delta_e
       
               ! update success counter
               
               nsuccess(is,this_box)%cluster_translate = nsuccess(is,this_box)%cluster_translate + 1
               nequil_success(is,this_box)%cluster_translate = nequil_success(is,this_box)%cluster_translate + 1
       
               IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)
               IF (ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
               IF (ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)

       
            ELSE
               reject_type = -3
       
               IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(im,is,this_box, &
                                              iclus_N, iclus_mol, iclus_is)
       
               IF ((int_charge_sum_style(this_box) == charge_ewald) .AND. (has_charge(is))) THEN
                  ! Also reset the old cos_sum and sin_sum for reciprocal space vectors. Note
                  ! that old vectors were set while difference in ewald reciprocal energy was computed.
                  !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)           
                  cos_sum(:,this_box) = cos_sum_old(:,this_box)
                  sin_sum(:,this_box) = sin_sum_old(:,this_box)
                  cos_mol(1:nvecs(this_box),:) = cos_mol_old(:,:)
                  sin_mol(1:nvecs(this_box),:) = sin_mol_old(:,:)
                  !$OMP END PARALLEL WORKSHARE
                  DEALLOCATE(cos_mol_old,sin_mol_old)
               END IF

            ENDIF
     
        END IF

    END IF

    ! update the cluster counters in anycase so that the clusters don't have to be found twice
    CALL Update_Cluster_Counters(2)

    IF ( .not. accept ) THEN
        ! reset the cluster parameters and positions
        cluster%clusmax = old_clusmax
        cluster%N = old_N

        cluster%clabel = old_clabel
        DEALLOCATE(old_clabel, old_N)

        DO imol = 1, iclus_N
    
            im = iclus_mol(imol)
            cs = iclus_is(imol)

            ! IF micelle type then some will not be in 'clustered'
            IF (cluster%criteria(2, int_micelle) .eqv. .TRUE.) THEN
                IF (im == 0) CYCLE
            END IF

            CALL Revert_Old_Cartesian_Coordinates(im,cs)

        END DO

    END IF

    IF ( MOD(ntrials(is,this_box)%cluster_translate,nupdate) == 0 ) THEN
       IF ( int_run_style == run_equil ) THEN 
          success_ratio = REAL(nequil_success(is,this_box)%cluster_translate,DP)/REAL(nupdate,DP)
       ELSE
          success_ratio = REAL(nsuccess(is,this_box)%cluster_translate,DP)/REAL(ntrials(is,this_box)%cluster_translate,DP)
       END IF
  
       WRITE(logunit,*)
       WRITE(logunit,'(A,I3,A,I1,A,F8.5)')'Success ratio, cluster translation of species ', is , ' in box ', this_box, ' : ', success_ratio
  
       IF ( int_run_style == run_equil ) THEN
  
          ! check if the acceptace is close to 0.5
  
           nequil_success(is,this_box)%cluster_translate = 0
  
           IF  ( success_ratio < 0.00005 ) THEN
               max_clus_disp(is,this_box) = 0.1_DP*max_clus_disp(is,this_box)
           ELSE
               ! minimum max_disp for this species
               IF (has_charge(is)) THEN
                  rcut_small = MIN(rcut_vdw(this_box),rcut_coul(this_box))
               ELSE
                  rcut_small = rcut_vdw(this_box)
               END IF
               max_clus_disp(is,this_box) = MIN(rcut_small,2.0_DP*success_ratio*max_clus_disp(is,this_box))
           END IF
  
           WRITE(logunit,'(A,I3,A,I1,A,F8.5)') 'Maximum width, cluster translation of species ', is,' in box ', this_box, ' : ', max_clus_disp(is,this_box)
          
       END IF
  
    END IF

    DEALLOCATE( iclus_is, iclus_mol )


  END SUBROUTINE Translate_Cluster

  SUBROUTINE Find_Clusters(this_box, count_or_move)

    !*********************************************************************************
    !
    ! Hoshen-Kopelman algorithm for efficiently identifying clusters
    !
    ! 2/19/15  : Andrew P. Santos
    !*********************************************************************************

    INTEGER, INTENT(IN) :: this_box, count_or_move
    INTEGER :: imol, jmol, iatom, i, im, jm, ic
    INTEGER :: is, js, is_clus, js_clus, start, N
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: neigh_list

    cluster%clusmax = 0

    ! cluster label of labels
    ! counts how many atoms are in the cluster
    cluster%N = 0
    cluster%clabel = 0
    DO is = 1, cluster%n_species_type(count_or_move)
        is_clus = cluster%species_type(count_or_move, is)

        IF ( cluster%criteria(count_or_move, int_micelle) .eqv. .TRUE.) THEN
            IF ( is_clus /= cluster%micelle_species) CYCLE
        END IF

        DO imol = 1, nmols(is_clus,this_box)
            !Get a list of the neighbors to an atom in the frame
            im = locate(imol, is_clus)
            ! Make sure that the molecule exists in the simulation
            IF( .NOT. molecule_list(im,is_clus)%live ) CYCLE 

            DO js = 1, cluster%n_species_type(count_or_move)

                js_clus = cluster%species_type(count_or_move, js)

                IF ( cluster%criteria(count_or_move, int_micelle) .eqv. .TRUE.) THEN
                    IF ( js_clus /= cluster%micelle_species) CYCLE
                END IF

                IF (js_clus == is_clus) THEN
                    start = imol+1
                ELSE
                    start = 1
                END IF

                ALLOCATE(neigh_list(nmolecules(js_clus)))
                neigh_list = .FALSE.
                DO jmol = start, nmols(js_clus,this_box)
                    jm = locate(jmol, js_clus)
                    ! Make sure that the molecule exists in the simulation
                    IF( .NOT. molecule_list(jm,js_clus)%live ) CYCLE 

                    neigh_list(jm) = Neighbor(count_or_move, jm, im, js_clus, is_clus)
                END DO

                ! Update cluster label list
                CALL Update_Labels(im, is_clus, js_clus, neigh_list)
                DEALLOCATE(neigh_list)
            END DO

        END DO
    END DO

    IF ( cluster%criteria(count_or_move, int_micelle) .eqv. .TRUE. ) CALL Micelle_Association(this_box, count_or_move)

    !write(*,*) 'c/m', count_or_move
    !write(*,*) 'N', cluster%N(1:200)
    !write(*,*) 'clabel', cluster%clabel
    cluster%n_oligomers = 0
    cluster%n_clusters = 0
    DO is = 1, cluster%n_species_type(count_or_move)
        is_clus = cluster%species_type(count_or_move, is)
        IF ( cluster%criteria(count_or_move, int_micelle) .eqv. .TRUE.) THEN
            ! IF micelle type then skip those that are associated
            IF ( is_clus /= cluster%micelle_species) CYCLE
        END IF

        DO imol = 1, nmols(is_clus,this_box)
            im = locate(imol, is_clus)
            ! Make sure that the molecule exists in the simulation
            IF( .NOT. molecule_list(im,is_clus)%live ) CYCLE 

            ! update the clabels so that they point to the correct cluster not the (-) or it
            DO WHILE ( cluster%N(cluster%clabel(im, is_clus)) < 0 )
                cluster%clabel(im, is_clus) = -cluster%N( cluster%clabel(im, is_clus) )
            END DO

        END DO
    END DO

    IF (count_or_move == 1) THEN
        CALL Update_Cluster_Counters(count_or_move)
    ENDIF

  END SUBROUTINE Find_Clusters

  SUBROUTINE Update_Cluster_Counters(count_or_move)
    INTEGER, INTENT(IN) :: count_or_move
    INTEGER :: is, ic, is_clus

    DO is = 1, cluster%n_species_type(count_or_move)
        is_clus = cluster%species_type(count_or_move, is)
        ic = 1
        DO WHILE ( cluster%N(ic) /= 0 )
    
            IF (cluster%N(ic) > 0) THEN
                ! Tally the aggregation number distribution
                cluster%M( cluster%N(ic) ) = cluster%M( cluster%N(ic) ) + 1
    
                ! Tally up the Nmols of oligomers and clustered
                IF (cluster%N(ic) <= cluster%M_olig(is_clus)) THEN
                    cluster%n_oligomers = cluster%n_oligomers + cluster%N(ic)
                ELSE
                    cluster%n_clusters = cluster%n_clusters + cluster%N(ic)
                END IF
            END IF
    
            ic = ic + 1
        END DO
    END DO
    !write(*,*) 'N', cluster%N(1:200)
    !write(*,*) 'clabel', cluster%clabel(1:500, 1)
    !write(*,*) 'M', cluster%M
    !write(*,*) 'olig, clus', cluster%n_oligomers, cluster%n_clusters

  END SUBROUTINE Update_Cluster_Counters

  SUBROUTINE Update_Labels(imol, is, js, neigh_list)

    !*********************************************************************************
    !
    ! Part of the Hoshen-Kopelman, where the labels of atoms/molecules in each cluster
    ! is updated
    !
    ! 2/19/15  : Andrew P. Santos
    !*********************************************************************************
    INTEGER, INTENT(IN) :: imol, is, js
    LOGICAL, INTENT(IN), DIMENSION(nmolecules(js)) :: neigh_list
    INTEGER :: iclus, nclus, ineigh
    INTEGER :: max_clus, min_clus

    DO ineigh = 1, nmolecules(js)
        IF (neigh_list(ineigh) .eqv. .FALSE.) CYCLE

        ! current site' cluster label
        iclus = cluster%clabel(imol, is)
        ! neighboring site's cluster label
        nclus = cluster%clabel(ineigh, js)
        ! IF the neighbor site has been assigned
        IF (nclus /= 0) THEN
            ! IF the current site has yet to be assigned
            IF (iclus == 0) THEN
                ! assign the cluster to the occupied neighbor's cluster
                cluster%clabel(imol, is) = nclus
                DO WHILE (cluster%N(nclus) < 0)
                    nclus = -cluster%N(nclus)
                END DO

                cluster%N(nclus) = cluster%N(nclus) + 1

            ! IF the current site is defined
            ELSE
                ! assign both labels with the larger value to the lower
                DO WHILE (cluster%N(nclus) < 0)
                    nclus = -cluster%N(nclus)
                END DO

                DO WHILE (cluster%N(iclus) < 0)
                    iclus = -cluster%N(iclus)
                END DO

                min_clus = min(nclus, iclus)
                max_clus = max(nclus, iclus)
                ! assign the cluster lower label value
                cluster%clabel(imol, is) = min_clus
                IF (min_clus /= max_clus) THEN
                    cluster%N(min_clus) = cluster%N(min_clus) + cluster%N(max_clus)
                    cluster%N(max_clus) = -min_clus
                END IF

            END IF

        ! IF the neighbor bead has not been assigned
        ELSE 
            ! IF the current site has yet to be assigned
            IF (iclus == 0) THEN
                ! assign site a new cluster label
                cluster%clusmax = cluster%clusmax + 1
                cluster%clabel(imol, is) = cluster%clusmax
                iclus = cluster%clusmax
                cluster%N(cluster%clusmax) = 1
            END IF

            cluster%clabel(ineigh, js) = iclus

            DO WHILE (cluster%N(iclus) < 0)
                iclus = -cluster%N(iclus)
            END DO

            cluster%N(iclus) = cluster%N(iclus) + 1
        END IF
    END DO

    IF (cluster%clabel(imol, is) == 0) THEN
        cluster%clusmax = cluster%clusmax + 1
        cluster%clabel(imol, is) = cluster%clusmax
        cluster%N(cluster%clusmax) = 1
    END IF

  END SUBROUTINE Update_Labels

  FUNCTION Neighbor(count_or_move, test_part, cur_part, test_type, cur_type)

    !*********************************************************************************
    !
    ! Return T/F if two atoms/molecules are neighbors
    !
    ! 2/19/15  : Andrew P. Santos
    !*********************************************************************************

    LOGICAL :: Neighbor
    INTEGER, INTENT(IN) :: test_part, cur_part, test_type, cur_type, count_or_move
    REAL(DP) :: rxij, ryij, rzij, rijsq, rxijp, ryijp, rzijp
    INTEGER :: test_atom, cur_atom, n_accept

    Neighbor = .FALSE.

    IF (.not. ANY(cluster%species_type(count_or_move,:) == test_type)) THEN
        RETURN
    ELSE IF( .NOT. molecule_list(test_part,test_type)%live ) THEN
        RETURN
    END IF
      
    IF (cluster%criteria(count_or_move, int_com)) THEN
        ! Get the positions of the COM of the two molecule species
        rxijp = molecule_list(test_part,test_type)%xcom - molecule_list(cur_part, cur_type)%xcom
        ryijp = molecule_list(test_part,test_type)%ycom - molecule_list(cur_part, cur_type)%ycom
        rzijp = molecule_list(test_part,test_type)%zcom - molecule_list(cur_part, cur_type)%zcom
        
        ! Now get the minimum image separation 
        CALL Minimum_Image_Separation(1, rxijp, ryijp, rzijp, rxij, ryij, rzij)

        rijsq = rxij*rxij + ryij*ryij + rzij*rzij
        
        IF (rijsq < cluster%min_distance_sq(count_or_move, cur_type, test_type, 0, 0 ))THEN
           Neighbor = .TRUE.
           RETURN
        END IF
                      
    END IF

    IF (cluster%criteria(count_or_move, int_type)) THEN
        DO test_atom = 1 , natoms(test_type)

            DO cur_atom = 1 , natoms(cur_type)
                IF (cluster%min_distance_sq(count_or_move, test_type, cur_type, test_atom, cur_atom) < 0.000001) CYCLE

                ! Get the positions of the COM of the two molecule species
                rxijp = atom_list(test_atom, test_part, test_type)%rxp - &
                        atom_list(cur_atom, cur_part, cur_type)%rxp
                ryijp = atom_list(test_atom, test_part, test_type)%ryp - &
                        atom_list(cur_atom, cur_part, cur_type)%ryp
                rzijp = atom_list(test_atom, test_part, test_type)%rzp - &
                        atom_list(cur_atom, cur_part, cur_type)%rzp
                
                ! Now get the minimum image separation 
                CALL Minimum_Image_Separation(1, rxijp, ryijp, rzijp, rxij, ryij, rzij)
        
                rijsq = rxij*rxij + ryij*ryij + rzij*rzij
                
                IF ( rijsq < cluster%min_distance_sq(count_or_move, test_type, cur_type, test_atom, cur_atom) ) THEN
                   Neighbor = .TRUE.
                   RETURN
                END IF

            END DO
        END DO
    END IF
    IF (cluster%criteria(count_or_move, int_skh)) THEN
        ! Get the positions of the COM of the two molecule species
        rxijp = molecule_list(test_part,test_type)%xcom - molecule_list(cur_part, cur_type)%xcom
        ryijp = molecule_list(test_part,test_type)%ycom - molecule_list(cur_part, cur_type)%ycom
        rzijp = molecule_list(test_part,test_type)%zcom - molecule_list(cur_part, cur_type)%zcom
        
        ! Now get the minimum image separation 
        CALL Minimum_Image_Separation(1, rxijp, ryijp, rzijp, rxij, ryij, rzij)

        rijsq = rxij*rxij + ryij*ryij + rzij*rzij
        
        n_accept = 0.0 
        IF (rijsq < cluster%r1_sq(count_or_move, cur_type, 0)) THEN
           Neighbor = .TRUE.
           RETURN
        ELSE IF (rijsq < cluster%r2_sq(count_or_move, cur_type, 0)) THEN
           n_accept = n_accept + 1.5
        ELSE IF (rijsq < cluster%r3_sq(count_or_move, cur_type, 0)) THEN
           n_accept = n_accept + 1.0
        END IF

        DO test_atom = 1 , natoms(test_type)
            IF (cluster%r3_sq(count_or_move, test_type, test_atom) < 0.000001) CYCLE

            DO cur_atom = 1 , natoms(cur_type)
                IF (cluster%r3_sq(count_or_move, cur_type, cur_atom) < 0.000001) CYCLE

                ! Get the positions of the COM of the two molecule species
                rxijp = atom_list(test_atom, test_part, test_type)%rxp - &
                        atom_list(cur_atom, cur_part, cur_type)%rxp
                ryijp = atom_list(test_atom, test_part, test_type)%ryp - &
                        atom_list(cur_atom, cur_part, cur_type)%ryp
                rzijp = atom_list(test_atom, test_part, test_type)%rzp - &
                        atom_list(cur_atom, cur_part, cur_type)%rzp
                
                ! Now get the minimum image separation 
                CALL Minimum_Image_Separation(1, rxijp, ryijp, rzijp, rxij, ryij, rzij)
        
                rijsq = rxij*rxij + ryij*ryij + rzij*rzij
                
                IF (rijsq < cluster%r1_sq(count_or_move, cur_type, 0)) THEN
                   Neighbor = .TRUE.
                   RETURN
                ELSE IF (rijsq < cluster%r2_sq(count_or_move, cur_type, 0)) THEN
                   n_accept = n_accept + 1.5
                ELSE IF (rijsq < cluster%r3_sq(count_or_move, cur_type, 0)) THEN
                   n_accept = n_accept + 1.0
                ELSE IF (n_accept >= 3.0) THEN
                   Neighbor = .TRUE.
                   RETURN
                END IF

            END DO
        END DO
    END IF

  END FUNCTION Neighbor

  SUBROUTINE Micelle_Association(this_box, c_or_m)

    !*********************************************************************************
    !
    ! Add molecules which are "associated" with a cluster, but cannot connect clusters
    !
    ! 2/19/15  : Andrew P. Santos
    !*********************************************************************************

    INTEGER, INTENT(IN) :: this_box, c_or_m
    INTEGER :: i, j, as, am, cm, cs, as_clus, nclus
    INTEGER :: ierr, line_nbr, nbr_entries
    REAL(DP) :: rxij, ryij, rzij, rijsq, rxijp, ryijp, rzijp

    cluster%criteria(c_or_m, int_com) = .FALSE.
    cluster%criteria(c_or_m, int_type) = .TRUE.

    cs = cluster%micelle_species
    
    DO i = 1, nmolecules(cs)
        cm = locate(i, cs)
        IF( .NOT. molecule_list(cm, cs)%live ) CYCLE

        DO as = 1, cluster%n_species_type(c_or_m)
            as_clus = cluster%species_type(c_or_m, as)

            IF ( as_clus == cluster%micelle_species) CYCLE

            DO j = 1, nmolecules(as_clus)
                am = locate(j, as_clus)
                IF( .NOT. molecule_list(am, as_clus)%live ) CYCLE
    
                IF (Neighbor(c_or_m, am, cm, as_clus, cs) .eqv. .true.) THEN
                    ! Add it to the clabel and N
                    cluster%clabel(am, as_clus) = cluster%clabel(cm, cs)
                    nclus = cluster%clabel(cm, cs)

                    DO WHILE (cluster%N(nclus) < 0)
                        nclus = -cluster%N(nclus)
                    END DO

                    cluster%N(nclus) = cluster%N(nclus) + 1

                END IF
    
            END DO
        END DO
    END DO

    cluster%criteria(c_or_m, int_com) = .TRUE.
    cluster%criteria(c_or_m, int_type) = .FALSE.
                    
  END SUBROUTINE Micelle_Association

  SUBROUTINE Compute_Cluster_Nonbond_Inter_Energy(nclus_mol,clus_mol,clus_is,E_inter_vdw,E_inter_qq,overlap)

    !**************************************************************************************************
    ! This subroutine computes interatomic LJ and charge interactions as well as virials associated
    ! with these interactions. 
    !
    ! CALLS
    ! 
    ! Minimum_Image_Separation
    ! Pair_Energy
    !
    ! CALLED BY
    !
    ! Translate_Cluster
    !
    ! Written by Jindal Shah on 12/07/07
    ! Adapted by Andrew Santos on 05/16
    !***************************************************************************************************

    IMPLICIT NONE

    !$ include 'omp_lib.h'

    INTEGER, INTENT(IN):: nclus_mol
    INTEGER, INTENT(IN):: clus_mol(:), clus_is(:)
    REAL(DP), INTENT(OUT) :: E_inter_vdw, E_inter_qq
    LOGICAL :: overlap
    !---------------------------------------------------------------------------------------------------

    INTEGER  :: ispecies, imolecule, this_box, this_locate, im_test, is_test, imc_test
    
    REAL(DP) :: Eij_vdw, Eij_qq
    REAL(DP) :: rcom, rx, ry, rz

    LOGICAL :: get_interaction

    INTEGER :: im, ic_mol, is

    LOGICAL :: my_overlap, shared_overlap

    E_inter_vdw = 0.0_DP
    E_inter_qq = 0.0_DP
    overlap = .FALSE.
    my_overlap = .FALSE.
    shared_overlap = .false.
    
    !print*, 'im, is, E_inter_vdw, E_inter_qq'
    !print*, 'clumol', clus_mol
    clusMolLoop: DO ic_mol = 1, nclus_mol

        im = clus_mol(ic_mol)
        IF (im == 0) CYCLE

        is = clus_is(ic_mol)
        this_box = molecule_list(im,is)%which_box

        speciesLoop: DO ispecies = 1, nspecies
           
           !$OMP PARALLEL DO DEFAULT(SHARED) &
           !$OMP PRIVATE(imolecule,this_locate,get_interaction) &
           !$OMP PRIVATE(rx,ry,rz,rcom,Eij_vdw,Eij_qq) &
           !$OMP PRIVATE(imc_test, is_test, im_test) &
           !$OMP SCHEDULE(DYNAMIC) &
           !$OMP REDUCTION(+:E_inter_vdw,E_inter_qq) & 
           !$OMP REDUCTION(.OR.:my_overlap)  
    
           
           moleculeLoop: DO imolecule = 1, nmolecules(ispecies)
              
              IF(shared_overlap) CYCLE
              
              this_locate = locate(imolecule,ispecies)
              
              IF (molecule_list(this_locate,ispecies)%which_box /= this_box) CYCLE moleculeLOOP
              
              ! make sure that the molecule is currently part of the simulation
              
              IF(.NOT. molecule_list(this_locate,ispecies)%live) CYCLE moleculeLOOP
              
              ! Skip those that are included in the clusters
              
              DO imc_test = 1, nclus_mol
                  im_test = clus_mol(imc_test)
                  is_test = clus_is(imc_test)
                  IF( im_test == this_locate .AND. ispecies == is_test) CYCLE moleculeLOOP
              END DO
              
              ! Determine if any atoms of these two molecules will interact
              CALL Check_Interaction(im,is,this_locate,ispecies,get_interaction,rcom,rx,ry,rz) 
    
              IF (.NOT. get_interaction) CYCLE moleculeLOOP       
              
              
              CALL Compute_Molecule_Pair_Interaction(im,is,this_locate,ispecies,this_box, &
                   Eij_vdw,Eij_qq,my_overlap)
              
              IF( my_overlap) shared_overlap = .true.
              
              E_inter_vdw = E_inter_vdw + Eij_vdw
              E_inter_qq  = E_inter_qq + Eij_qq
              
           END DO moleculeLoop
           !$OMP END PARALLEL DO
           
           IF(shared_overlap) THEN
              overlap = .true.
              RETURN
           ENDIF
           
        !print*, im, ispecies, E_inter_vdw, E_inter_qq
        END DO speciesLoop
    END DO clusMolLoop
    
  END SUBROUTINE Compute_Cluster_Nonbond_Inter_Energy

  SUBROUTINE Compute_Cluster_Ewald_Reciprocal_Energy_Difference(nclus_mol, clus_mol,clus_is,move_flag,this_box,v_recip_difference)
    !************************************************************************************************
    ! The subroutine computes the difference in Ewald reciprocal space energy for a given move.
    !/
    ! We will develop this routine for a number of moves.
    !
    ! Translation of COM
    ! Rotation about COM
    ! Angle Distortion
    ! Rigid Dihedral rotation
    ! Molecule insertion
    ! Molecule Deletion
    !***********************************************************************************************

    USE Type_Definitions
    USE Run_Variables
    
    IMPLICIT NONE

!    !$ include 'omp_lib.h'

    INTEGER, INTENT(IN):: nclus_mol
    INTEGER, INTENT(IN):: clus_mol(:), clus_is(:)
    INTEGER, INTENT(IN) :: move_flag, this_box
    REAL(DP), INTENT(OUT) :: v_recip_difference

    ! Local variables
    
    INTEGER :: i, ia, im, ic_mol, is, this_locate
    INTEGER, ALLOCATABLE :: sum_nmolec(:)

    REAL(DP) :: hdotr_new

    REAL(DP) :: cos_sum_im, sin_sum_im

    ! storage stuff

    ALLOCATE(sum_nmolec(nspecies))

    sum_nmolec(1) = 0

    DO is = 2, nspecies

       sum_nmolec(is) = SUM(nmolecules(1:is-1))

    END DO

    ! get the location of im 

    v_recip_difference = 0.0_DP

    !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
    cos_sum_old(1:nvecs(this_box),this_box) = cos_sum(1:nvecs(this_box),this_box)
    sin_sum_old(1:nvecs(this_box),this_box) = sin_sum(1:nvecs(this_box),this_box)
    !$OMP END PARALLEL WORKSHARE

    IF ( move_flag == int_cluster) THEN


       ! Cluster move. Therefore, the contribution of cos(hdotr) and
       ! sin(hdotr) of the old coordinates will be subtracted off for each of reciprocal vectors
       ! and corresponding terms for the new coordinates are added for all molecules in cluster

       ! Note that the flag INTRA will refer to any of the moves that correspond to the 
       ! intramolecular DOF change.
       
           
       !$OMP PARALLEL DO DEFAULT(SHARED) &
       !$OMP PRIVATE(i,ia,cos_sum_im,sin_sum_im) &
       !$OMP PRIVATE(hdotr_new) &
       !$OMP SCHEDULE(STATIC) &
       !$OMP REDUCTION(+:v_recip_difference)
       DO i = 1, nvecs(this_box)

          clusMolLoop: DO ic_mol = 1, nclus_mol

             im = clus_mol(ic_mol)
             IF (im == 0) CYCLE

             is = clus_is(ic_mol)

             this_locate = locate(im,is) + sum_nmolec(is)

             cos_sum_im = 0.0_DP
             sin_sum_im = 0.0_DP

             DO ia = 1,natoms(is)
   
                IF (nonbond_list(ia,is)%charge == 0) CYCLE
   
                ! let us compute the old and new hdotr
   
                hdotr_new = hx(i,this_box) * atom_list(ia,im,is)%rxp + &
                            hy(i,this_box) * atom_list(ia,im,is)%ryp + &
                            hz(i,this_box) * atom_list(ia,im,is)%rzp
                
                cos_sum_im = cos_sum_im + nonbond_list(ia,is)%charge * DCOS(hdotr_new)
                sin_sum_im = sin_sum_im + nonbond_list(ia,is)%charge * DSIN(hdotr_new)
   
             END DO
   
             cos_sum(i,this_box) = cos_sum(i,this_box) + cos_sum_im - cos_mol(i,this_locate)
             sin_sum(i,this_box) = sin_sum(i,this_box) + sin_sum_im - sin_mol(i,this_locate)
   
             ! set the molecules cos and sin terms to the one calculated here
             cos_mol(i,this_locate) = cos_sum_im
             sin_mol(i,this_locate) = sin_sum_im

          END DO clusMolLoop

          v_recip_difference = v_recip_difference + cn(i,this_box) * (cos_sum(i,this_box) * &
                             cos_sum(i,this_box) + sin_sum(i,this_box) * sin_sum(i,this_box))
   
       END DO
       !$OMP END PARALLEL DO

       v_recip_difference = v_recip_difference*charge_factor(this_box) - &
                            energy(this_box)%ewald_reciprocal

       RETURN

    END IF

  END SUBROUTINE Compute_Cluster_Ewald_Reciprocal_Energy_Difference
  !********************************************************************************************

END MODULE Cluster_Routines
