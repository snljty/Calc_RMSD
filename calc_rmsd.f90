! Calculates RMSD along a trajectory (multi-frame) file based on a reference file.
!
! Algorithm:
! let P be the ref, Q be the traj, both dimension(n, d), where n is the amount of atoms, d is the number of dimensions.
! first translate geometry center of P and Q to the origin point of the cartesian coordinates:
! P -= sum(P) / n, Q -= sum(Q) / n
! then rotate to align
! let the covariance matrix H be P^T @ Q, then do sigular value decomposition to H:
! H = U @ S @ V^T, after the SVD, then let sign be det(U @ V^T) // abs(det(U @ V^T)), then let 
!          1   0   0
! R = U @ (0   1   0  ) @ V^T, then R is the rotation matrix.
!          0   0  sign
! assign Q = Q @ R^T, then the two molecules are aligned.
! The RMSD is defined as sqrt(sum_i{sum((P_i - Q_i) ** 2) / n)
! The average of RMSD along the trajectory is E_RMSD = sum(RMSD) / f, where f is the amount of frames, 
! and the sample standard deviation is defined as sqrt(sum((RMSD - E_RMSD) ** 2) / (f - 1)).

# define stdin_unit 5
# define stdout_unit 6
# define stderr_unit 0

# define exit_success 0
# define exit_failure 1

program main
    implicit none
    integer :: argc, iarg
    integer, dimension(:), allocatable :: argl
    character(len=:), dimension(:), allocatable :: argv
    integer, parameter :: num_dims = 3
    integer :: num_atoms_ref, num_atoms_traj
    character(len=128) :: buf
    character(len=2), dimension(:), allocatable :: elems_ref, elems_traj
    real(kind=8), dimension(:, :), allocatable :: coords_ref, coords_traj
    character(len=128) :: ifl_ref_name, ifl_traj_name
    character(len=128), parameter :: ofl_result_name = "RMSD.txt"
    integer, parameter :: ifl_ref_unit = 11, ifl_traj_unit = 13
    integer, parameter :: ofl_result_unit = 12
    integer :: ind_atom
    integer :: io_status
    integer :: ind_frame, num_frames
    real(kind=8), dimension(:), allocatable :: rmsd_results
    real(kind=8), dimension(:, :), allocatable :: coords_tmp
    real(kind=8), external :: calc_rmsd
    real(kind=8) :: rmsd_ave

    argc = command_argument_count()
    allocate(argl(0:argc))
    do iarg = 0, argc
        call get_command_argument(iarg, length = argl(iarg))
    end do
    allocate(character(maxval(argl)) :: argv(0:argc))
    deallocate(argl)
    do iarg = 0, argc
        call get_command_argument(iarg, argv(iarg))
    end do

    if (argc > 2) then
        write(stderr_unit, "(a,i0,a)") "Error! At most 2 command arguments should be provided, but got ", argc, "."
        stop exit_failure
    else if (argc == 2) then
        ifl_ref_name = trim(argv(1))
        ifl_traj_name = trim(argv(2))
    else if (argc == 1) then
        write(*, "(a,1x,a,1x,a,1x,a)") "Usage:", trim(argv(0)), "[reference.xyz]", "[trajectory.xyz]"
        write(*, "(a)") "The trajectory xyz file could be multi-frams."
        stop exit_success
    else
        write(*, "(a)") "Input name of the reference xyz file:"
        read(*, "(a)") ifl_ref_name
        if (ifl_ref_name(1:1) == """") then
            ifl_ref_name(len_trim(ifl_ref_name):) = ""
            ifl_ref_name = ifl_ref_name(2:)
        end if
        write(*, "(a)") "Input name of the trajectory xyz file (which could be multi-frame):"
        read(*, "(a)") ifl_traj_name
        if (ifl_traj_name(1:1) == """") then
            ifl_traj_name(len_trim(ifl_traj_name):) = ""
            ifl_traj_name = ifl_traj_name(2:)
        end if
    end if
    deallocate(argv)

    if (index(ifl_ref_name, ".xyz", back = .true.) + len(".xyz") - 1 /= len_trim(ifl_ref_name)) then
        write(stderr_unit, "(a)") "Error! Suffix of the reference file should be "".xyz""!"
        stop exit_failure
    end if
    if (index(ifl_traj_name, ".xyz", back = .true.) + len(".xyz") - 1 /= len_trim(ifl_traj_name)) then
        write(stderr_unit, "(a)") "Error! Suffix of the trajectory file should be "".xyz""!"
        stop exit_failure
    end if

    open(ifl_ref_unit, file = trim(ifl_ref_name), status = "old", action = "read")
    read(ifl_ref_unit, "(a)") buf
    read(buf, *) num_atoms_ref
    allocate(elems_ref(num_atoms_ref))
    allocate(coords_ref(num_atoms_ref, num_dims))
    allocate(elems_traj(num_atoms_ref))
    allocate(coords_traj(num_atoms_ref, num_dims))
    allocate(coords_tmp(num_atoms_ref, num_dims))
    read(ifl_ref_unit, "(a)") buf
    do ind_atom = 1, num_atoms_ref
        read(ifl_ref_unit, "(a)") buf
        read(buf, *) elems_ref(ind_atom), coords_ref(ind_atom, :)
    end do
    call move_to_center(num_atoms_ref, coords_ref)
    close(ifl_ref_unit)

    ind_frame = 0
    open(ifl_traj_unit, file = trim(ifl_traj_name), status = "old", action = "read")
    do while (.true.)
        read(ifl_traj_unit, "(a)", iostat = io_status) buf
        if (io_status /= 0) then
            exit
        end if
        if (buf == "") then
            exit
        end if
        ind_frame = ind_frame + 1
        read(buf, *) num_atoms_traj
        if (num_atoms_traj /= num_atoms_ref) then
            write(stderr_unit, "(a,i0,a,i0,a,i0,a)") "Error! Frame ", ind_frame, &
                ": number of atoms mismatches with the reference (", num_atoms_traj, " and ", num_atoms_ref, ")."
            stop exit_failure
        end if
        read(ifl_traj_unit, "(a)") buf
        do ind_atom = 1, num_atoms_traj
            read(ifl_traj_unit, "(a)") buf
        end do
    end do
    close(ifl_traj_unit)
    num_frames = ind_frame
    if (num_frames == 1) then
        write(*, "(a)") "There is only one frame in the trajectory file."
    else
        write(*, "(a,i0,a)") "Total ", num_frames, " frames in the trajectory file."
    end if
    allocate(rmsd_results(num_frames))

    open(ifl_traj_unit, file = trim(ifl_traj_name), status = "old", action = "read", position = "rewind")
    do ind_frame = 1, num_frames
        read(ifl_traj_unit, "(a)") buf
        ! read(buf, *) num_atoms_traj
        read(ifl_traj_unit, "(a)") buf
        do ind_atom = 1, num_atoms_traj
            read(ifl_traj_unit, "(a)") buf
            read(buf, *) elems_traj(ind_atom), coords_traj(ind_atom, :)
        end do
        call move_to_center(num_atoms_traj, coords_traj)
        call rotate_to_align(num_atoms_ref, coords_ref, coords_traj, coords_tmp)
        rmsd_results(ind_frame) = calc_rmsd(num_atoms_ref, coords_ref, coords_traj)
    end do
    close(ifl_traj_unit)

    open(ofl_result_unit, file = trim(ofl_result_name), status = "old", action = "read", iostat = io_status)
    if (io_status == 0) then
        close(ofl_result_unit)
        write(stderr_unit, "(a,a,a)") "Warning: file """, trim(ofl_result_name), """ already exists."
        write(stderr_unit, "(a,a,a,a)") "It will be renamed to """, trim(ofl_result_name), ".bak", """."
        open(ofl_result_unit, file = trim(ofl_result_name) // ".bak", status = "old", action = "read", iostat = io_status)
        if (io_status == 0) then
            write(stderr_unit, "(a,a,a)") "Warning: file """, trim(ofl_result_name) // ".bak", """ already exists."
            write(stderr_unit, "(a)") "It will be removed."
            close(ofl_result_unit, status = "delete")
        end if
        call rename(trim(ofl_result_name), trim(ofl_result_name) // ".bak")
    end if
    rmsd_ave = sum(rmsd_results) / num_frames
    if (num_frames == 1) then
        write(*, "(a,f0.3)") "RMSD: ", rmsd_results(1)
    else
        write(*, "(a,f0.3,a,f0.3,a,f0.3,a,f0.3)") "RMSDs:    minimum = ", minval(rmsd_results), &
        ", maximum = ", maxval(rmsd_results), ", average = ", rmsd_ave , &
        ", sample standard deviation = ", sqrt(sum((rmsd_results - rmsd_ave) ** 2) / (num_frames - 1))
    end if
    open(ofl_result_unit, file = trim(ofl_result_name), status = "replace", action = "write", position = "rewind")
    write(ofl_result_unit, "(f12.6)") rmsd_results
    close(ofl_result_unit)
    write(*, "(a,a,a)") "The results have been saved to """, trim(ofl_result_name), """."

    deallocate(elems_ref)
    deallocate(coords_ref)
    deallocate(elems_traj)
    deallocate(coords_traj)
    deallocate(coords_tmp)
    deallocate(rmsd_results)

    stop
end program main

subroutine move_to_center(num_atoms, coords)
    implicit none
    integer, parameter :: num_dims = 3
    integer, intent(in) :: num_atoms
    real(kind=8), dimension(num_atoms, num_dims), intent(inout) :: coords
    real(kind=8), dimension(num_dims) :: ave
    integer :: ind_atom

    ave = sum(coords, dim = 1) / num_atoms
    do ind_atom = 1, num_atoms
        coords(ind_atom, :) = coords(ind_atom, :) - ave
    end do

    return
end subroutine move_to_center

subroutine rotate_to_align(num_atoms, coords1, coords2, coords0)
    ! assumed they are already moved to center
    implicit none
    integer, parameter :: num_dims = 3
    integer, intent(in) :: num_atoms
    real(kind=8), dimension(num_atoms, num_dims), intent(in) :: coords1
    real(kind=8), dimension(num_atoms, num_dims), intent(inout) :: coords2
    real(kind=8), dimension(num_atoms, num_dims), intent(out) :: coords0 ! temporary
    real(kind=8), dimension(num_dims, num_dims) :: h, u, v_t
    real(kind=8), dimension(num_dims) :: s
    real(kind=8), dimension(:), allocatable :: work
    integer :: lwork
    integer :: info
    integer, dimension(num_dims) :: ipiv
    real(kind=8), dimension(num_atoms, num_dims) :: r

    call dgemm("T", "N", num_dims, num_dims, num_atoms, 1.d0, coords1, num_atoms, coords2, num_atoms, 0.d0, h, num_dims)
    lwork = -1
    allocate(work(- lwork))
    call dgesvd("A", "A", num_dims, num_dims, h, num_dims, s, u, num_dims, v_t, num_dims, work, lwork, info)
    lwork = nint(work(1))
    deallocate(work)
    allocate(work(lwork))
    call dgesvd("A", "A", num_dims, num_dims, h, num_dims, s, u, num_dims, v_t, num_dims, work, lwork, info)
    lwork = 0
    deallocate(work)
    ! as h is destroied, below h is a temporary matrix for determining the sign of the determination
    call dgemm("N", "N", num_dims, num_dims, num_dims, 1.d0, u, num_dims, v_t, num_dims, 0.d0, h, num_dims)
    call dgetrf(num_dims, num_dims, h, num_dims, ipiv, info)
    if (h(1, 1) * h(2, 2) * h(3, 3) < 0.d0) then ! det(U @ V^T) < 0
        v_t(num_dims, :) = - v_t(num_dims, :)
    end if
    call dgemm("N", "N", num_dims, num_dims, num_dims, 1.d0, u, num_dims, v_t, num_dims, 0.d0, r, num_dims)
    call dlacpy("A", num_atoms, num_dims, coords2, num_atoms, coords0, num_atoms)
    call dgemm("N", "T", num_atoms, num_dims, num_dims, 1.d0, coords0, num_atoms, r, num_dims, 0.d0, coords2, num_atoms)

    return
end subroutine rotate_to_align

real(kind=8) function calc_rmsd(num_atoms, coords1, coords2)
    ! assumed they are already aligned
    implicit none
    integer, parameter :: num_dims = 3
    integer, intent(in) :: num_atoms
    real(kind=8), dimension(num_atoms, num_dims), intent(in) :: coords1
    real(kind=8), dimension(num_atoms, num_dims), intent(in) :: coords2

    calc_rmsd = sqrt(sum((coords1 - coords2) ** 2) / num_atoms)

    return
end function calc_rmsd


