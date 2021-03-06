Usage: .\calc_rmsd.exe [reference.xyz] [trajectory.xyz]

Calculates RMSD along a trajectory (multi-frame) file based on a reference file.

Algorithm:
Kabsch's algorithm of aligning is applied.
let P be the ref, Q be the traj, both dimension(n, d), where n is the amount of atoms, d is the number of dimensions.
first translate geometry center of P and Q to the origin point of the cartesian coordinates:
P -= sum(P) / n, Q -= sum(Q) / n
then rotate to align
let the covariance matrix H be P^T @ Q, then do sigular value decomposition to H:
H = U @ S @ V^T, after the SVD, then let sign be sign(det(U @ V^T)), and let 
         1   0   0
R = U @ (0   1   0  ) @ V^T, then R is the rotation matrix.
         0   0  sign
assign Q = Q @ R^T, then the two molecules are aligned.
The RMSD is defined as sqrt(sum_i{sum((P_i - Q_i) ** 2) / n)
The average of RMSD along the trajectory is E_RMSD = sum(RMSD) / f, where f is the amount of frames, 
and the sample standard deviation is defined as sqrt(sum((RMSD - E_RMSD) ** 2) / (f - 1)).

The examples are ref.xyz as the reference file and traj.xyz as the trajectory file, 
you can get the same result of RMSD using VMD. For how to do this in VMD you can check http://sobereva.com/290.

