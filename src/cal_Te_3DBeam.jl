function cal_Te_3DBeam(xi,yi,zi,xj,yj,zj,xk,yk,zk)

    vec_ij = [ xj - xi
               yj - yi
               zj - zi ]
    e1 = (vec_ij)/norm(vec_ij)

    vec_ik = [ xk - xi
               yk - yi
               zk - zi ]
    z_bar = cross(e1,vec_ik)
    e3 = z_bar/norm(z_bar)

    e2 = cross(e1,e3)

    T_blk = [ e1  e2  e3 ]

    Zero_blk = zeros(3,3)

    Te = transpose([    T_blk  Zero_blk  Zero_blk  Zero_blk
                     Zero_blk     T_blk  Zero_blk  Zero_blk
                     Zero_blk  Zero_blk     T_blk  Zero_blk
                     Zero_blk  Zero_blk  Zero_blk     T_blk ])

    return Te

end
