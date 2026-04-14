import numpy as np
import math

def AHBASE_load_auxbasis(par,aux):
    global nprim, ncoeff, prim_order, prim_expon
    nprim = 0
    ncoeff = 0
    
    prim_expon = []
    prim_order = []
    index = 0

    for name in par:
        aux.seek(0)
        if 'auxbasis' in name:
            auxname = name[9:]
            for basis in aux:
                if auxname in basis:
                    index = index+1
                    nshell = int(aux.readline())
                    for i in range(nshell):
                        oprimshell = aux.readline().split()
                        nprimshell = int(oprimshell[0])
                        ordprim = int(oprimshell[1])
                        if ordprim == 0:
                            order = 1
                        elif ordprim == 1:
                            order = 4
                        elif ordprim == -1:
                            order = 3
                        elif ordprim == -2:
                            order = 6
                        elif ordprim == 2:
                            order = 10
                        elif ordprim == 3:
                            order = 20
                        else:
                            print('bad ordprim value',ordprim)
                        nprim = nprim + nprimshell
                        ncoeff = ncoeff + order*nprimshell
                        for j in range(nprimshell):
                            ex = float(aux.readline())
                            prim_expon.append(ex)
                            prim_order.append(order)
    return 

def fac(n):
    if n < 0:
        return 1
    return math.factorial(n)

def dfac(n):
    if n <= 1:
        return 1.0
    if n%2 == 1:
        k = (n+1)//2
        return math.factorial(2*k) / (2**k * math.factorial(k))
    else:
        k = n//2
        return 2**k * math.factorial(k)
    
def auxfun(n,i,j):
    if (2*(n-i-j)-1) < 0:
        tmp = 1
    else:
        tmp = dfac(2*(n-i-j)-1)
    auxfun = (-1)**(i+j)*tmp/(fac(i)*fac(j)*fac(n-2*i)*fac(n-2*j))
    return auxfun
           
    
def normalize_aux(nx,ny,nz,ex):
    lc = nx+ny+nz
    lnorm = np.sqrt(2)*(np.pi/ex)**(2.5)/(4*ex)**lc
    factor = lnorm*(fac(nx)*fac(ny)*fac(nz))**2
    
    norm = 0
    for iax in range(1+int((nx)/2)):
        for jax in range(1+int((nx)/2)):
            wax = auxfun(nx,iax,jax)
            for iay in range(1+int((ny)/2)):
                for jay in range(1+int((ny)/2)):
                    way = auxfun(ny,iay,jay)
                    for iaz in range(1+int((nz)/2)):
                        for jaz in range(1+int((nz)/2)):
                            waz = auxfun(nz,iaz,jaz)
                            norm = norm + (wax*way*waz/(2*(lc-iax-iay-iaz-jax-jay-jaz)+1))
    norm = 1/np.sqrt(norm*factor)
    return norm

def aux_elec(nx,ny,nz,ex,norm):
    lc = nx+ny+nz
    factor = np.sqrt(np.pi**3/(2**lc*ex**(lc+3)))
   
    neaux = 0
    if (nx % 2 == 0) and (ny % 2 == 0) and (nz % 2 ==0):
        if nx == 0:
            tx = 1
        else:
            tx = dfac(nx-1)     
        if ny == 0:
            ty = 1
        else:
            ty = dfac(ny-1)
        if nz == 0:
            tz = 1
        else:
            tz = dfac(nz-1)
        neaux = neaux+(tx*ty*tz*norm*factor)
    return neaux

def AHBASE_calcnorms():
    global prim_norm, prim_elec, norm_offset
    norm_offset = np.zeros(nprim)
    prim_norm = np.zeros(ncoeff)
    prim_elec = np.zeros(ncoeff)
    
    norm_offset = np.zeros(nprim, dtype=int)
    for n in range(1,nprim):
        norm_offset[n] = norm_offset[n-1] + prim_order[n-1]

    jcoeff = -1
    for i in range(nprim):
        order = prim_order[i]
        ex = prim_expon[i]
        norm = 0
        neaux = 0
        if order == 1:
            nx, ny, nz = 0, 0, 0
            norm = normalize_aux(nx,ny,nz,ex)
            neaux = aux_elec(nx,ny,nz,ex,norm)
            jcoeff = jcoeff + 1
            prim_norm[jcoeff] = norm
            prim_elec[jcoeff] = neaux
        
        elif order == 3:
            nx, ny, nz = 1, 0, 0
            norm = normalize_aux(nx,ny,nz,ex)
            neaux = aux_elec(nx,ny,nz,ex,norm) 
            for k in range(3):
                jcoeff = jcoeff + 1
                prim_norm[jcoeff] = norm
                prim_elec[jcoeff] = neaux
            
        elif order == 4:
            # s normalization 
            nx, ny, nz = 0, 0, 0
            norm = normalize_aux(nx,ny,nz,ex)
            neaux = aux_elec(nx,ny,nz,ex,norm)
            jcoeff = jcoeff + 1
            prim_norm[jcoeff] = norm
            prim_elec[jcoeff] = neaux
        
            # p normalization
            nx, ny, nz = 1, 0, 0
            norm = normalize_aux(nx,ny,nz,ex)
            neaux = aux_elec(nx,ny,nz,ex,norm)
            for k in range(3):
                jcoeff = jcoeff + 1
                prim_norm[jcoeff] = norm
                prim_elec[jcoeff] = neaux
            
        elif order == 6:
            # d normalization 
            nx, ny, nz = 2, 0, 0
            norm = normalize_aux(nx,ny,nz,ex)
            neaux = aux_elec(nx,ny,nz,ex,norm)
            for k in range(3):
                jcoeff = jcoeff + 1
                prim_norm[jcoeff] = norm
                prim_elec[jcoeff] = neaux
       
            nx, ny, nz = 1, 1, 0
            norm = normalize_aux(nx,ny,nz,ex)
            neaux = aux_elec(nx,ny,nz,ex,norm)
            for k in range(3):
                jcoeff = jcoeff + 1
                prim_norm[jcoeff] = norm
                prim_elec[jcoeff] = neaux
        
        elif order == 10:
            # s normalization
            nx, ny, nz = 0, 0, 0
            norm = normalize_aux(nx,ny,nz,ex)
            neaux = aux_elec(nx,ny,nz,ex,norm)
            jcoeff = jcoeff + 1
            prim_norm[jcoeff] = norm
            prim_elec[jcoeff] = neaux
        
            # p normalization
            nx, ny, nz = 1, 0, 0
            norm = normalize_aux(nx,ny,nz,ex)
            neaux = aux_elec(nx,ny,nz,ex,norm)
            for k in range(3):
                jcoeff = jcoeff + 1
                prim_norm[jcoeff] = norm
                prim_elec[jcoeff] = neaux
        
            # d normalization
            nx, ny, nz = 2, 0, 0
            norm = normalize_aux(nx,ny,nz,ex)
            neaux = aux_elec(nx,ny,nz,ex,norm)
            for k in range(3):
                jcoeff = jcoeff + 1
                prim_norm[jcoeff] = norm
                prim_elec[jcoeff] = neaux
        
            nx, ny, nz = 1, 1, 0
            norm = normalize_aux(nx,ny,nz,ex)
            neaux = aux_elec(nx,ny,nz,ex,norm)
            for k in range(3):
                jcoeff = jcoeff + 1
                prim_norm[jcoeff] = norm
                prim_elec[jcoeff] = neaux

        elif order == 20:
            # s normalization
            nx, ny, nz = 0, 0, 0
            norm = normalize_aux(nx,ny,nz,ex)
            neaux = aux_elec(nx,ny,nz,ex,norm)
            jcoeff = jcoeff + 1
            prim_norm[jcoeff] = norm
            prim_elec[jcoeff] = neaux

            # p normalization
            nx, ny, nz = 1, 0, 0
            norm = normalize_aux(nx,ny,nz,ex)
            neaux = aux_elec(nx,ny,nz,ex,norm)
            for k in range(3):
                jcoeff = jcoeff + 1
                prim_norm[jcoeff] = norm
                prim_elec[jcoeff] = neaux
        
            # d normalization
            nx, ny, nz = 2, 0, 0
            norm = normalize_aux(nx,ny,nz,ex)
            neaux = aux_elec(nx,ny,nz,ex,norm)
            for k in range(3):
                jcoeff = jcoeff + 1
                prim_norm[jcoeff] = norm
                prim_elec[jcoeff] = neaux
        
            nx, ny, nz = 1, 1, 0
            norm = normalize_aux(nx,ny,nz,ex)
            neaux = aux_elec(nx,ny,nz,ex,norm)
            for k in range(3):
                jcoeff = jcoeff + 1
                prim_norm[jcoeff] = norm
                prim_elec[jcoeff] = neaux
            
            # f^3 normalization
            nx, ny, nz = 3, 0, 0
            norm = normalize_aux(nx,ny,nz,ex)
            neaux = aux_elec(nx,ny,nz,ex,norm)
            for k in range(3):
                jcoeff = jcoeff + 1
                prim_norm[jcoeff] = norm
                prim_elec[jcoeff] = neaux

            # f^2,1 normalization
            nx, ny, nz = 2, 1, 0
            norm = normalize_aux(nx,ny,nz,ex)
            neaux = aux_elec(nx,ny,nz,ex,norm)
            for k in range(6):
                jcoeff = jcoeff + 1
                prim_norm[jcoeff] = norm
                prim_elec[jcoeff] = neaux
        
            # f^1,1,1 normalization
            nx, ny, nz = 1, 1, 1
            norm = normalize_aux(nx,ny,nz,ex)
            neaux = aux_elec(nx,ny,nz,ex,norm)
            for k in range(6):
                jcoeff = jcoeff + 1
                prim_norm[jcoeff] = norm
                prim_elec[jcoeff] = neaux
        else:
            print('AH_loadbasis: can only handle up to F auxis, exiting')
    return 

def AHTYPE_load_site_info(par,aux):
    global nsites, ncoeff, off_prim, site_num_prims, basis_index, off_basis_prims
    global site_prim_order, site_prim_expon, coeff_offset
    nsites = 0
    sitetype = []
    sitename = []

    par.seek(0)
    for name in par:
        if 'sitetype' in name:
            x = name.split()
            sitename.append(x[1])
            sitetype.append(x[2])
            nsites = nsites+1

    auxtype = []
    num_prims = []
    par.seek(0)
    for name in par:
        counter = 0
        if 'auxbasis' in name:
            auxname = name[9:]
            y = auxname.split()
            auxtype.append(y[2])
            aux.seek(0)
            for basis in aux:
                if auxname in basis:
                    nshell = int(aux.readline())
                    for i in range(nshell):
                        oprimshell = aux.readline().split()
                        nprimshell = int(oprimshell[0])
                        ordprim = int(oprimshell[1])
                        for j in range(nprimshell):
                            aux.readline()
                            counter = counter+1
            num_prims.append(counter)

    basis_index = []
    for i in sitename:
        counter = 0
        for j in auxtype:
            i = i[:1]
            counter = counter+1
            if i == j:
                basis_index.append(counter)

    site_nprim = 0
    site_num_prims = np.zeros(nsites, dtype=int)
    for n in range(nsites):
        site_num_prims[n] = num_prims[basis_index[n]-1]
        site_nprim = site_nprim + num_prims[basis_index[n]-1]
    
    off_prim = np.zeros(nsites, dtype=int)
    for n in range(1,nsites):
        off_prim[n] = off_prim[n-1] + site_num_prims[n-1]

    nbasis = len(num_prims)
    off_basis_prims = np.zeros(nbasis, dtype=int)
    for n in range(1,nbasis):
        off_basis_prims[n] = off_basis_prims[n-1] + num_prims[n-1]

    site_prim_order = np.zeros(site_nprim, dtype=int)
    site_prim_expon = np.zeros(site_nprim)
    for n in range(nsites):
        off = off_prim[n]
        num = site_num_prims[n]
        ind_basis = basis_index[n]
        offb = off_basis_prims[ind_basis-1]
        for k in range(num):
            site_prim_order[off+k] = prim_order[offb+k]
            site_prim_expon[off+k] = prim_expon[offb+k]
    ncoeff = 0
    for n in range(site_nprim):
        ncoeff = ncoeff + site_prim_order[n]

    coeff_offset = np.zeros(site_nprim,dtype=int)
    for n in range(1,site_nprim):
        coeff_offset[n] = coeff_offset[n-1] + site_prim_order[n-1]
    return

def AHCOEFF_copy_norms_from_basis():
    global cart_coeff_norms, aux_elec
    cart_coeff_norms = np.zeros(ncoeff)
    aux_elec = np.zeros(ncoeff)
    for n in range(nsites):
        off = off_prim[n]
        num = site_num_prims[n]
        ind_basis = basis_index[n]
        offb = off_basis_prims[ind_basis-1]
        for k in range(num):
            order = site_prim_order[off+k]
            coff = coeff_offset[off+k]
            bcoff = norm_offset[offb+k]
            for j in range(order):
                cart_coeff_norms[coff+j] = prim_norm[bcoff+j]
                aux_elec[coff+j] = prim_elec[bcoff+j]
    return

def AHFRAME_get_molframe(par,flabel,coord):
    global names, numbers, molframe, num_deflist, frame_deflist, jsite, local_crds
    num_deflist = 0

    numbers = ["","",""]
    names = ["","",""]
    jsite = []

    par.seek(0)
    Lines = par.readlines()

    for line in Lines:
        if "molframe" in line:
            sline = line.split()
            for n in range(nsites):
                if flabel[n] == sline[1]:
                    molframe = n
        if 'frame' in line and not "molframe" in line:
            sline = line.split()
            num_pt1 = int(sline[2])
            num_pt2 = int(sline[3])
            num_deflist = num_deflist + num_pt1 + num_pt2 
            for n in range(nsites):
                tmp1 = []
                tmp2 = []
                if flabel[n] == sline[1]:
                    tmp1.append(int(sline[2]))
                    tmp1.append(int(sline[3]))
                    tmp2.append(sline[4])
                    for i in range(int(sline[2])):
                        tmp2.append(sline[5+i])
                    numbers[n] = tmp1
                    names[n] = tmp2

    names=names[0]+names[1]+names[2]
    jsite = []
    for label in names:
        for n in range(nsites):
            if flabel[n] == label:
                jsite.append(n+1)
    frame_deflist = np.zeros(5*num_deflist,dtype=int)
   
    site_crds = np.zeros(3*nsites)
    for j in range(nsites):
        for k in range(3):
            site_crds[3*j+k] = coord[j][k]
    
    p1 = np.zeros(3)
    j=0 # aqui debe ir un ciclo
    for k in range(3):
        p1[k] = site_crds[3*j+k]
    
    local_crds = np.zeros((3*nsites))
    for j in range(nsites):
        for k in range(3):
            pass
            local_crds[3*j+k] = site_crds[3*j+k] - p1[k]
    return

def AHFRAME_load_deflist():
    global jsite
    # first pass get num frame def list
    num_deflist = 0
    for n in range(nsites):
        num_pt1 = numbers[n][0]
        num_pt2 = numbers[n][1]
        if n == molframe:
            tmp_p2 = np.zeros(num_pt1, dtype=int)
            tmp_p3 = np.zeros(num_pt2, dtype=int)
            tot_def = num_pt1
            tot_def2 = num_pt2

        for k in range(num_pt1):
            frame_deflist[5*(num_deflist)+0] = n+1
            frame_deflist[5*(num_deflist)+1] = jsite[num_deflist]
            frame_deflist[5*(num_deflist)+2] = n+1
            frame_deflist[5*(num_deflist)+3] = 1
            frame_deflist[5*(num_deflist)+4] = num_pt1
            if n == molframe:
                tmp_p2[k] = jsite[num_deflist]
            num_deflist = num_deflist +1 
        for k in range(num_pt2):
            frame_deflist[5*(num_deflist)+0] = n+1
            frame_deflist[5*(num_deflist)+1] = jsite[num_deflist]
            frame_deflist[5*(num_deflist)+2] = n+1
            frame_deflist[5*(num_deflist)+3] = 2
            frame_deflist[5*(num_deflist)+4] = num_pt2
            if n == molframe:
                tmp_p3[k] = jsite[num_deflist]
            num_deflist = num_deflist +1

    # now calculate the vectors that define molecular frame
    tot_def = 2
    tot_def2 = 1
    tiny = 1e-9
    
    bohr = 0.529177249

    ua = np.zeros((tot_def,3))
    v1 = np.zeros(3)
    for j in range(tot_def):
        jsite = int(tmp_p2[j])-1
        for k in range(3):
            ua[j,k] = local_crds[3*jsite+k]
        wt = np.sqrt(ua[j,0]*ua[j,0]+ua[j,1]*ua[j,1]+ua[j,2]*ua[j,2])
        
        if wt > tiny:
            ua[j,:] = ua[j,:]/wt
        else:
            print("less than tiny")
            ua[j,:] = 0
        for k in range(3):
            v1[k] = v1[k] + ua[j,k]
    wt = np.sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2])
    if wt > tiny:
        v1[:] = v1[:]/wt
    else:
        v1[:] = 0
    
    ub = np.zeros((tot_def2,3))
    u2 = np.zeros(3)
    for j in range(tot_def2):
        jsite = int(tmp_p3[j])-1
        for k in range(3):
            ub[j,k] = local_crds[3*jsite+k]
        wt = np.sqrt(ub[j,0]*ub[j,0]+ub[j,1]*ub[j,1]+ub[j,2]*ub[j,2])
        if wt > tiny:
            ub[j,:] = ub[j,:]/wt
        else:
            ub[j,:] = 0
        for k in range(3):
            u2[k] = u2[k] + ub[j,k]
    
    v2 = np.zeros(3)
    vdu = np.dot(v1,u2)
    for j in range(3):
        v2[j] = u2[j] - vdu*v1[j]
    wt = np.sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2])
    if wt > tiny:
        v2[:] = v2[:]/wt
    else:
        v2[:] = 0
    v3 = np.cross(v1,v2)

    vec = np.zeros(3)
    newcoord=np.zeros((nsites,3))
    for j in range(nsites):
        for k in range(3):
                vec[k] = local_crds[3*j+k]
        r1 = np.dot(vec,v1)
        r2 = np.dot(vec,v2)
        r3 = np.dot(vec,v3)
        local_crds[3*j+k] = (v1[k]*r1+v2[k]*r2+v3[k]*r3)
    
    for k in range(3*nsites):
        local_crds[k] = local_crds[k]*bohr
    
    return

def AHFRAME_build_frames(coord):
    global frame
    nframes = 3  # cambiar despues
    epsilon = 1e-9
    tiny = 1e-9

    # clear the frame def pts
    frame_def_pts = np.zeros(3*2*nsites)
    frame_axis = np.zeros(3)

    frame_axis[0] = 3   #z then x then y
    frame_axis[1] = 1   #z then x then y
    frame_axis[2] = 2   #z then x then y

    # first the def points (ATMSITE_build_frame_def_pts)

    frdefpt = np.zeros((3,2,nframes))
    framedeflist = np.reshape(frame_deflist,(5,7),order='F')

    for n in range(num_deflist):
        i = framedeflist[0,n]-1   # center i gives tail of vector
        j = framedeflist[1,n]-1   # center j gives head of vector
        k = framedeflist[2,n]-1   # the frame number this vector helps define
        l = framedeflist[3,n]-1   # the frame defpoint (1 or 2) within frame
        m = framedeflist[4,n]   # m is number of such (unit) vectors 
  
        if i >= 0 and j >= 0:
            dx = coord[j,0] - coord[i,0]
            dy = coord[j,1] - coord[i,1]
            dz = coord[j,2] - coord[i,2]
            wt = m*np.sqrt(dx*dx+dy*dy+dz*dz)
        
            if abs(wt) < epsilon:
                frdefpt[0,l,k] = frdefpt[0,l,k]
                frdefpt[1,l,k] = frdefpt[1,l,k]
                frdefpt[2,l,k] = frdefpt[2,l,k]
            else:
                frdefpt[0,l,k] = frdefpt[0,l,k] + dx/wt
                frdefpt[1,l,k] = frdefpt[1,l,k] + dy/wt
                frdefpt[2,l,k] = frdefpt[2,l,k] + dz/wt     

    for i in range(3):
        for j in range(2):
            for k in range(nframes):
                if abs(frdefpt[i,j,k]) < tiny: frdefpt[i,j,k] = 0

    # next complete the frames (ATMSITE_defpoints_to_frames)

    frame = np.zeros((3,3,nframes))

    k1 = frame_axis[0]
    k2 = frame_axis[1]
    k3 = frame_axis[2]

    u = np.zeros(3)
    v = np.zeros(3)
    w = np.zeros(3)

    for n in range(nframes):
        # u is unit vector in direction of primary def pt
        for i in range(3):
            u[i] = frdefpt[i,0,n]
        siz = np.linalg.norm(u)

        for i in range(3):
            if abs(siz) < epsilon:
                u[i] = 0
            else:
                u[i] = u[i]/siz
            
        # v is unit vector given by component of secondary pt orthog to u
        for i in range(3):
            v[i] = frdefpt[i,1,n]
        dot = np.dot(u,v)
    
        for i in range(3):
            v[i] = v[i]-dot*u[i]
        siz = np.linalg.norm(v)
    
        for i in range(3):
            if abs(siz) < epsilon:
                v[i] = 0
            else:
                v[i] = v[i]/siz
    
        # w is u cross v
        w = np.cross(u,v)    

        # now define frame
        k1 = 2
        k2 = 0
        k3 = 1

        for i in range(3):
            if abs(u[i]) < epsilon:
                u[i] = 0
            if abs(v[i]) < epsilon:
                v[i] = 0
            if abs(w[i]) < epsilon:
                w[i] = 0
    
        # now define frame
        for i in range(3):
            frame[i,k1,n] = u[i]
            frame[i,k2,n] = v[i]
            frame[i,k3,n] = w[i]

    return

def XFORM_MPOLE_matrix(A_xy,Mpole_xy,order):
    qind1 = [1,2,3,1,1,2]
    qind2 = [1,2,3,2,3,3]
    oind1 = [1,2,3,1,1,2,2,3,3,1]
    oind2 = [1,2,3,1,1,2,2,3,3,2]
    oind3 = [1,2,3,2,3,1,3,1,2,3]
    hind1 = [1,2,3,1,1,2,2,3,3,1,1,2,1,2,3]
    hind2 = [1,2,3,1,1,2,2,3,3,1,1,2,1,2,3]
    hind3 = [1,2,3,1,1,2,2,3,3,2,3,3,2,1,1]
    hind4 = [1,2,3,2,3,1,3,1,2,2,3,3,3,3,2]
    
    # CHARGE case
    Mpole_xy[0,0] = 1
    if order == 1: return
    for j in range(1,order):
        Mpole_xy[1,j] = 0
        Mpole_xy[j,1] = 0
    
    # DIPOLES
    # D'_j
    for j in range(1,4):
        for k in range(1,4):
            Mpole_xy[j,k] = A_xy[j-1,k-1]
            
    if order == 4: return
    for j in range(1,4):
        for k in range(4,order):
            Mpole_xy[j,k] = 0
            Mpole_xy[k,j] = 0
    
    # QUADRUPOLES
    # Q'_kk
    for ind1 in range(3):
        k = qind1[ind1]-1
        for ind2 in range(6):
            i = qind1[ind2]-1
            j = qind2[ind2]-1
            Mpole_xy[ind1+4,ind2+4] = A_xy[k,i]*A_xy[k,j]
    # Q'_kl
    for ind1 in range(3,6):
        k = qind1[ind1]-1
        l = qind2[ind1]-1
        for ind2 in range(6):
            i = qind1[ind2]-1
            j = qind2[ind2]-1
            Mpole_xy[ind1+4,ind2+4] = A_xy[k,i]*A_xy[l,j] + A_xy[k,j]*A_xy[l,i]
    if order == 10: return
    for j in range(4,10):
        for k in range(10,order):
            Mpole_xy[k,j] = 0
            Mpole_xy[j,k] = 0
    
    # OCTUPOLES
    # O'_lll
    for ind1 in range(3):
        l = oind1(ind1)-1
        for ind2 in range(10):
            i = oind1[ind2]-1
            j = oind2[ind2]-1
            k = oind3[ind2]-1
            Mpole_xy[ind1+10,ind2+10] = A_xy[l,i]*A_xy[l,j]*A_xy[l,k]
    # O'_llm
    for ind1 in range(3,9):
        l = oind1[ind1]-1
        m = oind3[ind1]-1
        for ind2 in range(9):
            i = oind1[ind2]-1
            k = oind3[ind2]-1
            Mpole_xy[ind1+10,ind2+10] = (A_xy(l,i)*A_xy(l,i)*A_xy(m,k)
                                         + 2*A_xy(l,i)*A_xy(m,i)*A_xy(l,k))
        Mpole_xy[ind1+10,20] = (A_xy[l,1]*A_xy[l,2]*A_xy[m,3]
                                + A_xy[l,1]*A_xy[m,2]*A_xy[l,3]
                                + A_xy[m,1]*A_xy[l,2]*A_xy[l,3])
    # O'_123
    Mpole_xy[20,11] = 6*A_xy[1,1]*A_xy[2,1]*A_xy[3,1]
    Mpole_xy[20,12] = 6*A_xy[1,2]*A_xy[2,2]*A_xy[3,2]
    Mpole_xy[20,13] = 6*A_xy[1,3]*A_xy[2,3]*A_xy[3,3]  
    for ind2 in range(3,9):
        i = oind1(ind2)-1
        k = oind3(ind2)-1
        Mpole_xy[20,10+ind2] = 2*(A_xy[1,i]*A_xy[2,i]*A_xy[3,k]
                                  + A_xy[1,i]*A_xy[2,k]*A_xy[3,i]
                                  + A_xy[1,k]*A_xy[2,i]*A_xy[3,i])
    Mpole_xy[20,20] = (A_xy[1,1]*A_xy[2,2]*A_xy[3,3]
                       + A_xy[1,1]*A_xy[3,2]*A_xy[2,3]
                       + A_xy[2,1]*A_xy[1,2]*A_xy[3,3]
                       + A_xy[2,1]*A_xy[3,2]*A_xy[1,3]
                       + A_xy[3,1]*A_xy[1,2]*A_xy[2,3]
                       + A_xy[3,1]*A_xy[2,2]*A_xy[1,3])
    if order == 20: return
    for j in range(10,20):
        for k in range(11,20):
            Mpole_xy[k,j] = 0
            Mpole_xy[j,k] = 0
    
    # HEXADECAPOLES
    # H'_mmmm
    for ind1 in range(3):
        m = hind1(ind1)
        for ind2 in range(15):
            i = hind1[ind2]-1
            j = hind2[ind2]-1
            k = hind3[ind2]-1
            l = hind4[ind2]-1
            Mpole_xy[ind1+20,ind2+20] = A_xy[m,i]*A_xy[m,j]*A_xy[m,k]*A_xy[m,l]
    # H'_mmmn
    for ind1 in range(3,9):
        m = hind1[ind1]-1
        n = hind4[ind1]-1
        # can put iiii & iiij together
        for ind2 in range(9):
            i = hind1[ind2]-1
            j = hind4[ind2]-1
            Mpole_xy[ind1+20,ind2+20] = (A_xy[m,i]*A_xy[m,i]*A_xy[m,i]*A_xy[n,j]
                                         + 3*A_xy[m,i]*A_xy[m,i]*A_xy[n,i]*A_xy[m,j])
        # can put iijj & iijk together
        for ind2 in range(9,15):
            i = hind1[ind2]-1
            j = hind3[ind2]-1
            k = hind4[ind2]-1
            Mpole_xy[ind1+20,ind2+20] = (A_xy[m,i]*A_xy[m,i]*A_xy[m,j]*A_xy[n,k]
                                         + A_xy[m,i]*A_xy[m,i]*A_xy[m,k]*A_xy[n,j]
                                         + 2*A_xy[m,i]*A_xy[m,j]*A_xy[m,k]*A_xy[n,i])
    # H'_mmnn
    for ind1 in range(9,12):
        m = hind1[ind1]-1
        n = hind3[ind1]-1
        # can put iiii & iiij together
        for ind2 in range(9):
            i = hind1[ind2]-1
            j = hind4[ind2]-1
            Mpole_xy[ind1+20,ind2+20] = 3*(A_xy[m,i]*A_xy[m,i]*A_xy[n,i]*A_xy[n,j]
                                           + A_xy[m,i]*A_xy[m,j]*A_xy[n,i]*A_xy[n,i])
        # can put iijj & iijk together
        for ind2 in range(9,15):
            i = hind1[ind2]-1
            j = hind3[ind2]-1
            k = hind4[ind2]-1
            Mpole_xy[ind1+20,ind2+20] = (A_xy[m,i]*A_xy[m,i]*A_xy[n,j]*A_xy[n,k]
                                         + A_xy[m,j]*A_xy[m,k]*A_xy[n,i]*A_xy[n,i]
                                         + 2*A_xy[m,i]*A_xy[m,j]*A_xy[n,i]*A_xy[n,k]
                                         + 2*A_xy[m,i]*A_xy[m,k]*A_xy[n,i]*A_xy[n,j])
    # H'_mmnp
    for ind1 in range(12,15):
        m = hind1[ind1]-1
        n = hind3[ind1]-1
        p = hind4[ind1]-1
        # can put iiii & iiij together
        for ind2 in range(9):
            i = hind1[ind2]-1
            j = hind4[ind2]-1
            Mpole_xy[ind1+20,ind2+20] = (3*A_xy[m,i]*A_xy[m,i]*A_xy[n,i]*A_xy[p,j]
                                         + 3*A_xy[m,i]*A_xy[m,i]*A_xy[n,j]*A_xy[p,i]
                                         + 6*A_xy[m,i]*A_xy[m,j]*A_xy[n,i]*A_xy[p,i])
        # can put iijj & iijk together
        for ind2 in range(9,15):
            i = hind1[ind2]-1
            j = hind3[ind2]-1
            k = hind4[ind2]-1
            Mpole_xy[ind1+20,ind2+20] = (A_xy[m,i]*A_xy[m,i]*A_xy[n,j]*A_xy[p,k]
                                         + A_xy[m,i]*A_xy[m,i]*A_xy[n,k]*A_xy[p,j]
                                         + 2 *(A_xy[m,i]*A_xy[m,j]*A_xy[n,i]*A_xy[p,k]
                                               + A_xy[m,i]*A_xy[m,j]*A_xy[n,k]*A_xy[p,i]
                                               + A_xy[m,i]*A_xy[m,k]*A_xy[n,i]*A_xy[p,j]
                                               + A_xy[m,i]*A_xy[m,k]*A_xy[n,j]*A_xy[p,i]
                                               + A_xy[m,j]*A_xy[m,k]*A_xy[n,i]*A_xy[p,i]))
    return

def XFORM_MPOLE(Mpole_xy,Mpole_in,Mpole_out,order):
    if order == 0: return
    Mpole_out[0] = Mpole_xy[0,0]*Mpole_in[0]
    if order == 1: return

    # DIPOLES
    if order == 3:
        for i in range(3):
            Mpole[i] = 0
            for j in range(3):
                Mpole_out[i] = Mpole_out[i] + Mpole_xy[i+1,j+1]*Mpole_in[j]
        return
    for i in range(1,4):
        Mpole_out[i] = 0
        for j in range(1,4):
            Mpole_out[i] = Mpole_out[i] + Mpole_xy[i,j]*Mpole_in[j]
    if order == 4: return
    
    # QUADRUPOLES
    for i in range(4,10):
        Mpole_out[i] = 0
        for j in range(4,10):
            Mpole_out[i] = Mpole_out[i] + Mpole_xy[i,j]*Mpole_in[j]
    if order == 10: return
    
    # OCTUPOLES
    for i in range(10,20):
        Mpole_out[i] = 0
        for j in range(10,20):
            Mpole_out[i] = Mpole_out[i] + Mpole_xy[i,j]*Mpole_in[j]
    if order == 20: return
        
    # HEXADECAPOLES
    for i in range(20,35):
        Mpole_out[i] = 0
        for j in range(20,35):
            Mpole_out[i] = Mpole_out[i] + Mpole_xy[i,j]*Mpole_in[j]
    if order == 35: return
    
    return


def reordenar_inverso(aux_coefs):
    gem = np.copy(aux_coefs)
    psi4 = np.copy(aux_coefs)

    gem[11] = psi4[9]
    gem[12] = psi4[10]
    gem[9] = psi4[11]
    gem[13] = psi4[12]
    gem[10] = psi4[13]

    gem[21] = psi4[19]
    gem[22] = psi4[20]
    gem[19] = psi4[21]
    gem[23] = psi4[22]
    gem[20] = psi4[23]

    gem[31] = psi4[29]
    gem[32] = psi4[30]
    gem[29] = psi4[31]
    gem[33] = psi4[32]
    gem[30] = psi4[33]

    gem[41] = psi4[39]
    gem[42] = psi4[40]
    gem[39] = psi4[41]
    gem[43] = psi4[42]
    gem[40] = psi4[43]

    gem[54] = psi4[52]
    gem[55] = psi4[53]
    gem[52] = psi4[54]
    gem[56] = psi4[55]
    gem[53] = psi4[56]

    gem[67] = psi4[65]
    gem[68] = psi4[66]
    gem[65] = psi4[67]
    gem[69] = psi4[68]
    gem[66] = psi4[69]

    return gem

def cart_to_gherm(aux_coefs):

    coff = 0
    cart_coefs = np.zeros(ncoeff)
    global_herm_coeffs = np.zeros(ncoeff)
 
    for i in range(ncoeff):
        cart_coefs[i] = aux_coefs[i] * cart_coeff_norms[i]

    for i, order in enumerate(site_prim_order):
        alpha = site_prim_expon[i]
        factor = (np.pi/alpha)**(3/2)
        if order == 1:
            global_herm_coeffs[coff] = cart_coefs[coff] * factor
        elif order == 10:
            sum_terms = (cart_coefs[coff+4] + cart_coefs[coff+5] + cart_coefs[coff+6])/(2*alpha)
            global_herm_coeffs[coff] = cart_coefs[coff] + sum_terms
            for k in range(2,5):
                global_herm_coeffs[coff+k-1] = cart_coefs[coff+k-1] / (2*alpha)
            for k in range(5,11):
                global_herm_coeffs[coff+k-1] = cart_coefs[coff+k-1] / ((2*alpha)**2)
            for k in range(1,11):
                global_herm_coeffs[coff+k-1] = global_herm_coeffs[coff+k-1] * factor

        coff = coff + order

    return global_herm_coeffs

def AHFRAME_global_to_local(global_herm_coeffs):
    local_herm_coeffs = np.zeros(ncoeff)
    fr_valid = True

    counter = -1
    for n in range(3):
        A = frame[:,:,n]
        At = np.transpose(A)
        num = site_num_prims[n]
        for k in range(num):
            counter = counter + 1
            order = site_prim_order[counter]
            scoff = coeff_offset[counter]
            Mpole_xy = np.zeros((20,20))
            XFORM_MPOLE_matrix(At, Mpole_xy, order)
            if order == 1:
                local_herm_coeffs[scoff] = global_herm_coeffs[scoff]
            elif order > 1 and not fr_valid:
                local_herm_coeffs[scoff] = global_herm_coeffs[scoff]
            else:
                XFORM_MPOLE(Mpole_xy,global_herm_coeffs[scoff:scoff+order],
                        local_herm_coeffs[scoff:scoff+order],order)

    for valor in local_herm_coeffs:
       print(valor)

    return local_herm_coeffs

