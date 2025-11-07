*KEYWORD
*CONTROL_CONTACT
$#  slsfac    rwpnal    islchk    shlthk    penopt    thkchg     orien    enmass
       0.0       0.0         0         0         0         0         0         1
$#  usrstr    usrfrc     nsbcs    interm     xpene     ssthk      ecdt   tiedprj
         0         0         1         0       4.0         0         0         0
$#   sfric     dfric       edc       vfc        th     th_sf    pen_sf      
       0.0       0.0       0.0       0.0       0.0       0.0       0.0
$#  ignore    frceng   skiprwg    outseg   spotstp   spotdel   spothin       
         2         0         0         0         0         0       0.0
$#    isym    nserod    rwgaps    rwgdth     rwksf      icov    swradf    ithoff
         0         0         1       0.0       1.0         0       1.0         0
$#  shledg    pstiff    ithcnt    tdcnof     ftall    unused    shltrw    igactc
         0         0         0         0         0                 0.0         0
*CONTROL_ENERGY
$#    hgen      rwen    slnten     rylen     irgen     
         2         2         1         1         2
*CONTROL_MPP_DECOMPOSITION_DISTRIBUTE_SPH_ELEMENTS
*CONTROL_MPP_DECOMPOSITION_OUTDECOMP
$#   itype     
         1
*CONTROL_MPP_IO_NODUMP
*CONTROL_SPH
$#    ncbs     boxid        dt      idim   nmneigh      form     start      maxv
         1         01.00000E20         3       150         0       0.01.00000E15
$#    cont     deriv       ini     ishow     ierod     icont     iavis     isymp
         0         0         0         1         1         0         0       100
$#    ithk     istab        ql         -   sphsort       
         1         0      0.01                   0
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas     nosol     
5.01000E-3         0       0.5       0.01.000000E8         0
*CONTROL_TIMESTEP
$#  dtinit    tssfac      isdo    tslimt     dt2ms      lctm     erode     ms1st
       0.0       0.9         0       0.0       0.0         0         1         0
$#  dt2msf   dt2mslc     imscl    unused    unused     rmscl    unused      ihdo
       0.0         0         0                           0.0                   0
*DATABASE_GLSTAT
$#      dt    binary      lcur     ioopt     
1.00000E-7         0         0         1
*DATABASE_MATSUM
$#      dt    binary      lcur     ioopt     
1.00000E-7         0         0         1
*DATABASE_NODOUT
$#      dt    binary      lcur     ioopt   option1   option2       
2.00000E-4         0         0         1       0.0         0
*DATABASE_RCFORC
$#      dt    binary      lcur     ioopt     
1.00000E-7         0         0         1
*DATABASE_SLEOUT
$#      dt    binary      lcur     ioopt     
1.00000E-7         0         0         1
*DATABASE_BINARY_D3PLOT
$#      dt      lcdt      beam     npltc    psetid      
0.20000E-3         0         0         0         1
$#   ioopt      rate    cutoff    window      type      pset    
         0       0.0       0.0       0.0         0         0
*DATABASE_HISTORY_NODE_SET
$#     id1       id2       id3       id4       id5       id6       id7       id8
         3         0         0         0         0         0         0         0
*BOUNDARY_SPC_SET
$#    nsid       cid      dofx      dofy      dofz     dofrx     dofry     dofrz
         1         0         1         1         1         1         1         1
*CONTACT_ERODING_NODES_TO_SURFACE_MPP_ID
$#     cid                                                                 title
         1SphereToTarget
$#  ignore     bcket    lcbckt    ns2trk   inititr    parmax    unused    cparm8
         0         1         0         3         2    1.0005                   0
$#    ssid      msid     sstyp     mstyp    sboxid    mboxid       spr       mpr
         2         2         4         3         0         0         0         0
$#      fs        fd        dc        vc       vdc    penchk        bt        dt
      0.18      0.18       0.0       0.0      10.0         0       0.01.00000E20
$#     sfs       sfm       sst       mst      sfst      sfmt       fsf       vsf
       1.0       1.0       0.0       0.0       1.0       1.0       1.0       1.0
$#    isym    erosop      iadj    
         0         1         1
*CONTACT_ERODING_SINGLE_SURFACE_MPP_ID
$#     cid                                                                 title
         2SphereSelfContact
$#  ignore     bcket    lcbckt    ns2trk   inititr    parmax    unused    cparm8
         0         1         0         3         2    1.0005                   0
$#    ssid      msid     sstyp     mstyp    sboxid    mboxid       spr       mpr
         2         0         3         0         0         0         0         0
$#      fs        fd        dc        vc       vdc    penchk        bt        dt
      0.18      0.18       0.0       0.0      10.0         0       0.01.00000E20
$#     sfs       sfm       sst       mst      sfst      sfmt       fsf       vsf
       1.0       1.0       0.0       0.0       1.0       1.0       1.0       1.0
$#    isym    erosop      iadj    
         0         1         1
$#    soft    sofscl    lcidab    maxpar     sbopt     depth     bsort    frcfrq
         2       0.1         0     1.025       2.0        35         1         1
*SET_PART_LIST_TITLE
SubstratePart
$#     sid       da1       da2       da3       da4    solver      
         1       0.0       0.0       0.0       0.0MECH
$#    pid1      pid2      pid3      pid4      pid5      pid6      pid7      pid8
         1         0         0         0         0         0         0         0
*PART
$#                                                                         title
Ti64_Target
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         1         4         4         0         0         0         0         0
*SECTION_SPH_TITLE
SPH_Section
$#   secid      cslh      hmin      hmax    sphini     death     start   sphkern
         4       1.2       0.2       2.0       0.01.00000E20       0.0         0
*MAT_ISOTROPIC_ELASTIC_PLASTIC_TITLE
ElasticPerfPlasticTarget
$#     mid        ro         g      sigy      etan      bulk    
         44.43000E-6      44.0      0.88       0.0     120.0
*PART
$#                                                                         title
Impactor
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         2         3         3         0         0         0         0         0
*SECTION_SOLID_PERI
$#   secid    elform      
         3        48
$       dr     ptype
$#      dr     ptype     
       0.8         0
*MAT_ELASTIC_PERI
$#     mid        ro         e        gt        gs  
         32.44000E-6      72.0   1.00E-5   1.00E20
*INITIAL_VELOCITY_GENERATION
$#nsid/pid      styp     omega        vx        vy        vz     ivatn      icid
         2         2       0.0       0.0       0.0    -160.0         0         0
$#      xc        yc        zc        nx        ny        nz     phase    irigid
       0.0       0.0       0.0       0.0       0.0       0.0         0         0
*DEFINE_CURVE
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         1         0       0.0       0.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0                 0.0
                 1.0                10.0
*END