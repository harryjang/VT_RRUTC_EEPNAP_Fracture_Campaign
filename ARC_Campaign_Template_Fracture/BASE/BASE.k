*KEYWORD
$========$=========$=========$=========$=========$=========$=========$=========$
*TITLE
BASE
$========$=========$=========$=========$=========$=========$=========$=========$
*PARAMETER
$#   prmrn      valn     
Rtheta1   THETA1
Rtheta2   THETA2
Rtheta3   THETA3
Rxtrans   0.0
Rytrans   0.0
$========$=========$=========$=========$=========$=========$=========$=========$
*INCLUDE
../control_parameters.k
../target_geometry.k
$========$=========$=========$=========$=========$=========$=========$=========$
*DEFINE_TRANSFORMATION
$   TRANID
         1
$   OPTION        A1        A2        A3        A4        A5        A6        A7
$               COSX      COSY      COSZ      XLOC      YLOC      ZLOC     ANGLE
    ROTATE         1         0         0       0.0       0.0   0.05574   &theta1
    ROTATE         0         1         0       0.0       0.0   0.05574   &theta2
    ROTATE         0         0         1       0.0       0.0   0.05574   &theta3
$             XTRANS    YTRANS    ZTRANS
    TRANSL   &xtrans   &ytrans  Z_SHIFT
$--------$---------$---------$---------$---------$---------$---------$---------$
*INCLUDE_TRANSFORM
$ Filename
../projectile_geometry.k
$   IDNOFF    IDEOFF    IDPOFF    IDMOFF    IDSOFF    IDFOFF    IDDOFF
         0         0         0         0         0         0         0
$   IDROFF
         0
$   FCTMAS    FCTTIM    FCTLEN    FCTTEM   INCOUT1
         0         0         0         0         0
$   TRANID
         1
$========$=========$=========$=========$=========$=========$=========$=========$
*END