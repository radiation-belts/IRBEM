!***************************************************************************************************
! This file is part of IRBEM-LIB.
!
!    IRBEM-LIB is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    IRBEM-LIB is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with IRBEM-LIB.  If not, see <http://www.gnu.org/licenses/>.
!
!-----------------------------------------------------------------------------
c
c     drift_loss_cone
c
c     Given only a spacecraft location (plus epoch and field model), return
c     the two pitch angles that partition local pitch angle space into the
c     three usual classes:
c
c        alpha < alpha_blc                 : bounce loss cone
c                                            (mirrors below stop_alt on the
c                                             LOCAL field line, lost within
c                                             a quarter bounce)
c        alpha_blc < alpha < alpha_dlc     : drift loss cone / quasi-trapped
c                                            (survives locally, but mirrors
c                                             below stop_alt somewhere else
c                                             along the drift orbit - on Earth
c                                             this is essentially always the
c                                             South Atlantic Anomaly)
c        alpha > alpha_dlc                 : stably trapped
c
c     Both boundaries are returned as equatorial pitch angles (referenced to
c     Bmin on the local field line) and as local pitch angles at the
c     spacecraft.  Because the drift orbit depends only on Bmirror and the
c     field model, the result is independent of species and energy, and it is
c     symmetric about 90 deg (alpha and 180-alpha share a mirror field, hence
c     a drift shell); the loss cone around 180 deg is 180-alpha_xxx.
c
c     Method
c     ------
c     alpha_blc comes from find_foot_opt: trace the local field line to
c     stop_alt in both hemispheres and take the WEAKER of the two foot point
c     fields, Bm_blc = min(B_north,B_south).  A particle mirrors below
c     stop_alt as soon as Bmirror > Bm_blc.
c
c     alpha_dlc is found by bisection, using trace_drift_bounce_orbit_opt
c     (from drift_bounce_orbit.f) to trace the full drift-bounce orbit at
c     each trial angle and comparing the returned hmin (minimum geodetic
c     altitude anywhere on the drift orbit) against stop_alt.
c
c     The search is anchored at the minimum-B point of the local field line,
c     NOT at the spacecraft, and its parameter is the pitch angle there - the
c     equatorial pitch angle.  This matters: a particle seen at the
c     spacecraft has Bmirror >= B0, so a search over local pitch angle spans
c     only Bmirror in [B0, infinity) and simply cannot reach the more
c     equatorially mirroring half of the distribution.  Off the magnetic
c     equator that is most of it - at 60 deg latitude and 800 km altitude,
c     B0/Bmin is about 59, so no local pitch angle corresponds to an
c     equatorial pitch angle above 7.5 deg.  Anchoring at Bmin makes the
c     whole range Bmirror >= Bmin reachable.
c
c     The bracket is [alpha_blc_eq, 90 deg], valid because hmin increases
c     monotonically with the angle and hmin <= stop_alt at alpha_blc_eq by
c     construction.  At 90 deg the particle mirrors at the minimum-B point
c     itself, the most deeply trapped orbit the field line supports, so a
c     returned alpha_dlc_eq of 90 is a tested result and not a saturation.
c     Bisection stops when the bracket is narrower than dlc_tol (an input,
c     in degrees of equatorial pitch angle) and returns its midpoint, so the
c     result lies on a grid of spacing up to dlc_tol and is within dlc_tol/2
c     of the true boundary.  Neighbouring locations can therefore differ by
c     up to dlc_tol even where the true boundary is smooth.  The local pitch
c     angle can be coarser than that near 90 deg, where
c     d(alpha_loc)/d(alpha_eq) diverges.
c
c     Cost: about 12 full drift-shell traces per point at dlc_tol = 0.1,
c     roughly 1 s at options(3)=0, options(4)=1 on a current CPU.  Each
c     halving of dlc_tol adds one more trace.
c
c     Notes on resolution
c     -------------------
c     hmin is sampled at every one of the Nder = 25*(options(4)+1) longitudes
c     that trace_drift_bounce_orbit_opt visits (check_hmin is called
c     independently of whether the orbit is stored), so options(4) controls
c     how finely the SAA is sampled even though only 25 longitudes are ever
c     stored.  options(4)=0 gives 14.4 deg longitude steps, which can miss
c     the bottom of the SAA; options(4)=1..3 (7.2 to 3.6 deg) is recommended
c     for drift loss cone work.
c
c     Failure handling: the tracer stops closing drift-bounce orbits once the
c     mirror points sink to roughly -100 km altitude, i.e. well inside the
c     loss cone, so a failed trial angle is treated as lost for bracketing
c     purposes.  The trapped end of the bracket (alpha = 90) is always
c     verified first, so a shell that will not close for any pitch angle -
c     an open field line, or magnetopause shadowing at high L - is reported
c     as a failure rather than silently treated as a loss.  See iflag.
c
c--------------------------------------------------------------------------

      SUBROUTINE drift_loss_cone1(kext,options,sysaxes,
     &     iyearsat,idoy,UT,xIN1,xIN2,xIN3,maginput,stop_alt,dlc_tol,
     &     alpha_blc_eq,alpha_dlc_eq,alpha_blc_loc,alpha_dlc_loc,
     &     BLOCAL,BMIN,Lm,Lstar,hmin,hmin_lon,iflag)

c     INPUTS - all with their usual IRBEM meaning, except:
c     REAL*8  stop_alt   - geodetic altitude (km) below which a particle is
c                          considered lost to the atmosphere (typically 100)
c     REAL*8  dlc_tol    - bisection tolerance (deg) on alpha_dlc_eq.  The
c                          result is within dlc_tol/2 of the true boundary.
c                          0.1 is enough for most uses; 0.01 costs about
c                          three more drift-shell traces per point and is
c                          worth it when comparing the boundary between
c                          nearby locations.  Must lie in (0, 90)
c
c     OUTPUTS
c     REAL*8  alpha_blc_eq  - bounce loss cone, equatorial pitch angle (deg)
c     REAL*8  alpha_dlc_eq  - drift  loss cone, equatorial pitch angle (deg)
c     REAL*8  alpha_blc_loc - bounce loss cone, local pitch angle (deg)
c     REAL*8  alpha_dlc_loc - drift  loss cone, local pitch angle (deg).
c                             Either comes back as 90 when the corresponding
c                             mirror field is below B0: such particles mirror
c                             above the spacecraft and never reach it, so
c                             every pitch angle observable there is inside
c                             the cone.  Read the equatorial angle for the
c                             boundary itself in that case
c     REAL*8  BLOCAL        - |B| at the spacecraft (nT)
c     REAL*8  BMIN          - |B| at the minimum-B point of the local field
c                             line (nT); the reference for the equatorial
c                             pitch angles
c     REAL*8  Lm, Lstar     - McIlwain L and L* of the marginally trapped
c                             drift shell, i.e. of the alpha_dlc boundary
c     REAL*8  hmin          - lowest geodetic altitude (km) anywhere on that
c                             marginally trapped drift orbit.  It should come
c                             back just above stop_alt; how far above is a
c                             useful check on convergence
c     REAL*8  hmin_lon      - geodetic longitude (deg) at which hmin occurs.
c                             This is where these particles are lost, and for
c                             Earth it should land in the South Atlantic
c     INTEGER*4 iflag       - status:
c                              0 = normal
c                              1 = the whole field line is in the loss cone
c                                  (alpha_dlc_eq = 90): even the particle
c                                  mirroring at Bmin reaches below stop_alt
c                                  somewhere on its drift orbit
c                              2 = drift loss cone not resolvably wider than
c                                  the bounce loss cone at this location
c                                  (typical when the spacecraft is itself
c                                  over the South Atlantic Anomaly)
c                              3 = alpha_dlc found but BMIN unavailable, so
c                                  only the local angles are valid
c                              4 = boundary found, but one or more trial
c                                  angles inside the loss cone could not be
c                                  traced and were treated as lost.  The
c                                  boundary is still bracketed to within
c                                  dlc_tol; treat as baddata if you want the
c                                  strict policy
c                             -1 = field model / coordinate setup failed
c                             -2 = could not find a foot point at stop_alt in
c                                  either hemisphere; nothing is computed
c                             -3 = the drift-bounce orbit does not close even
c                                  for equatorially mirroring particles
c                                  (open field line, magnetopause shadowing);
c                                  nothing is computed
c                             -4 = stop_alt out of range
c                             -5 = dlc_tol out of range

      IMPLICIT NONE
      INCLUDE 'variables.inc'
c
c     inputs
      INTEGER*4    kext,options(5)
      INTEGER*4    sysaxes
      INTEGER*4    iyearsat
      INTEGER*4    idoy
      REAL*8       UT
      REAL*8       xIN1,xIN2,xIN3
      REAL*8       maginput(25)
      REAL*8       stop_alt,dlc_tol
c
c     outputs
      REAL*8       alpha_blc_eq,alpha_dlc_eq
      REAL*8       alpha_blc_loc,alpha_dlc_loc
      REAL*8       BLOCAL,BMIN
      REAL*8       Lm,Lstar
      REAL*8       hmin,hmin_lon
      INTEGER*4    iflag
c
c     internal
      INTEGER*4    k_ext,k_l,kint
      INTEGER*4    opt(5)
      INTEGER*4    t_resol,r_resol,Ifail
      REAL*8       xGEO(3)
      REAL*8       alti,lati,longi
      REAL*8       R0
c     R0: radius (RE) below which a particle is considered lost and the
c     tracer gives up.  It must sit well BELOW the stop_alt surface, or the
c     drift-bounce orbit stops closing before hmin can register a value
c     smaller than stop_alt.  0.85 RE is about 960 km below sea level: deep
c     enough that low drift shells (L < 1.15, whose mirror points plunge
c     several hundred km "below" the surface in the SAA) still close, and
c     still far outside the core, where the IGRF expansion breaks down.
c     Raising it towards 1.0 makes more points fail with iflag=-3; there is
c     no accuracy penalty for lowering it, only run time.
      PARAMETER   (R0 = 0.85D0)
c
      COMMON /magmod/k_ext,k_l,kint
      INTEGER*4 int_field_select, ext_field_select
c
c     initialize outputs
      alpha_blc_eq  = baddata
      alpha_dlc_eq  = baddata
      alpha_blc_loc = baddata
      alpha_dlc_loc = baddata
      BLOCAL        = baddata
      BMIN          = baddata
      Lm            = baddata
      Lstar         = baddata
      hmin          = baddata
      hmin_lon      = baddata
      iflag         = 0
c
      IF (stop_alt.LT.0.D0 .OR. stop_alt.GE.6378.0D0*500.0D0) THEN
         iflag = -4
         RETURN
      ENDIF
c     no lower bound beyond > 0 is needed: MAXIT caps the bisection, and
c     past about 1D-10 deg the bracket simply stops shrinking
      IF (.NOT.(dlc_tol.GT.0.D0 .AND. dlc_tol.LT.90.D0)) THEN
         iflag = -5
         RETURN
      ENDIF
c
c     local copy of options so we never modify the caller's array
      DO Ifail = 1,5
         opt(Ifail) = options(Ifail)
      ENDDO
      IF (opt(1).EQ.0) opt(1) = 1   ! L* is required by the tracer
      IF (opt(3).LT.0 .OR. opt(3).GT.9) opt(3) = 0
      t_resol = opt(3)+1
      r_resol = opt(4)+1
      k_l     = opt(1)
c
      kint  = int_field_select ( opt(5) )
      k_ext = ext_field_select ( kext )
c
      CALL INITIZE

      CALL init_fields ( kint, iyearsat, idoy, ut, opt(2) )

      CALL get_coordinates ( sysaxes, xIN1, xIN2, xIN3,
     &     alti, lati, longi, xGEO )

      CALL set_magfield_inputs ( k_ext, maginput, ifail )

      IF ( ifail.LT.0 ) THEN
         iflag = -1
         RETURN
      ENDIF

      IF (k_ext .EQ. 13 .OR. k_ext .EQ. 14) THEN
         CALL INIT_TS07D_COEFFS(iyearsat,idoy,ut,ifail)
         CALL INIT_TS07D_TLPR
         IF ( ifail.LT.0 ) THEN
            iflag = -1
            RETURN
         ENDIF
      ENDIF
c
      CALL drift_loss_cone_opt(t_resol,r_resol,xGEO,stop_alt,dlc_tol,R0,
     &     alpha_blc_eq,alpha_dlc_eq,alpha_blc_loc,alpha_dlc_loc,
     &     BLOCAL,BMIN,Lm,Lstar,hmin,hmin_lon,iflag)
c
      END                       ! end subroutine drift_loss_cone1

c     --------------------------------------------------------------------

      SUBROUTINE drift_loss_cone_opt(t_resol,r_resol,xGEO,stop_alt,
     &     dlc_tol,R0,
     &     alpha_blc_eq,alpha_dlc_eq,alpha_blc_loc,alpha_dlc_loc,
     &     BLOCAL,BMIN,Lm,Lstar,hmin,hmin_lon,iflag)
c
c     Core routine.  Assumes the field model has already been initialized
c     (INITIZE / init_fields / set_magfield_inputs) and that xGEO is a GEO
c     cartesian position in RE.  Arguments as in drift_loss_cone1, plus R0.
c     dlc_tol is not range checked here; drift_loss_cone1 does that.
c
      IMPLICIT NONE
      INCLUDE 'variables.inc'
c
c     inputs
      INTEGER*4  t_resol,r_resol
      REAL*8     xGEO(3),stop_alt,dlc_tol,R0
c
c     outputs
      REAL*8     alpha_blc_eq,alpha_dlc_eq
      REAL*8     alpha_blc_loc,alpha_dlc_loc
      REAL*8     BLOCAL,BMIN
      REAL*8     Lm,Lstar
      REAL*8     hmin,hmin_lon
      INTEGER*4  iflag
c
c     internal
      INTEGER*4  Ifail,ilost,iter,nfail
      INTEGER*4  MAXIT
      PARAMETER  (MAXIT = 40)
      REAL*8     Bvec(3),B0
      REAL*8     xeq(3),xanch(3),Banch,Bm_dlc
      REAL*8     alpha_blc_anch,alpha_dlc_anch
      REAL*8     XFOOT(3),BFOOT(3),BFN,BFS,Bm_blc
      REAL*8     alo,ahi,amid,sn
      REAL*8     Lm_t,Lstar_t,hmin_t,hminlon_t
      LOGICAL    Beq_ok
      REAL*8     pi,rad
      COMMON /rconst/rad,pi
c
      alpha_blc_eq  = baddata
      alpha_dlc_eq  = baddata
      alpha_blc_loc = baddata
      alpha_dlc_loc = baddata
      BLOCAL        = baddata
      BMIN          = baddata
      Lm            = baddata
      Lstar         = baddata
      hmin          = baddata
      hmin_lon      = baddata
      iflag         = 0
      nfail         = 0
c
c     ---- local field ----------------------------------------------------
      CALL CHAMP(xGEO,Bvec,B0,Ifail)
      IF (Ifail.LT.0 .OR. B0.LE.0.D0) THEN
         iflag = -1
         RETURN
      ENDIF
      BLOCAL = B0
c
c     ---- minimum B on the local field line --------------------------------
c     The pitch angle search is anchored here, not at the spacecraft.  A
c     particle observed at the spacecraft has Bmirror >= B0, so a search
c     parameterised by the LOCAL pitch angle can never reach mirror fields
c     below B0 and cannot answer the question for the more equatorially
c     mirroring part of the distribution.  Anchoring at the minimum-B point
c     makes the whole range Bmirror >= Bmin reachable, which is the full set
c     of particles bound to this field line.
      CALL loc_equator_opt(xGEO,BMIN,xeq)
      Beq_ok = (BMIN.NE.baddata) .AND. (BMIN.GT.0.D0)
      IF (Beq_ok) THEN
         IF (BMIN.GT.B0) BMIN = B0   ! numerical guard
         Banch    = BMIN
         xanch(1) = xeq(1)
         xanch(2) = xeq(2)
         xanch(3) = xeq(3)
      ELSE
c        fall back to searching from the spacecraft; equatorial pitch angles
c        cannot be reported and the search cannot reach Bmirror < B0
         BMIN     = baddata
         Banch    = B0
         xanch(1) = xGEO(1)
         xanch(2) = xGEO(2)
         xanch(3) = xGEO(3)
      ENDIF
c
c     ---- bounce loss cone from the two foot points at stop_alt -----------
      CALL find_foot_opt(xGEO,stop_alt, 1,XFOOT,BFOOT,BFN)
      CALL find_foot_opt(xGEO,stop_alt,-1,XFOOT,BFOOT,BFS)

      Bm_blc = baddata
      IF (BFN.NE.baddata) Bm_blc = BFN
      IF (BFS.NE.baddata) THEN
         IF (Bm_blc.EQ.baddata) THEN
            Bm_blc = BFS
         ELSE IF (BFS.LT.Bm_blc) THEN
            Bm_blc = BFS         ! weaker foot point sets the loss cone
         ENDIF
      ENDIF
      IF (Bm_blc.EQ.baddata .OR. Bm_blc.LE.0.D0) THEN
         iflag = -2
         RETURN
      ENDIF
c
      CALL dlc_angles(Bm_blc,Banch,B0,BMIN,Beq_ok,
     &     alpha_blc_anch,alpha_blc_loc,alpha_blc_eq)
c
      IF (Bm_blc.LE.Banch) THEN
c        nothing on this field line mirrors above stop_alt at all
         alpha_dlc_anch = 90.D0
         alpha_dlc_loc  = 90.D0
         IF (Beq_ok) alpha_dlc_eq = 90.D0
         iflag = 1
         RETURN
      ENDIF
c
c     ---- drift loss cone -------------------------------------------------
c     Bracket [alpha_blc_anch , 90] in pitch angle AT THE ANCHOR.  hmin
c     increases monotonically with the angle, and hmin <= stop_alt at
c     alpha_blc_anch by construction, so the hmin = stop_alt crossing lies
c     inside the bracket.  At the upper end the particle mirrors at the
c     minimum-B point itself, which is the most deeply trapped orbit this
c     field line supports; if even that is lost, the whole field line is in
c     the drift loss cone and 90 deg is a tested result, not a clamp.
c
c     A trial angle whose drift-bounce orbit fails to close is treated as
c     LOST for bracketing.  That is justified here: the tracer only gives up
c     once mirror points on the drift orbit sink to roughly -100 km, i.e.
c     several hundred km BELOW any sensible stop_alt, so the failure region
c     is strictly inside the loss cone.  It is not justified if the shell
c     opens for some other reason (magnetopause shadowing at high L), which
c     is why the trapped end of the bracket is verified first: if the shell
c     will not close even for equatorially mirroring particles, nothing is
c     computed.  nfail>0 is reported via iflag=4 so that these points can be
c     filtered if desired.
c
      ahi = 90.D0
      alo = alpha_blc_anch
c
      CALL drift_loss_cone_test(t_resol,r_resol,xanch,ahi,R0,stop_alt,
     &     ilost,Lm_t,Lstar_t,hmin_t,hminlon_t)
      IF (ilost.EQ.-1) THEN
         iflag = -3
         RETURN
      ENDIF
      Lm       = Lm_t
      Lstar    = Lstar_t
      hmin     = hmin_t
      hmin_lon = hminlon_t
      IF (ilost.EQ.1) THEN
c        even the particle mirroring at Bmin dips below stop_alt somewhere
         alpha_dlc_anch = 90.D0
         alpha_dlc_loc  = 90.D0
         IF (Beq_ok) alpha_dlc_eq = 90.D0
         iflag = 1
         RETURN
      ENDIF
c
c     is there any drift loss cone at all beyond the bounce loss cone?
      CALL drift_loss_cone_test(t_resol,r_resol,xanch,alo,R0,stop_alt,
     &     ilost,Lm_t,Lstar_t,hmin_t,hminlon_t)
      IF (ilost.EQ.-1) nfail = nfail+1
      IF (ilost.EQ.0) THEN
         alpha_dlc_anch = alpha_blc_anch
         alpha_dlc_loc  = alpha_blc_loc
         alpha_dlc_eq   = alpha_blc_eq
         Lm       = Lm_t
         Lstar    = Lstar_t
         hmin     = hmin_t
         hmin_lon = hminlon_t
         iflag    = 2
         RETURN
      ENDIF
c
      DO iter = 1,MAXIT
         IF ((ahi-alo).LE.dlc_tol) GOTO 20
         amid = 0.5D0*(alo+ahi)
         CALL drift_loss_cone_test(t_resol,r_resol,xanch,amid,R0,
     &        stop_alt,ilost,Lm_t,Lstar_t,hmin_t,hminlon_t)
         IF (ilost.EQ.0) THEN
c           still trapped: this is the new marginally trapped shell
            ahi      = amid
            Lm       = Lm_t
            Lstar    = Lstar_t
            hmin     = hmin_t
            hmin_lon = hminlon_t
         ELSE
            alo = amid
            IF (ilost.EQ.-1) nfail = nfail+1
         ENDIF
      ENDDO
 20   CONTINUE
c
      alpha_dlc_anch = 0.5D0*(alo+ahi)
      sn = SIN(alpha_dlc_anch*rad)
      Bm_dlc = Banch/(sn*sn)
      CALL dlc_angles(Bm_dlc,Banch,B0,BMIN,Beq_ok,
     &     alpha_dlc_anch,alpha_dlc_loc,alpha_dlc_eq)
      IF (.NOT.Beq_ok) iflag = 3
      IF (nfail.GT.0 .AND. iflag.EQ.0) iflag = 4
c
      END                       ! end subroutine drift_loss_cone_opt

c     --------------------------------------------------------------------

      SUBROUTINE dlc_angles(Bm,Banch,B0,BMIN,Beq_ok,
     &     alpha_anch,alpha_loc,alpha_eq)
c
c     Express a mirror field Bm as a pitch angle at the anchor point, at the
c     spacecraft, and at the minimum-B point.
c
c     alpha_loc comes back as 90 deg when Bm < B0.  That is not a clamp: such
c     a particle mirrors above the spacecraft and never reaches it, so every
c     pitch angle that CAN be observed at the spacecraft lies inside the cone.
c
      IMPLICIT NONE
      INCLUDE 'variables.inc'
      REAL*8     Bm,Banch,B0,BMIN
      LOGICAL    Beq_ok
      REAL*8     alpha_anch,alpha_loc,alpha_eq
      REAL*8     pi,rad
      COMMON /rconst/rad,pi
c
      alpha_anch = ASIN(MIN(1.D0,SQRT(Banch/Bm)))/rad
      IF (Bm.LE.B0) THEN
         alpha_loc = 90.D0
      ELSE
         alpha_loc = ASIN(MIN(1.D0,SQRT(B0/Bm)))/rad
      ENDIF
      IF (Beq_ok) THEN
         alpha_eq = ASIN(MIN(1.D0,SQRT(BMIN/Bm)))/rad
      ELSE
         alpha_eq = baddata
      ENDIF
c
      END                       ! end subroutine dlc_angles


c     --------------------------------------------------------------------

      SUBROUTINE drift_loss_cone_test(t_resol,r_resol,xGEO,alpha,R0,
     &     stop_alt,ilost,Lm,Lstar,hmin,hmin_lon)
c
c     Trace the full drift-bounce orbit of a particle with local pitch angle
c     alpha (deg) at xGEO and decide whether it reaches below stop_alt.
c
c     ilost =  1  lost   (hmin < stop_alt somewhere on the drift orbit)
c           =  0  trapped (hmin >= stop_alt all the way round)
c           = -1  the drift-bounce orbit could not be closed; caller decides
c
      IMPLICIT NONE
      INCLUDE 'variables.inc'
c
      INTEGER*4  t_resol,r_resol
      REAL*8     xGEO(3),alpha,R0,stop_alt
      INTEGER*4  ilost
      REAL*8     Lm,Lstar,hmin,hmin_lon
c
c     internal
      INTEGER*4  Ilflag
      INTEGER*4  ind(25)
      REAL*8     alpha1(25),Bmir1(25),xmir(3,25)
      REAL*8     BL,XJ,Bmin_t,Bmir_t
c     these are large; keep them out of the stack
      REAL*8     BLOC(1000,25),POS(3,1000,25)
      SAVE       BLOC,POS
c
      COMMON /flag_L/Ilflag
c
      ilost    = -1
      Lm       = baddata
      Lstar    = baddata
      hmin     = baddata
      hmin_lon = baddata
c
      alpha1(1) = alpha
      CALL find_bm_nalpha(xGEO,1,alpha1,BL,Bmir1,xmir)
      IF (Bmir1(1).EQ.baddata) RETURN
c     note: find_bm_nalpha returns Bmir=0 for alpha=90 and xmir=xGEO;
c     trace_drift_bounce_orbit_opt recomputes Bmir at xmir, so this is fine.
c
      Ilflag = 0                ! never reuse the previous drift shell guess
      CALL trace_drift_bounce_orbit_opt(t_resol,r_resol,xmir,R0,
     &     Lm,Lstar,XJ,BLOC,Bmin_t,Bmir_t,POS,ind,hmin,hmin_lon)
c
      IF (Lstar.EQ.baddata .OR. hmin.EQ.baddata) THEN
         Lm       = baddata
         Lstar    = baddata
         hmin     = baddata
         hmin_lon = baddata
         RETURN
      ENDIF
c
      IF (hmin.LT.stop_alt) THEN
         ilost = 1
      ELSE
         ilost = 0
      ENDIF
c
      END                       ! end subroutine drift_loss_cone_test

c     --------------------------------------------------------------------

      SUBROUTINE drift_loss_cone_multi(ntime,kext,options,sysaxes,
     &     iyearsat,idoy,UT,xIN1,xIN2,xIN3,maginput,stop_alt,dlc_tol,
     &     alpha_blc_eq,alpha_dlc_eq,alpha_blc_loc,alpha_dlc_loc,
     &     BLOCAL,BMIN,Lm,Lstar,hmin,hmin_lon,iflag)
c
c     Loop drift_loss_cone1 over ntime positions.  Array arguments follow the
c     usual IRBEM convention (fixed size NTIME_MAX, maginput(25,NTIME_MAX)).
c
      IMPLICIT NONE
      INCLUDE 'variables.inc'
      INCLUDE 'ntime_max.inc'
c
      INTEGER*4  ntime,kext,options(5),sysaxes
      INTEGER*4  iyearsat(ntime_max),idoy(ntime_max)
      REAL*8     UT(ntime_max)
      REAL*8     xIN1(ntime_max),xIN2(ntime_max),xIN3(ntime_max)
      REAL*8     maginput(25,ntime_max)
      REAL*8     stop_alt,dlc_tol
c
      REAL*8     alpha_blc_eq(ntime_max),alpha_dlc_eq(ntime_max)
      REAL*8     alpha_blc_loc(ntime_max),alpha_dlc_loc(ntime_max)
      REAL*8     BLOCAL(ntime_max),BMIN(ntime_max)
      REAL*8     Lm(ntime_max),Lstar(ntime_max)
      REAL*8     hmin(ntime_max),hmin_lon(ntime_max)
      INTEGER*4  iflag(ntime_max)
c
      INTEGER*4  isat
c
      DO isat = 1,ntime
         IF (xIN1(isat).EQ.baddata .AND. xIN2(isat).EQ.baddata
     &        .AND. xIN3(isat).EQ.baddata) THEN
            alpha_blc_eq(isat)  = baddata
            alpha_dlc_eq(isat)  = baddata
            alpha_blc_loc(isat) = baddata
            alpha_dlc_loc(isat) = baddata
            BLOCAL(isat)        = baddata
            BMIN(isat)          = baddata
            Lm(isat)            = baddata
            Lstar(isat)         = baddata
            hmin(isat)          = baddata
            hmin_lon(isat)      = baddata
            iflag(isat)         = -1
         ELSE
            CALL drift_loss_cone1(kext,options,sysaxes,
     &           iyearsat(isat),idoy(isat),UT(isat),
     &           xIN1(isat),xIN2(isat),xIN3(isat),
     &           maginput(1,isat),stop_alt,dlc_tol,
     &           alpha_blc_eq(isat),alpha_dlc_eq(isat),
     &           alpha_blc_loc(isat),alpha_dlc_loc(isat),
     &           BLOCAL(isat),BMIN(isat),Lm(isat),Lstar(isat),
     &           hmin(isat),hmin_lon(isat),iflag(isat))
         ENDIF
      ENDDO
c
      END                       ! end subroutine drift_loss_cone_multi
